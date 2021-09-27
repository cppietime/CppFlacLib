/*
flacsubframe.cpp
*/

#include <iostream>

#include <vector>
#include "flacutil.hpp"

#include "bitutil.hpp"

namespace Flac {
    
    FlacSubframe::FlacSubframe(
        std::vector<sample_t>& data, const FlacEncodeOptions& options) :
            options {options}
    {
        bool constant = true;
        for (size_t i = 1; i < data.size(); i++) {
            if (data[i] != data[i - 1]) {
                constant = false;
                break;
            }
        }
        if (constant) {
            params.type = CONSTANT;
            sample = data[0];
            return;
        }
        
        FlacSubframeParams bestParams;
        bestParams.type = VERBATIM;
        bestParams.bitSize = data.size() * options.bitsPerSample;
        
        int lpcMethod = options.flags & lpcMethodMask;
        
        if (lpcMethod != lpcMethodNone) {
            auto fixedResidue = calcResidueFixed(data, options.maxPred);
            for (int pred = options.minPred; pred <= options.maxPred; pred++) {
                FlacSubframeParams fixParams(fixedResidue[pred - 1], options, pred);
                // std::cout << "Fix param " << pred << " takes " << fixParams.bitSize << " bits\n";
                if (fixParams.bitSize < bestParams.bitSize) {
                    bestParams = fixParams;
                }
            }
        }
        else if (lpcMethod != lpcMethodNone) {
            auto lpcCoeffs = calcLpcCoeffs(data, options.maxPred, options.bitsPerSample);
            FlacSubframeParams lpcParams;
            if (lpcMethod == lpcMethodBinary) {
                lpcParams = binaryIdealLpcOrder(lpcCoeffs, data, options);
            }
            else if (lpcMethod == lpcMethodBruteForce) {
                lpcParams = naiveIdealLpcOrder(lpcCoeffs, data, options);
            }
            else if (lpcMethod != lpcMethodEstimate) {
                int levels = 1;
                switch (lpcMethod) {
                    case lpcMethodLevel2:
                        levels = 2;
                        break;
                    case lpcMethodLevel4:
                        levels = 4;
                        break;
                    case lpcMethodLevel8:
                        levels = 8;
                        break;
                }
                lpcParams = linearIdealLpcOrder(lpcCoeffs, data, options, levels);
            }
            // std::cout << "LPC order " << lpcParams.predOrder << " takes " << lpcParams.bitSize << " bits\n";
            if (lpcParams.bitSize < bestParams.bitSize) {
                bestParams = lpcParams;
            }
        }
        
        params = bestParams;
        
        switch (params.type) {
            case VERBATIM:
                // std::cout << "Verbatim subframe\n";
                rawData = data;
                break;
            case FIXED:
                // std::cout << "Fixed LPC of order " << params.predOrder << std::endl;
                zippedResidue = calcResidueFixed(data, params.predOrder).back();
                rawData = std::vector<sample_t>(&data[0], &data[params.predOrder]);
                break;
            case LPC:
                zippedResidue = calcResidueLpc(data, params.coefficients, params.lpcShift);
                rawData = std::vector<sample_t>(&data[0], &data[params.predOrder]);
                break;
        }
    }
    
    void FlacSubframe::writeResidue(BitBuffer::BitBufferOut& bbo)
    {
        int bigK = 0;
        for (auto it = params.kParams.begin(); it < params.kParams.end(); it++) {
            if (*it > 14) {
                bigK = 1;
            }
        }
        bbo.write(bigK, 2);
        bbo.write(params.partOrder, 4);
        size_t start = 0;
        size_t partSize = (zippedResidue.size() + params.predOrder) >> params.partOrder;
        size_t end = partSize - params.predOrder;
        for (size_t p = 0; p < (1 << params.partOrder); p++) {
            int k = params.kParams[p];
            k = std::max(1, k);
            bbo.write(k, bigK + 4);
            while (start < end) {
                auto code = riceEncode(zippedResidue[start++], k);
                while (code.first > 32) {
                    bbo.write(0, 32);
                    code.first--;
                }
                bbo.write(0, code.first);
                bbo.write(1, 1);
                bbo.write(code.second, k);
            }
            end += partSize;
        }
    }
    
    void FlacSubframe::writeTo(BitBuffer::BitBufferOut& bbo)
    {
        switch (params.type) {
            case CONSTANT:
                bbo.write(0, 8);
                bbo.write(sample, options.bitsPerSample);
                break;
            case VERBATIM:
                bbo.write(1 << 1, 8);
                for (auto it = rawData.begin(); it != rawData.end(); it++) {
                    bbo.write(*it, options.bitsPerSample);
                }
                break;
            case FIXED: {
                bbo.write((8 | params.predOrder) << 1, 8);
                for (auto it = rawData.begin(); it != rawData.end(); it++) {
                    bbo.write(*it, options.bitsPerSample);
                }
                writeResidue(bbo);
                break;
            }
            case LPC:
                bbo.write((32 | (params.predOrder - 1)) << 1, 8);
                for (auto it = rawData.begin(); it != rawData.end(); it++) {
                    bbo.write(*it, options.bitsPerSample);
                }
                bbo.write(options.bitsPerCoefficient - 1, 4);
                bbo.write(params.lpcShift, 5);
                for (auto it = params.coefficients.begin(); it != params.coefficients.end(); it++) {
                    bbo.write(*it, options.bitsPerCoefficient);
                }
                writeResidue(bbo);
                break;
        }
    }
    
}