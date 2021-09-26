/*
flacencode.cpp
*/

#include <cstdlib>
#include <vector>
#include <iostream>
#include "flacutil.hpp"

namespace Flac {
    
    std::vector<std::vector<zip_t>> sumsForEachPartition(
        const std::vector<zip_t>& residue,
        int maxK, int pred, int part)
    {
        std::vector<std::vector<zip_t>> sums;
        size_t nParts = 1 << part;
        size_t n = residue.size() + pred;
        size_t partSize = n >> part;
        for (int k = 0; k <= maxK; k++) {
            std::vector<zip_t> fork;
            size_t start = 0;
            size_t end = partSize - pred;
            for (size_t p = 0; p < nParts; p++) {
                zip_t sum = 0;
                for (size_t i = start; i < end; i++) {
                    if (k == 0) {
                        sum += residue[i];
                    }
                    else {
                        sum += k + 1 + (residue[i] >> k);
                    }
                }
                // std::cout << k << ", " << pred << ", " << part << ": " << sum << '(' << start << ',' << end << ')' << std::endl;
                fork.push_back(sum);
                start = end;
                end += partSize;
            }
            sums.push_back(fork);
        }
        return sums;
    }
    
    void decreaseParititionOfSums(std::vector<std::vector<zip_t>>& sums, size_t maxK)
    {
        size_t originalSize = sums[0].size();
        if (originalSize == 0 || maxK == 0) {
            return;
        }
        for (size_t k = 0; k <= maxK; k++) {
            for (size_t i = 0; (i << 1) < originalSize; i++) {
                sums[k][i] = sums[k][i << 1] + sums[k][(i << 1) + 1];
            }
            for (size_t i = 0; (i << 1) < originalSize; i++) {
                sums[k].pop_back();
            }
        }
    }
    
    std::pair<size_t, int> idealK(
        const std::vector<std::vector<zip_t>>& sums, size_t n,
        int partition, bool useEstimation)
    {
        size_t bestBits;
        int bestK;
        for (int k = 0; k < sums.size(); k++) {
            size_t bits = useEstimation ? sizeBitsOfSum(sums[0][partition], n, k)
                : sums[k][partition];
            if (bits < bestBits || k == 0) {
                bestBits = bits;
                bestK = k;
            }
        }
        return std::pair<size_t, int>(bestBits, bestK);
    }
    
    std::tuple<size_t, int, std::vector<int>> idealRiceEncoding(
        const std::vector<zip_t>& data,
        int pred, int minPart, int maxPart,
        bool useEstimation)
    {
        size_t bestBits;
        int bestPartition;
        std::vector<int> bestK;
        
        auto sums = sumsForEachPartition(data, useEstimation ? 0 : FLAC_MAX_K, pred, maxPart);
        for (int part = maxPart; part >= minPart; part--) {
            size_t bits = 0;
            size_t nParts = 1 << part;
            std::vector<int> k;
            for (size_t i = 0; i < nParts; i++) {
                auto bitsK = idealK(sums, data.size() - pred, i, useEstimation);
                bits += bitsK.first;
                k.push_back(bitsK.second);
            }
            if (part == maxPart || bits < bestBits) {
                bestBits = bits;
                bestPartition = part;
                bestK = k;
            }
            decreaseParititionOfSums(sums, useEstimation ? 0 : FLAC_MAX_K);
        }
        return std::tie(bestBits, bestPartition, bestK);
    }
    
    FlacSubframeParams::FlacSubframeParams(
                const std::vector<zip_t>& residue,
                const FlacEncodeOptions& options,
                const std::vector<lpc_t>& fCoef,
                int lpcShift) : 
        lpcShift {lpcShift}
    {
            type = LPC;
            coefficients = fCoef;
            predOrder = coefficients.size();
            auto params = idealRiceEncoding(
                residue, predOrder,
                options.minPart, options.maxPart,
                (options.flags & riceMethodMask) == riceMethodEstimate);
            partOrder = std::get<1>(params);
            bitSize =
                ((options.bitsPerSample + options.bitsPerCoefficient) * predOrder + 4 + 5)
                + 2 + 4 + (5 << partOrder)
                + std::get<0>(params);
            kParams = std::get<2>(params);
    }
            
    FlacSubframeParams::FlacSubframeParams(
            const std::vector<zip_t>& residue,
            const FlacEncodeOptions& options,
            int predOrder) :
        predOrder {predOrder}
    {
        type = FIXED;
        auto params = idealRiceEncoding(
            residue, predOrder,
            options.minPart, options.maxPart,
            (options.flags & riceMethodMask) == riceMethodEstimate);
        partOrder = std::get<1>(params);
        bitSize =
            options.bitsPerSample * predOrder
            + 2 + 4 + (5 << partOrder)
            + std::get<0>(params);
        kParams = std::get<2>(params);
    }
    
    FlacSubframe::FlacSubframe(std::vector<sample_t>& data, const FlacEncodeOptions& options)
    {
        bool constant = true;
        for (size_t i = 1; i < data.size(); i++) {
            if (data[i] != data[i - 1]) {
                constant = false;
                break;
            }
        }
        if (constant) {
            std::cout << "Constant subframe\n";
            params.type = CONSTANT;
            sample = data[0];
            return;
        }
        
        FlacSubframeParams bestParams;
        bestParams.type = VERBATIM;
        bestParams.bitSize = data.size() * options.bitsPerSample;
        std::cout << "Verbatim takes " << bestParams.bitSize << " bits\n";
        
        int lpcMethod = options.flags & lpcMethodMask;
        
        if (lpcMethod == lpcMethodFixed) {
            auto fixedResidue = calcResidueFixed(data, options.maxPred);
            for (int pred = options.minPred; pred <= options.maxPred; pred++) {
                FlacSubframeParams fixParams(fixedResidue[pred - 1], options, pred);
                std::cout << "Fix param " << pred << " takes " << fixParams.bitSize << " bits\n";
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
            std::cout << "LPC order " << lpcParams.predOrder << " takes " << lpcParams.bitSize << " bits\n";
            if (lpcParams.bitSize < bestParams.bitSize) {
                bestParams = lpcParams;
            }
        }
        
        params = bestParams;
        
        switch (params.type) {
            case VERBATIM:
                std::cout << "Verbatim subframe\n";
                rawData = data;
                break;
            case FIXED:
                std::cout << "Fixed LPC of order " << params.predOrder << std::endl;
                zippedResidue = calcResidueFixed(data, params.predOrder).back();
                rawData = std::vector<sample_t>(&data[0], &data[params.predOrder]);
                break;
            case LPC:
                zippedResidue = calcResidueLpc(data, params.coefficients, params.lpcShift);
                rawData = std::vector<sample_t>(&data[0], &data[params.predOrder]);
                break;
        }
    }
    
}