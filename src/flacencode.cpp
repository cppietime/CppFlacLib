/*
flacencode.cpp
*/

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <endian.h>
#include "flacutil.hpp"

#include "bitutil.hpp"

namespace Flac {
    
    std::vector<std::vector<usum_t>> sumsForEachPartition(
        const std::vector<zip_t>& residue,
        int maxK, int pred, int part)
    {
        std::vector<std::vector<usum_t>> sums;
        sums.reserve(maxK + 1);
        
        size_t nParts = 1 << part;
        size_t n = residue.size() + pred;
        size_t partSize = n >> part;
        
        std::vector<usum_t> fork;
        fork.reserve(nParts);
        
        for (int k = 0; k <= maxK; k++) {
            // std::vector<zip_t> fork;
            fork.clear();
            size_t start = 0;
            size_t end = partSize - pred;
            for (size_t p = 0; p < nParts; p++) {
                usum_t sum = 0;
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
    
    void decreaseParititionOfSums(std::vector<std::vector<usum_t>>& sums, size_t maxK)
    {
        size_t originalSize = sums[0].size();
        if (originalSize == 0) {
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
        const std::vector<std::vector<usum_t>>& sums, size_t n,
        int partition, bool useEstimation)
    {
        if (useEstimation) {
            int k = estimateK(sums[0][partition], n);
            size_t bitSize = sizeBitsOfSum(sums[0][partition], n, k);
            // std::cout << "Estimated k=" << k << " with size " << bitSize << std::endl;
            return std::pair<size_t, int>(bitSize, k);
        }
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
        bool useEstimation, int maxK)
    {
        /* Find highest candidate partition order given data.size() and pred */
        int pmax = BitManip::lsbSet(data.size() + pred);
        while ((data.size() + pred) >> pmax <= pred && pmax > 0) {
            pmax--;
        }
        minPart = std::min(minPart, pmax);
        maxPart = std::min(maxPart, pmax);
        
        size_t bestBits;
        int bestPartition;
        std::vector<int> bestK;
        
        auto sums = sumsForEachPartition(data, useEstimation ? 0 : maxK, pred, maxPart);
        for (int part = maxPart; part >= minPart; part--) {
            size_t nParts = 1 << part;
            size_t bits = 5 * nParts;
            std::vector<int> k;
            k.reserve(nParts);
            for (size_t i = 0; i < nParts; i++) {
                auto bitsK = idealK(sums, data.size(), i, useEstimation);
                bits += bitsK.first;
                k.push_back(bitsK.second);
            }
            if (part == maxPart || bits < bestBits) {
                bestBits = bits;
                bestPartition = part;
                bestK = k;
            }
            decreaseParititionOfSums(sums, useEstimation ? 0 : maxK);
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
                (options.flags & riceMethodMask) == riceMethodEstimate,
                options.maxK);
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
            (options.flags & riceMethodMask) == riceMethodEstimate,
            options.maxK);
        partOrder = std::get<1>(params);
        bitSize =
            options.bitsPerSample * predOrder
            + 2 + 4 + (5 << partOrder)
            + std::get<0>(params);
        kParams = std::get<2>(params);
    }
    
}