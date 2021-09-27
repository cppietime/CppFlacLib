/*
lpc.cpp
*/

#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "flacutil.hpp"

namespace Flac {
    
    const static float FLAC_LPC_ZERO = 0.001;
    const static int FLAC_MAX_SHIFT = 15;
    
    std::vector<float> autocorrelation(
        const std::vector<float>& data, size_t maxLag, size_t bitsPerSample)
    {
        std::vector<float> ac;
        ac.reserve(maxLag + 1);
        for (size_t i = 0 ; i <= maxLag; i++) {
            float aci = 1;
            for (size_t j = i; j < data.size(); j++) {
                aci += (float)(data[j] * data[j - i]);
            }
            ac.push_back(aci);
        }
        return ac;
    }
    
    std::vector<std::vector<float>> calcLpcCoeffs(
        const std::vector<sample_t>& data, int maxOrder, size_t bitsPerSample)
    {
        maxOrder = std::min(maxOrder, (int)data.size() - 2);
        
        std::vector<float> windowed;
        windowed.reserve(data.size());
        float norm = 1.0 / (1 << (bitsPerSample - 1));
        float center = (data.size() - 1) / 2.0;
        for (size_t i = 0; i < data.size(); i++) {
            float factor = (i - center) / (center + 1);
            factor = 1 - factor * factor;
            windowed.push_back(data[i] * norm * factor);
        }
        
        std::vector<float> autocorr = autocorrelation(windowed, maxOrder, bitsPerSample);
        norm = autocorr[0];
        for (auto it = autocorr.begin(); it != autocorr.end(); it++) {
            // *it /= norm;
            std::cout << *it << ", ";
        }
        std::cout << std::endl;
        
        std::vector<std::vector<float>> coeffs;
        coeffs.reserve(maxOrder);
        float dac0 = (std::abs(autocorr[0]) <= FLAC_LPC_ZERO) ?
            (1.0 / FLAC_LPC_ZERO) : 1.0 / autocorr[0];
        std::vector<float> x(1, autocorr[1] / autocorr[0]);
        x.reserve(maxOrder);
        coeffs.push_back(x); /* 1st order coeff */
        
        std::vector<float> fvec(1, 1 / autocorr[0]); /* 1st order forward/backward vector */
        fvec.reserve(maxOrder);
        std::vector<float> fnext;
        fnext.reserve(maxOrder);
        
        for (size_t order = 1; order < maxOrder; order++) {
            float error = 0;
            for (size_t i = 0; i < order; i++) {
                /* Dot product of backward vector and autocorr */
                error += fvec[i] * autocorr[order - i];
            }
            
            float den = 1 - error * error;
            if (std::abs(den) < FLAC_LPC_ZERO) {
                den = FLAC_LPC_ZERO;
            }
            float alpha = 1 / den;
            float beta = alpha * -error;
            
            /* The next forward vector, or
            [fvec.[0]] * alpha + [0.reversed_fvec] * beta */
            fnext.clear();
            fnext.push_back(alpha * fvec[0]);
            for (size_t i = 1; i < order; i++) {
                fnext.push_back(alpha * fvec[i] + beta * fvec[order - i]);
            }
            fnext.push_back(beta * fvec[0]);
            fvec = fnext;
            
            error = 0;
            for (size_t i = 0; i < order; i++) {
                /* Dot product of previous x and reversed autocorr */
                error += x[i] * autocorr[order - i];
            }
            
            /* x = [previous_x.[0]] + reverse(f) * (y[orer] - error) */
            for (size_t i = 0; i < order; i++) {
                x[i] += fvec[order - i] * (autocorr[order + 1] - error);
            }
            x.push_back(fvec[0] * (autocorr[order + 1] - error));
            coeffs.push_back(x);
            
            /* Sanity check */
            // float sanity = 0;
            // for (size_t j = 0; j <= order; j++) {
                // float att = 0;
                // for (size_t i = 0; i <= order; i++) {
                    // att += (autocorr[std::abs((int)i - (int)j)] * x[i]);
                // }
                // att -= autocorr[j + 1];
                // sanity += std::abs(att);
            // }
            // std::cout << (sanity) << " sanity\n";
            
            for (auto it = x.begin(); it != x.end(); it++) {
                std::cout << *it << ", ";
            }
            std::cout << std::endl;
        }
        return coeffs;
    }
    
    std::pair<std::vector<lpc_t>, int> quantizeLpcCoeffs(
        const std::vector<float>& fCoef, size_t bitsPerCoefficient)
    {
        int bitsNeeded = 0;
        for (auto it = fCoef.begin(); it != fCoef.end(); it++) {
            bitsNeeded = std::max(bitsNeeded, (int)std::ceil(log2(std::abs(*it + 1))));
        }
        bitsNeeded++;
        int shift = bitsPerCoefficient - bitsNeeded;
        shift = std::min(shift, FLAC_MAX_SHIFT);
        // shift = 4;
        float factor = (shift >= 0) ? (1 << shift) : (1.0 / (1 << -shift));
        float error = 0;
        std::vector<lpc_t> quantized;
        quantized.reserve(fCoef.size());
        for (auto it = fCoef.begin(); it != fCoef.end(); it++) {
            float f = *it * factor;
            error += f;
            lpc_t q = std::floor(f);
            error -= q;
            quantized.push_back(q);
            // std::cout << *it << " quantized to " << q << " >> " << shift << std::endl;
        }
        return std::pair<std::vector<lpc_t>, int>(quantized, shift);
    }
    
    std::vector<std::vector<zip_t>> calcResidueFixed(
        const std::vector<sample_t>& data, int maxOrder)
    {
        maxOrder = std::min(maxOrder, (int)data.size() - 1);
        std::vector<std::vector<sample_t>> residues(maxOrder);
        std::vector<std::vector<zip_t>> zipped(maxOrder);
        for (size_t i = 1; i < data.size(); i++) {
            residues[0].push_back(data[i] - data[i - 1]);
            zipped[0].push_back(signZip(residues[0].back()));
            for (size_t j = 2; j <= maxOrder; j++) {
                if (j > i) {
                    break;
                }
                residues[j - 1].push_back(residues[j - 2][i + 1 - j] - residues[j - 2][i - j]);
                zipped[j - 1].push_back(signZip(residues[j - 1].back()));
            }
        }
        return zipped;
    }
    
    std::vector<zip_t> calcResidueLpc(
        const std::vector<sample_t>& data,
        const std::vector<lpc_t>& coeffs, int shift)
    {
        size_t order  = coeffs.size();
        
        std::vector<zip_t> residue;
        residue.reserve(data.size() - order);
        
        for (size_t i = order; i < data.size(); i++) {
            sample_t pred = 0;
            for (size_t j = 0; j < order; j++) {
                pred += coeffs[j] * data[i - 1 - j];
            }
            pred >>= shift;
            sample_t diff = std::round(data[i] - pred);
            zip_t zipped = signZip(diff);
            residue.push_back(zipped);
        }
        return residue;
    }
    
    FlacSubframeParams naiveIdealLpcOrder(
        const std::vector<std::vector<float>>& coeffs,
        const std::vector<sample_t>& data,
        const FlacEncodeOptions& options)
    {
        auto quantized = quantizeLpcCoeffs(coeffs[options.minPred - 1], options.bitsPerCoefficient);
        auto residue = calcResidueLpc(data, quantized.first, quantized.second);
        FlacSubframeParams best(residue, options, quantized.first, quantized.second);
        // std::cout << "Order " << options.minPred << " has " << best.bitSize << " bits\n";
        
        for (int order = options.minPred; order < options.maxPred; order++) {
            quantized = quantizeLpcCoeffs(coeffs[order], options.bitsPerCoefficient);
            residue = calcResidueLpc(data, quantized.first, quantized.second);
            FlacSubframeParams candidate(residue, options, quantized.first, quantized.second);
            // std::cout << "Order " << (order + 1) << " has " << candidate.bitSize << " bits\n";
            if (candidate.bitSize < best.bitSize) {
                best = candidate;
            }
        }
        
        return best;
    }
    
    FlacSubframeParams binaryIdealLpcOrder(
        const std::vector<std::vector<float>>& coeffs,
        const std::vector<sample_t>& data,
        const FlacEncodeOptions& options)
    {
        int level = (options.minPred + options.maxPred + 1) / 2;
        auto quantized = quantizeLpcCoeffs(coeffs[level - 1], options.bitsPerCoefficient);
        auto residue = calcResidueLpc(data, quantized.first, quantized.second);
        FlacSubframeParams best(residue, options, quantized.first, quantized.second);
        
        for (int step = (options.maxPred + 1 - options.minPred) / 2; step > 0; step >>= 1) {
            int down = std::max(level - step, (int)options.minPred);
            int up = std::min((int)options.maxPred, level + step);
            int newLvl = level;
            
            if (down != level) {
                quantized = quantizeLpcCoeffs(coeffs[down - 1], options.bitsPerCoefficient);
                residue = calcResidueLpc(data, quantized.first, quantized.second);
                FlacSubframeParams candidate(residue, options, quantized.first, quantized.second);
                if (candidate.bitSize < best.bitSize) {
                    best = candidate;
                    newLvl = down;
                }
            }
            
            if (up != level) {
                quantized = quantizeLpcCoeffs(coeffs[up - 1], options.bitsPerCoefficient);
                residue = calcResidueLpc(data, quantized.first, quantized.second);
                FlacSubframeParams candidate(residue, options, quantized.first, quantized.second);
                if (candidate.bitSize < best.bitSize) {
                    best = candidate;
                    newLvl = up;
                }
            }
            
            level = newLvl;
        }
        
        return best;
    }
    
    FlacSubframeParams linearIdealLpcOrder(
        const std::vector<std::vector<float>>& coeffs,
        const std::vector<sample_t>& data,
        const FlacEncodeOptions& options,
        int levels)
    {
        int step = (options.maxPred + 1 - options.minPred) / (levels - 1);
        
        auto quantized = quantizeLpcCoeffs(coeffs[options.minPred - 1], options.bitsPerCoefficient);
        auto residue = calcResidueLpc(data, quantized.first, quantized.second);
        FlacSubframeParams best(residue, options, quantized.first, quantized.second);
        
        for (int order = options.minPred + step; order <= options.maxPred; order += step) {
            quantized = quantizeLpcCoeffs(coeffs[order - 1], options.bitsPerCoefficient);
            residue = calcResidueLpc(data, quantized.first, quantized.second);
            FlacSubframeParams candidate(residue, options, quantized.first, quantized.second);
            if (candidate.bitSize < best.bitSize) {
                best = candidate;
            }
        }
        
        return best;
    }
    
}