/*
lpc.cpp
*/

#include <vector>
#include <utility>
#include <cmath>
#include <iostream>
#include "flacutil.hpp"

namespace Flac {
    
    std::vector<float> autocorrelation(const std::vector<float>& data, size_t maxLag)
    {
        std::vector<float> ac;
        for (size_t i = 0 ; i <= maxLag; i++) {
            float aci = 0;
            for (size_t j = i; j < data.size(); j++) {
                aci += (float)(data[j] * data[j - i]) / 32767 / 32767;
            }
            ac.push_back(aci);
        }
        return ac;
    }
    
    /* TODO replace with real function. This is just a placeholder to get something */
    std::vector<std::vector<float>> calcLpcCoeffs(
        const std::vector<sample_t>& data, int maxOrder, size_t bitsPerSample)
    {
        std::vector<float> windowed;
        windowed.reserve(data.size());
        float norm = 1.0 / (1 << bitsPerSample);
        float center = (data.size() + 1) / 2.0;
        for (size_t i = 0; i < data.size(); i++) {
            float factor = (i - center) / (center + 0.5);
            factor = 1 - factor * factor;
            windowed.push_back(data[i] * norm * factor);
        }
        std::vector<std::vector<float>> coeffs;
        /* Should this apply a welch window? */
        std::vector<float> autocorr = autocorrelation(windowed, maxOrder);
        for (auto it = autocorr.begin(); it != autocorr.end(); it++) {
            std::cout << *it << ", ";
        }
        std::cout << std::endl;
        std::vector<float> x(1, autocorr[1] / autocorr[0]);
        std::cout << "Starting at " << x[0] << std::endl;
        coeffs.push_back(x); /* 1st order coeff */
        std::vector<float> fvec(1, 1 / autocorr[0]); /* 1st order forward/backward vector */
        for (size_t order = 1; order < maxOrder; order++) {
            float error = 0;
            for (size_t i = 0; i < order; i++) {
                /* Dot product of backward vector and autocorr */
                error += fvec[order - 1 - i] * autocorr[i + 1];
            }
            std::cout << "Error = " << error << std::endl;
            
            float alpha = 1 / (1 - error * error);
            float beta = alpha * -error;
            
            /* The next forward vector, or
            [fvec.[0]] * alpha + [0.reversed_fvec] * beta */
            std::vector<float> fnext;
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
            
            /* x = [previous_x.[0]] + beta * (y - error) */
            for (size_t i = 0; i < order; i++) {
                x[i] += beta * (autocorr[i + 1] - error);
            }
            x.push_back(beta * (autocorr[order + 1] - error));
            coeffs.push_back(x);
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
        float factor = (shift >= 0) ? (1 << shift) : (1.0 / (1 << -shift));
        float error = 0;
        std::vector<lpc_t> quantized;
        for (auto it = fCoef.begin(); it != fCoef.end(); it++) {
            float f = *it * factor;
            error += f;
            lpc_t q = std::floor(error);
            error -= q;
            quantized.push_back(q);
            // std::cout << *it << " quantized to " << q << " >> " << shift << std::endl;
        }
        return std::pair<std::vector<lpc_t>, int>(quantized, shift);
    }
    
    std::vector<std::vector<zip_t>> calcResidueFixed(
        const std::vector<sample_t>& data, int maxOrder)
    {
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
        std::vector<zip_t> residue;
        std::vector<float> fCoef;
        float factor = (shift >= 0) ? (1.0 / (1 << shift)) : (1 << -shift);
        for (auto it = coeffs.begin(); it != coeffs.end(); it++) {
            fCoef.push_back(*it * factor);
        }
        size_t order  = coeffs.size();
        for (size_t i = order; i < data.size(); i++) {
            float pred = 0;
            for (size_t j = 0; j < order; j++) {
                pred += fCoef[j] * data[i - 1 - j];
            }
            sample_t diff = data[i] - pred;
            zip_t zipped = signZip(diff);
            // std::cout << data[i] << "->" << zipped << std::endl;
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
        std::cout << "Order " << options.minPred << " has " << best.bitSize << " bits\n";
        
        for (int order = options.minPred; order < options.maxPred; order++) {
            quantized = quantizeLpcCoeffs(coeffs[order], options.bitsPerCoefficient);
            residue = calcResidueLpc(data, quantized.first, quantized.second);
            FlacSubframeParams candidate(residue, options, quantized.first, quantized.second);
            std::cout << "Order " << (order + 1) << " has " << candidate.bitSize << " bits\n";
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
            int down = std::max(level - step, options.minPred);
            int up = std::min(options.maxPred, level + step);
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