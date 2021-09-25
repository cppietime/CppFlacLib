/*
flacencode.cpp
*/

#include <cstdlib>
#include <vector>
#include "flacutil.hpp"

void Flac::decreaseParititionOfSums(std::vector<std::vector<sample_t>>& sums, size_t maxK)
{
    size_t originalSize = sums[0].size();
    if (originalSize == 0 || maxK == 0) {
        return;
    }
    for (size_t k = 0; k < maxK; k++) {
        for (size_t i = 0; (i << 1) < originalSize; i++) {
            sums[k][i] = sums[k][i << 1] + sums[k][(i << 1) + 1];
        }
        for (size_t i = 0; (i << 1) < originalSize; i++) {
            sums[k].pop_back();
        }
    }
}