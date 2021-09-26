/*
flactest.cpp
*/

#include <vector>
#include <iostream>
#include <cmath>
#include "flacutil.hpp"

int main()
{
    std::vector<Flac::sample_t> samples;
    for (size_t i = 0; i < 64; i++) {
        samples.push_back(std::int16_t(sin(i / 10.0) * 32760));
    }
    Flac::FlacEncodeOptions options;
    options.bitsPerSample = 16;
    options.bitsPerCoefficient = 12;
    options.minPred = 1;
    options.maxPred = 32;
    options.minPart = 0;
    options.maxPart = 0;
    options.flags = Flac::riceMethodExact | Flac::lpcMethodBruteForce;
    Flac::FlacSubframe subframe(samples, options);
    Flac::FlacSubframeParams params = subframe.params;
    std::cout << params.type << std::endl;
    return 0;
}