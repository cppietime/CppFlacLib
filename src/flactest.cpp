/*
flactest.cpp
*/

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdint>
#include <string>
#include "flacutil.hpp"

#include "bitutil.hpp"

int main()
{
    std::srand(std::time(nullptr));
    std::vector<Flac::sample_t> samples;
    for (size_t i = 0; i < 44100; i++) {
        // samples.push_back(i);
        samples.push_back(1 * (std::int16_t)(32767 * sin(i * 440.0 * 2 * M_PI / 44100.0)));
        // samples.push_back((std::int16_t)(32767 * sin(i * 440.0 * 2 * M_PI / 44100.0)));
    }
    Flac::FlacEncodeOptions options(1, 16, 44100);
    options.bitsPerCoefficient = 12;
    options.minPred = 16;
    options.maxPred = 32;
    options.minPart = 0;
    options.maxPart = 2;
    options.flags = Flac::riceMethodExact | Flac::lpcMethodBruteForce;
    // options.blockSize = 512;
    Flac::Flac flac(options);
    flac << samples;
    flac.finalize();
    
    std::stringstream sstr;
    flac.writeHeaderTo(sstr);
    while (!flac.empty()) {
        sstr << flac;
    }
    std::string str = sstr.str();
    
    std::ofstream out("test.flac", std::ios_base::out | std::ios_base::binary);
    for (size_t i = 0; i < str.size(); i++) {
        // std::cout << std::hex << std::uppercase << std::setw(2) << std::setfill('0') << (int)(unsigned char)str.data()[i] << ", ";
        out.put(str.data()[i]);
    }
    out.close();
    std::cout << std::endl;
    return 0;
}