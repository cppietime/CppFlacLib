/*
flac.cpp
*/

#include <iostream>
#include <sstream>
#include <endian.h>
#include "flacutil.hpp"

#include "bitutil.hpp"

namespace Flac {
    
    void Flac::processBuffer()
    {
        std::cout << "Processing " << buffer.size() << std::endl;
        if ((options.bitsPerSample & 7) != 0) {
            std::stringstream sstr;
            BitBuffer::BitBufferOut bbo(sstr);
            for (auto it = buffer.begin(); it != buffer.end(); it++) {
                bbo.write(*it, options.bitsPerSample);
            }
            bbo.flush();
            std::string str = sstr.str();
            md5.consume(str.data(), str.size());
        }
        else {
            std::vector<std::uint8_t> dbuf;
            for (auto it = buffer.begin(); it != buffer.end(); it++) {
                sample_t sample = *it;
                for (size_t i = 0; i < options.bitsPerSample; i += 8) {
                    dbuf.push_back(sample & 0xff);
                    sample >>= 8;
                }
            }
            md5 << dbuf;
        }
        frames.push(FlacFrame(buffer, options, numFrames++));
        buffer.clear();
    }
    
    void Flac::writeHeaderTo(std::ostream& stream)
    {
        stream.write("fLaC", 4);
        stream.put(0x80); /* STREAMINFO */
        BitBuffer::BitBufferOut bbo(stream);
        constexpr size_t headerLength = 2 + 2 + 3 + 3 + 8 + 16;
        bbo.write(headerLength, 24);
        bbo.write(options.blockSize, 16);
        bbo.write(options.blockSize, 16);
        bbo.write(0, 24);
        bbo.write(0, 24);
        bbo.write(options.sampleRate, 20);
        bbo.write(options.numChannels - 1, 3);
        bbo.write(options.bitsPerSample - 1, 5);
        std::cout << "BPS = " << options.bitsPerSample << std::endl;
        bbo.write(numSamples >> 32, 4);
        bbo.write(numSamples & 0xFFFFFFFF, 32);
        auto digest = md5.finalize();
        for (auto it = digest.begin(); it != digest.end(); it++) {
            bbo.write(*it, 8);
        }
        bbo.flush();
    }
    
    std::ostream& operator<<(std::ostream& stream, Flac &flac)
    {
        FlacFrame frame = flac.frames.front();
        flac.frames.pop();
        frame.writeTo(stream);
        return stream;
    }
    
}
