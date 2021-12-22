/*
flac.cpp
*/

#include <iostream>
#include <sstream>
#include "flacutil.hpp"

#include "bitutil.hpp"

namespace Flac {
    
    void Flac::processBuffer()
    {
        // if ((options.bitsPerSample & 7) != 0) {
            // std::stringstream sstr;
            // BitBuffer::BitBufferOut bbo(sstr);
            // for (auto it = buffer.begin(); it != buffer.end(); it++) {
                // bbo.write(*it, options.bitsPerSample);
            // }
            // bbo.flush();
            // std::string str = sstr.str();
            // md5.consume(str.data(), str.size());
        // }
        // else {
            size_t bytes = (options.bitsPerSample + 7) >> 3;
            std::vector<std::uint8_t> dbuf;
            for (auto it = buffer.begin(); it != buffer.end(); it++) {
                sample_t sample = *it;
                for (size_t i = 0; i < bytes; i += 1) {
                    dbuf.push_back(sample & 0xff);
                    sample >>= 8;
                }
            }
            md5 << dbuf;
        // }
        frames.push(FlacFrame(buffer, options, numFrames++));
        buffer.clear();
    }
    
    void Flac::writeHeaderTo(std::ostream& stream)
    {
        stream.write("fLaC", 4);
        streaminfoPos = stream.tellp();
        stream.put(0x80); /* STREAMINFO */
        BitBuffer::BitBufferOut bbo(stream);
        constexpr size_t headerLength = 2 + 2 + 3 + 3 + 8 + 16;
        bbo.write(headerLength, 24);
        bbo.write(options.blockSize, 16);
        bbo.write(options.blockSize, 16);
        bbo.write((finalized && frames.empty()) ? minFrame : 0, 24);
        bbo.write((finalized && frames.empty()) ? maxFrame : 0, 24);
        bbo.write(options.sampleRate, 20);
        bbo.write(options.numChannels - 1, 3);
        bbo.write(options.bitsPerSample - 1, 5);
        bbo.write(finalized ? (numSamples >> 32) : 0, 4);
        bbo.write(finalized ? (numSamples & 0xFFFFFFFF) : 0, 32);
        auto digest = md5.finalize();
        for (auto it = digest.begin(); it != digest.end(); it++) {
            bbo.write(*it, 8);
        }
        bbo.flush();
    }
    
    void Flac::rewriteParams(std::ostream& stream) const
    {
        if (streaminfoPos == -1) {
            return;
        }
        std::streampos store = stream.tellp();
        std::streampos target = streaminfoPos;
        target += 8;
        stream.flush();
        stream.seekp(target);
        
        BitBuffer::BitBufferOut bbo(stream);
        bbo.write((finalized && frames.empty()) ? minFrame : 0, 24);
        bbo.write((finalized && frames.empty()) ? maxFrame : 0, 24);
        bbo.write(options.sampleRate, 20);
        bbo.write(options.numChannels - 1, 3);
        bbo.write(options.bitsPerSample - 1, 5);
        bbo.write(finalized ? (numSamples >> 32) : 0, 4);
        bbo.write(finalized ? (numSamples & 0xFFFFFFFF) : 0, 32);
        bbo.flush();
        
        stream.flush();
        stream.seekp(store);
    }
    
    std::ostream& operator<<(std::ostream& stream, Flac &flac)
    {
        FlacFrame frame = flac.frames.front();
        flac.frames.pop();
        std::streampos base = stream.tellp();
        frame.writeTo(stream);
        int frameSize = stream.tellp() - base;
        flac.minFrame = std::min(frameSize, flac.minFrame);
        flac.maxFrame = std::max(frameSize, flac.maxFrame);
        return stream;
    }
    
}
