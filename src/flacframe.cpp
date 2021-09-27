/*
flacframe.cpp
*/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdint>
#include "flacutil.hpp"

namespace Flac {
    
    FlacFrame::FlacFrame(
        const std::vector<sample_t>& data,
        const FlacEncodeOptions& options,
        std::uint32_t frameNo) :
            options {options},
            blockSize {data.size() / options.numChannels},
            frameNo {frameNo}
    {
        /* Break up data into channels and compact for MD5 */
        std::vector<std::vector<sample_t>> channels(options.numChannels);
        for (auto it = channels.begin(); it != channels.end(); it++) {
            it->reserve(data.size() / options.numChannels);
        }
        for (size_t i = 0; i < data.size(); i++) {
            sample_t sample = data[i];
            channels[i % options.numChannels].push_back(sample);
        }
        
        /* Check decorrelation */
        if (options.numChannels == 2) {
            std::vector<sample_t> mid, side;
            mid.reserve(channels[0].size());
            side.reserve(channels[0].size());
            for (size_t i = 0; i < channels[0].size(); i++) {
                side.push_back(channels[0][i] - channels[1][i]);
                mid.push_back((channels[0][i] + channels[1][i]) >> 1);
            }
            std::vector<sample_t> *modes[4] = {&channels[0], &channels[1], &side, &mid};
            size_t bits[4];
            for (int i = 0; i < 4; i++) {
                auto residue = calcResidueFixed(*modes[i], 2)[1];
                zip_t sum = 0;
                for (auto it = residue.begin(); it != residue.end(); it++) {
                    sum += *it;
                }
                bits[i] = estimateK(sum, residue.size());
            }
            size_t score[4] =
                {bits[0] + bits[1], bits[0] + bits[2], bits[1] + bits[2], bits[2] + bits[3]};
            int bestScore = 0;
            for (int i = 1; i < 4; i++) {
                if (score[i] < score[bestScore]) {
                    bestScore = i;
                }
            }
            switch (bestScore) {
                case 0: { /* Left-Right */
                    subframes.push_back(FlacSubframe(channels[0], options));
                    subframes.push_back(FlacSubframe(channels[1], options));
                    channelMode = modeFor(2);
                    break;
                }
                case 1: { /* Left-Side */
                    FlacEncodeOptions sideOptions = options;
                    sideOptions.bitsPerSample++;
                    subframes.push_back(FlacSubframe(channels[0], options));
                    subframes.push_back(FlacSubframe(side, sideOptions));
                    channelMode = LEFT_SIDE;
                    break;
                }
                case 2: { /* Right-Side */
                    FlacEncodeOptions sideOptions = options;
                    sideOptions.bitsPerSample++;
                    subframes.push_back(FlacSubframe(channels[1], options));
                    subframes.push_back(FlacSubframe(side, sideOptions));
                    channelMode = RIGHT_SIDE;
                    break;
                }
                case 3: { /* Mid-Side */
                    FlacEncodeOptions sideOptions = options;
                    sideOptions.bitsPerSample++;
                    subframes.push_back(FlacSubframe(mid, options));
                    subframes.push_back(FlacSubframe(side, sideOptions));
                    channelMode = MID_SIDE;
                    break;
                }
            }
        }
        else {
            channelMode = modeFor(options.numChannels);
            for (auto it = channels.begin(); it != channels.end(); it++) {
                subframes.push_back(FlacSubframe(*it, options));
            }
        }
    }
    
    const static size_t blockSizeCodes[] = {
        SIZE_MAX,
        192,
        576, 576 << 1, 576 << 2, 576 << 3,
        SIZE_MAX, SIZE_MAX,
        0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000
    };
    
    const static size_t rateSizeCodes[] = {
        SIZE_MAX,
        88200,
        176400,
        192000,
        8000,
        16000,
        22050,
        24000,
        32000,
        44100,
        48000,
        96000,
        SIZE_MAX,
        SIZE_MAX,
        SIZE_MAX,
        SIZE_MAX
    };
    
    const static int bitsCodes[] = {
        -1,
        8, 12,
        -1,
        16, 20, 24,
        -1
    };
    
    void FlacFrame::writeTo(std::ostream& stream) const
    {
        std::stringstream sstr;
        BitBuffer::BitBufferOut bbo(sstr);
        bbo.write(0xFFF8, 16); /* Sync code */
        int sizeCode = -1;
        for (int i = 0; i < 16; i++) {
            if (blockSizeCodes[i] == blockSize) {
                sizeCode = i;
                break;
            }
        }
        if (sizeCode == -1) {
            if (blockSize <= 0x100) {
                sizeCode = 6;
            }
            else {
                sizeCode = 7;
            }
        }
        bbo.write(sizeCode, 4); /* Block size */
        int rateCode = -1;
        for (int i = 0; i < 16; i++) {
            if (rateSizeCodes[i] == options.sampleRate) {
                rateCode = i;
                break;
            }
        }
        if (rateCode == -1) {
            if (options.sampleRate < 256) {
                rateCode = 12;
            }
            else if (options.sampleRate < 65538) {
                rateCode = 13;
            }
            else {
                rateCode = 14;
            }
        }
        bbo.write(rateCode, 4); /* Sample rate */
        bbo.write(channelMode, 4); /* Channel assignments */
        int bitsCode = -1;
        for (int i = 0; i < 16; i++) {
            if (bitsCodes[i] == options.bitsPerSample) {
                bitsCode = i;
                break;
            }
        }
        bbo.write(bitsCode << 1, 4); /* Bits per sample */
        bbo.writeUtf8(frameNo); /* Frame number */
        
        /* Encoded block size */
        if (sizeCode == 6) {
            bbo.write(blockSize - 1, 8);
        }
        else if (sizeCode == 7) {
            bbo.write(blockSize - 1, 16);
        }
        
        /* Encoded sample rate */
        if (rateCode == 12) {
            bbo.write(options.sampleRate, 8);
        }
        else if (rateCode == 13) {
            bbo.write(options.sampleRate, 16);
        }
        else if (rateCode == 14) {
            bbo.write(options.sampleRate / 10, 16);
        }
        // bbo.flush();
        
        std::string str = sstr.str();
        std::uint8_t crc8 = Digest::crc8(str.data(), str.size());
        bbo.write(crc8, 8);
        // bbo.flush();
        
        for (auto it = subframes.begin(); it != subframes.end(); it++) {
            it->writeTo(bbo);
        }
        bbo.flush();
        
        str = sstr.str();
        std::uint16_t crc16 = Digest::crc16(str.data(), str.size());
        bbo.write(crc16, 16);
        bbo.flush();
        
        str = sstr.str();
        stream.write(str.data(), str.size());
    }
    
}
