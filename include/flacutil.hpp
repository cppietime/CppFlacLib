/*
flacutil.hpp
Yaakov Schectman, 2021
FLAC encoder
*/

#ifndef _FLACUTIL_HPP
#define _FLACUTIL_HPP

#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <climits>
#include <sstream>
#include <utility>
#include <tuple>
#include <vector>
#include <queue>
#include <iostream>

#include "bitutil.hpp"

// #define FLAC_MAX_LPC 32
// #define FLAC_MIN_LPC 1 // Not sure if I'll need this
// #define FLAC_MAX_K 30
// #define FLAC_MAX_FIXED 4

namespace Flac {
    
    constexpr unsigned int FLAC_MAX_LPC = 32;
    constexpr unsigned int FLAC_MIN_LPC = 1;
    constexpr unsigned int FLAC_MAX_K = 30;
    constexpr unsigned int FLAC_MAX_FIXED = 4;
    constexpr unsigned int FLAC_MAX_PARTORDER = 15;
    constexpr unsigned int FLAC_MAX_SINGLE_HZ = 65537;
    constexpr unsigned int FLAC_DEFAULT_BLOCKSIZE = 8192;
    constexpr unsigned int FLAC_DEFAULT_LPCBITS = 12;
    constexpr unsigned int FLAC_DEFAULT_MINPRED = 1;
    constexpr unsigned int FLAC_DEFAULT_MAXPRED = 32;
    constexpr unsigned int FLAC_DEFAULT_MINPART = 0;
    constexpr unsigned int FLAC_DEFAULT_MAXPART = 14;
    constexpr size_t FLAC_MIN_BPS = 8;
    
    using flag_t = std::uint8_t;
    const flag_t riceMethodMask = 0x01;
    const flag_t riceMethodEstimate = 0x00;
    const flag_t riceMethodExact = 0x01;
    const flag_t lpcMethodMask = 0x0E;
    const flag_t lpcMethodFixed = 0x00;
    const flag_t lpcMethodEstimate = 0x02;
    const flag_t lpcMethodBruteForce = 0x04;
    const flag_t lpcMethodBinary = 0x06;
    const flag_t lpcMethodLevel2 = 0x08;
    const flag_t lpcMethodLevel4 = 0x0A;
    const flag_t lpcMethodLevel8 = 0x0C;
    const flag_t lpcMethodNone = 0x0E;
    
    using lpc_t = std::int16_t;
    using sample_t = std::int32_t;
    using zip_t = std::uint32_t;
    using residue_t = std::uint16_t;
    using usum_t = std::uint32_t;
    
    struct FlacEncodeOptions {
        public:
            size_t blockSize;
            size_t bitsPerSample;
            size_t sampleRate;
            size_t bitsPerCoefficient;
            size_t numChannels;
            unsigned int minPred;
            unsigned int maxPred;
            unsigned int minPart;
            unsigned int maxPart;
            unsigned int maxK;
            flag_t flags;
            
            FlacEncodeOptions(
                size_t numChannels,
                size_t bitsPerSample,
                size_t sampleRate,
                size_t blockSize = FLAC_DEFAULT_BLOCKSIZE,
                size_t bitsPerCoefficient = FLAC_DEFAULT_LPCBITS,
                unsigned int minPred = FLAC_DEFAULT_MINPRED,
                unsigned int maxPred = FLAC_DEFAULT_MAXPRED,
                unsigned int minPart = FLAC_DEFAULT_MINPART,
                unsigned int maxPart = FLAC_DEFAULT_MAXPART,
                unsigned int maxK = FLAC_MAX_K,
                flag_t flags = lpcMethodFixed | riceMethodEstimate
            ) :
                numChannels {numChannels},
                bitsPerSample {std::max(FLAC_MIN_BPS, bitsPerSample)},
                sampleRate {sampleRate},
                blockSize {blockSize},
                bitsPerCoefficient {bitsPerCoefficient},
                minPred {minPred},
                maxPred {maxPred},
                minPart {minPart},
                maxPart {maxPart},
                maxK {maxK},
                flags {flags}
            {
                    if ((flags & lpcMethodMask) == lpcMethodFixed) {
                        maxPred = std::min(maxPred, FLAC_MAX_FIXED);
                    }
                    maxPart = std::min(maxPart, FLAC_MAX_PARTORDER);
                    minPart = std::min(maxPart, minPart);
                    maxK = std::min(maxK, FLAC_MAX_K);
                    minPred = std::min(minPred, FLAC_MAX_LPC);
                    maxPred = std::min(maxPred, FLAC_MAX_LPC);
                    minPart = std::min(minPart, FLAC_MAX_PARTORDER);
                    maxPart = std::min(maxPart, FLAC_MAX_PARTORDER);
                    if (sampleRate > FLAC_MAX_SINGLE_HZ) {
                        sampleRate = 10 * (sampleRate / 10);
                    }
            }
    };
    
    enum SubframeType {
        CONSTANT,
        VERBATIM,
        FIXED,
        LPC
    };
    
    struct FlacSubframeParams {
        public:
            SubframeType type;
            int predOrder;
            int partOrder;
            int lpcShift;
            size_t bitSize;
            size_t bitsPerSample;
            std::vector<int> kParams;
            std::vector<lpc_t> coefficients;
            
            FlacSubframeParams(){}
            
            /*
            Find the partition order, k params, and bit length for an ideal encoding using
            the provided options if LPC or FIXED (idealRiceEncoding),
            Otherwise, only make sure the type is set properly.
            */
            FlacSubframeParams(
                const std::vector<zip_t>& residue,
                const FlacEncodeOptions& options,
                const std::vector<lpc_t>& fCoef,
                int lpcShift = 0);
            
            FlacSubframeParams(
                const std::vector<zip_t>& residue,
                const FlacEncodeOptions& options,
                int predOrder);
    };
    
    class FlacSubframe {
        private:
            const FlacEncodeOptions options;
            FlacSubframeParams params;
            sample_t sample; // For constant frames
            std::vector<sample_t> rawData; // For verbatim frames
            std::vector<zip_t> zippedResidue;
            
            void writeResidue(BitBuffer::BitBufferOut& bbo) const;
        public:
            /*
            Construct from unencoded data with provided encoding options
            
            First, check if it should be encoded as a constant subblock.
            If not, and if using dynamic LPC, calculate LPC coefficients (calcLPC).
            Calculate and store the verbatim bit length (options.bitsPerSample * data.size()).
            If using fixed LPC, calculate and zip the residue for each order (calcResidueFixed),
            and find and store each order's bit length and parameters (idealRiceEncoding).
            If using dynamic LPC and estimating LPC order, use the highest order LPC that uses
            all of its coefficients,
            otherwise, using the lpcMethod, try each candidate LPC order storing their bit lengths
            and parameters (idealRiceEncoding).
            Store the parameters for whichever attempted method yields the least bit length
            */
            FlacSubframe(std::vector<sample_t>& data, const FlacEncodeOptions& options);
            
            inline const FlacSubframeParams& getParams()
            {
                return params;
            }
            
            /*
            Write this subframe, not necessarily byte-aligned, to a BitBufferOut
            */
            void writeTo(BitBuffer::BitBufferOut& bbo) const;
    };
    
    enum ChannelMode {
        LEFT_SIDE = 8,
        RIGHT_SIDE = 9,
        MID_SIDE = 10
    };
    
    /*
    Convert a number of channels to a channel assignment index
    */
    inline int modeFor(int numChannels)
    {
        return numChannels - 1;
    }
    
    class FlacFrame {
        private:
            const FlacEncodeOptions options;
            size_t blockSize;
            std::vector<FlacSubframe> subframes;
            int channelMode;
            std::uint32_t frameNo;
        public:
            FlacFrame(
                const std::vector<sample_t>& data,
                const FlacEncodeOptions& options,
                std::uint32_t frameNo);
            
            /*
            Encode this frame and write out to stream
            */
            void writeTo(std::ostream& stream) const;
    };
    
    class Flac {
        private:
            const FlacEncodeOptions& options;
            std::vector<sample_t> buffer;
            size_t blockSize;
            Digest::MD5Context md5;
            std::queue<FlacFrame> frames;
            size_t numFrames;
            size_t numSamples;
            int minFrame;
            int maxFrame;
            bool finalized;
            std::streampos streaminfoPos;
            
            /*
            Create a temporary LE buffer of buffer's samples for MD5 and process it
            */
            void processBuffer();
        public:
            Flac(const FlacEncodeOptions& options) :
                options {options},
                blockSize {options.blockSize},
                numFrames {0},
                numSamples {0},
                minFrame {INT_MAX},
                maxFrame {0},
                finalized {false},
                streaminfoPos{-1}
            {
                    buffer.reserve(blockSize * options.numChannels);
            }
        
            /*
            Load data into the buffer, processing into frames when the buffer fills up
            */
            template <class T>
            inline Flac& operator<<(std::vector<T> data)
            {
                finalized = false;
                for (auto it = data.begin(); it != data.end(); it++) {
                    buffer.push_back((sample_t)*it);
                    if (buffer.size() == blockSize * options.numChannels) {
                        processBuffer();
                    }
                }
                numSamples += data.size() / options.numChannels;
                return *this;
            }
            
            /*
            If the buffer is still holding anything, process it into a final frame
            */
            inline void finalize()
            {
                if (!buffer.empty()) {
                    processBuffer();
                }
                finalized = true;
            }
            
            /*
            The number of frames stored
            */
            inline size_t storedFrames() const
            {
                return frames.size();
            }
            
            /*
            Returns true iff frames is empty
            */
            inline bool empty() const
            {
                return frames.empty();
            }
            
            /*
            Write the fLaC FourCC and STREAMINFO metadata block to a stream
            */
            void writeHeaderTo(std::ostream& stream);
            
            /*
            If a position from unknown params in a STREAMINFO metadata block is stored,
            rewrite those params to the stream.
            This method assumes the stream passed to it is the same as the stream passed to
            writeHeaderTo
            */
            void rewriteParams(std::ostream& stream) const;
            
            /*
            Pop one frame from the queue and encode it into the provided stream
            */
            friend std::ostream& operator<<(std::ostream& stream, Flac &flac);
            
    };
    
    /* Does NOT write the entire FLAC file. Only pops and writes one frame at a time */
    std::ostream& operator<<(std::ostream& stream, Flac &flac);
    
    /*
    For a given prediction order pred and partition order part, return a nested vector X such that
    X[k][i] is the number of bits it would take to encode the i-th partition of residue with rice
    paramter k,
    unless maxK == 0, in which case
    X[0][i] is the sum of all residues in the i-th partition
    */
    std::vector<std::vector<usum_t>> sumsForEachPartition(
        const std::vector<zip_t>& residue,
        int maxK, int pred, int part);
    
    /*
    Change sums to be for one lower partition order
    */
    void decreasePartitionOfSums(std::vector<std::vector<usum_t>>& sums, int maxK);
    
    /*
    Zip signed samples to unsigned values
    */
    inline zip_t signZip(sample_t sample)
    {
        sample <<= 1;
        if (sample < 0) {
            return ~sample;
        }
        return sample;
    }
    
    /*
    Encode a value with a given rice encoding
    */
    inline std::pair<size_t, residue_t> riceEncode(zip_t sample, int k)
    {
        return std::pair<size_t, residue_t>(sample >> k, sample & ((1 << k) - 1));
    }
    
    /*
    Estimate how many bits it takes to rice-encode values by their sum
    */
    inline size_t sizeBitsOfSum(usum_t sum, size_t n, int k)
    {
        return (k + 1) * n + ((sum - (n >> 1)) >> k);
    }
    
    /*
    Estimate and return the ideal K param
    */
    inline int estimateK(usum_t sum, size_t n)
    {
        float mean = (float)sum / n + 0.5;
        int k = BitManip::msbSet((int)std::ceil(mean)) - 1;
        // std::cout << mean << "->" << k << std::endl;
        return std::max(k, 0);
    }
    
    /*
    Find which k minimizes count(sums[k][partition])
    With useEstimation, count(n) = sizeBitsOfSum
    Otherwise, count(n) = n
    
    Return <bit length, k>
    */
    std::pair<size_t, int> idealK(
        const std::vector<std::vector<usum_t>>& sums, size_t n,
        int partition, bool useEstimation);
    
    /*
    Returns <length in bits, partition order, k paramter for each partition>
    
    First, calculate the sums for each partition using the max order (sumsForEachPartition).
    Then, for each candidate partition order, in descending order, iterate over those sums to find
    the ideal k parameter for that partition, and estimate the number of bits needed to encode it.
    Decrease the partition order on the sums by (decreasePartitionOfSums).
    Return whichever partition order minimizes the resultant bits, as well as its k-vector, as
    <bit length, partition order, list of k params>
    */
    std::tuple<size_t, int, std::vector<int>> idealRiceEncoding(
        const std::vector<zip_t>& data,
        int pred, int minPart, int maxPart,
        bool useEstimation, int maxK);

    /*
    Calculate and return the autocorrelations of lags [0, maxLag]
    */
    std::vector<float> autocorrelation(
        const std::vector<float>& data, size_t maxLag, size_t bitsPerSample);

    /*
    Calculate the LPC coefficients of n samples, stored in data, using Levinson-Durbin recursion
    
    Note that in the returned vector X, X[i] holds the (i + 1)th order coefficients, as order 0
    has no coefficients
    */
    std::vector<std::vector<float>> calcLpcCoeffs(
        const std::vector<sample_t>& data, int maxOrder, size_t bitsPerSample);
    
    /*
    Quantize LPC coefficients with a given bit precision, returning the pair of
    <quantized integer coefficients, bit shift to normalize them>
    */
    std::pair<std::vector<lpc_t>, int> quantizeLpcCoeffs(
        const std::vector<float>& floats, size_t bitsPerCoefficient);
    
    /*
    Calculate the 1st through maxOrder-th order residues of data using fixed LPC coefficients
    
    In the returned vector X, X[i] holds the (i + 1)th order residue, as 0th order is just data
    */
    std::vector<std::vector<zip_t>> calcResidueFixed(
        const std::vector<sample_t>& data, int maxOrder);
    
    /*
    Calculate the LPC residue of data using provided coefficients
    Order can be inferred from coeffs.size()
    */
    std::vector<zip_t> calcResidueLpc(
        const std::vector<sample_t>& data,
        const std::vector<lpc_t>& coeffs, int shift);
    
#define FLAC_LPC_THRESH 0.10

    /*
    Find the greatest order LPC coefficients whose final coefficient has an absolute value of
    at least FLAC_LPC_THRESH
    
    Returns the order in [minOrder, maxOrder], not the index in [minOrder - 1, maxOrder)
    */
    inline int approxIdealLpcOrder(
        const std::vector<std::vector<float>>& coeffs,
        int minOrder, 
        int maxOrder)
    {
        for (size_t i = maxOrder; i >= minOrder; i--) {
            const std::vector<float>& orderCoeffs = coeffs[i - 1];
            if (std::abs(orderCoeffs[i - 1]) >= FLAC_LPC_THRESH) {
                return i;
            }
        }
        return minOrder;
    }
    
    /*
    Tries the specified number of levels linearly distributed over [minOrder, maxOrder] to find
    which yields the least bit size
    
    Start with level = minOrder, step = (maxOrder + 1 - minOrder) / levels,
    and repeat incrementing level by step
    */
    FlacSubframeParams linearIdealLpcOrder(
        const std::vector<std::vector<float>>& coeffs,
        const std::vector<sample_t>& data,
        const FlacEncodeOptions& options,
        int levels);
    
    /*
    Brute-force searches every prediction order in [minOrder, maxOrder] to find
    which yields the least bit size
    */
    FlacSubframeParams naiveIdealLpcOrder(
        const std::vector<std::vector<float>>& coeffs,
        const std::vector<sample_t>& data,
        const FlacEncodeOptions& options);
    
    /*
    Tries a binary search over [minOrder, maxOrder] to find
    which yields the least bit size
    
    Starts at level = floor((minOrder + maxOrder) / 2), step = ceil((maxOrder + 1 - minOrder) / 2)
    Finds ideal value out of level, level - step, and level + step, clamping to the bounds
    level is set to the preferred out of those three, and step is halved, until step is 0
    */
    FlacSubframeParams binaryIdealLpcOrder(
        const std::vector<std::vector<float>>& coeffs,
        const std::vector<sample_t>& data,
        const FlacEncodeOptions& options);

}

#endif