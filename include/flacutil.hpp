/*
flacutil.hpp
Yaakov Schectman, 2021
FLAC encoder
*/

#ifndef _FLACUTIL_HPP
#define _FLACUTIL_HPP

#include <cstdlib>
#include <cstdint>
#include <sstream>
#include <utility>
#include <tuple>

#define FLAC_MAX_LPC 32
#define FLAX_MAX_K 30

namespace Flac {
    
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
    
    using lpc_t = std::int16_t;
    using sample_t = std::int32_t;
    using zip_t = std::uint32_t;
    using residue_t = std::uint16_t;
    
    struct FlacEncodeOptions {
        public:
            size_t blockSize;
            size_t bitsPerSample;
            int minPred;
            int maxPred;
            int minPart;
            int maxPart;
            flag_t flags;
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
            std::vector<int> kParams;
            std::vector<lpc_t> coefficients;
    };
    
    class FlacSubframe {
        private:
            FlacSubframeParams params;
            std::vector<sample_t> rawData; // TODO should be reference?
            std::vector<zip_t> zippedResidue;
        public:
            /*
            Construct from unencoded data with provided encoding options
            
            First, check if it should be encoded as a constant subblock
            If not, and if using dynamic LPC, calculate LPC coefficients
            Calculate and store the verbatim bit length
            If using fixed LPC, calculate and zip the residue for each order,
            and find and store each order's bit length and parameters
            If using dynamic LPC and estimating LPC order, use the highest order LPC that uses
            all of its coefficients
            Otherwise, using the lpcMethod, try each candidate LPC order storing their bit lengths
            and parameters
            Store the parameters for whichever attempted method yields the least bit length
            */
            FlacSubframe(std::vector<sample_t>& data, FlacEncodeOptions& options);
    };
    
    /*
    For a given prediction order pred and partition order part, return a nested vector X such that
    X[k][i] is the number of bits it would take to encode the i-th partition of residue with rice
    paramter k,
    unless maxK == 0, in which case
    X[0][i] is the sum of all residues in the i-th parameter
    */
    std::vector<std::vector<zip_t>> sumsForEachPartition(
        std::vector<zip_t>& residue,
        int maxK, int pred, int part);
    
    /*
    Change sums to be for one lower partition order
    */
    void decreasePartitionOfSums(std::vector<std::vector<zip_t>>& sums, int maxK);
    
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
    inline size_t sizeBitsOfSum(zip_t sum, size_t n, int k)
    {
        return (k + 1) * n + (sum >> k);
    }
    
    /*
    Find which k minimizes count(sums[k][partition])
    With useEstimation, count(n) = sizeBitsOfSum
    Otherwise, count(n) = n
    
    Return <k, bit length>
    */
    std::pair<int, size_t> idealK(
        std::vector<std::vector<zip_t>>& sums,
        int partition, bool useEstimation);
    
    /*
    Returns <partition order, k paramter for each partition, length in bits>
    
    First, calculate the sums for each partition using the max order.
    Then, for each candidate partition order, in descending order, iterate over those sums to find
    the ideal k parameter for that partition, and estimate the number of bits needed to encode it.
    Return whichever partition order minimizes the resultant bits, as well as its k-vector, as
    <partition order, list of k params, bit length>
    */
    std::tuple<int, std::vector<int>, size_t> idealRiceEncoding(
        std::vector<zip_t>& data,
        int pred, int minPart, int maxPart,
        bool useEstimation);

    /*
    Calculate and return the autocorrelations of lags [0, maxLag]
    */
    std::vector<float> autocorrelation(std::vector<sample_t>& data, size_t maxLag);

    /*
    Calculate the LPC coefficients of n samples, stored in data, using Levinson-Durbin recursion
    */
    std::vector<std::vector<float>> calcLPC(std::vector<sample_t>& data);

}

#endif