#ifndef EP_H
#define SEP_H

//
// We put common functions that are used across the speck encoder here.
//

#include <array>
#include <cstddef>  // std::size_t
#include <cstdint>  // fixed width integers
#include <cstdlib>
#include <iterator>
#include <memory>
#include <utility>  // std::pair
#include <vector>

namespace EP {

using std::size_t;  // Seems most appropriate

//
// A few shortcuts
//
using vecd_type = std::vector<double>;
using vecf_type = std::vector<float>;
using vec8_type = std::vector<uint8_t>;
using dims_type = std::array<size_t, 3>;

//
// Helper functions
//
// Read from and write to a file
auto write_n_bytes(const char* filename, size_t n_bytes, const void* buffer) -> bool;
auto read_n_bytes(const char* filename, size_t n_bytes, void* buffer) -> bool;
template <typename T>
auto read_whole_file(const char* filename) -> std::vector<T>;

// Calculate a suite of statistics
// Note that arr1 is considered as the ground truth array, so it's the range of
// arr1 that is used internally for psnr calculations.
template <typename T>
void calc_stats(const T* arr1,
                const T* arr2,
                size_t len,
                T& rmse,
                T& linfty,
                T& psnr,
                T& arr1min,
                T& arr1max);

// Given a whole volume size and a desired chunk size, this helper function
// returns a list of chunks specified by 6 integers: 
// chunk[0], [2], [4]: starting index of this chunk;
// chunk[1], [3], [5]: length of this chunk
auto chunk_volume(const dims_type& vol_dim, const dims_type& chunk_dim)
    -> std::vector<std::array<size_t, 6>>;

// Gather a chunk from a bigger volume
template <typename T>
auto gather_chunk(const T* vol, dims_type vol_dim, const std::array<size_t, 6>& chunk) 
    -> std::vector<T>;

// Put this chunk to a bigger volume
template <typename T>
void scatter_chunk(std::vector<T>& big_vol,
                   dims_type vol_dim,
                   const std::vector<T>& small_vol,
                   const std::array<size_t, 6>& chunk);

};  // namespace EP

#endif
