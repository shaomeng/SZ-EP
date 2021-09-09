#include "EP.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numeric>


auto EP::read_n_bytes(const char* filename, size_t n_bytes, void* buffer) -> bool
{
  std::unique_ptr<std::FILE, decltype(&std::fclose)> fp(std::fopen(filename, "rb"), &std::fclose);

  if (!fp)
    return false;

  std::fseek(fp.get(), 0, SEEK_END);
  if (std::ftell(fp.get()) < n_bytes)
    return false;

  std::fseek(fp.get(), 0, SEEK_SET);
  if (std::fread(buffer, 1, n_bytes, fp.get()) != n_bytes)
    return false;

  return true;
}

template <typename T>
auto EP::read_whole_file(const char* filename) -> std::vector<T>
{
  std::vector<T> buf;

  std::unique_ptr<std::FILE, decltype(&std::fclose)> fp(std::fopen(filename, "rb"), &std::fclose);
  if (!fp)
    return buf;

  std::fseek(fp.get(), 0, SEEK_END);
  const size_t file_size = std::ftell(fp.get());
  const size_t num_vals = file_size / sizeof(T);
  std::fseek(fp.get(), 0, SEEK_SET);

  buf.resize(num_vals);
  size_t nread = std::fread(buf.data(), sizeof(T), num_vals, fp.get());
  if (nread != num_vals)
    buf.clear();

  return buf;
}
template auto EP::read_whole_file(const char*) -> std::vector<float>;
template auto EP::read_whole_file(const char*) -> std::vector<double>;
template auto EP::read_whole_file(const char*) -> std::vector<uint8_t>;

auto EP::write_n_bytes(const char* filename, size_t n_bytes, const void* buffer) -> bool
{
  std::unique_ptr<std::FILE, decltype(&std::fclose)> fp(std::fopen(filename, "wb"), &std::fclose);
  if (!fp)
    return false;

  if (std::fwrite(buffer, 1, n_bytes, fp.get()) != n_bytes)
    return false;
  else
    return true;
}

template <typename T>
void EP::calc_stats(const T* arr1,
                       const T* arr2,
                       size_t len,
                       T& rmse,
                       T& linfty,
                       T& psnr,
                       T& arr1min,
                       T& arr1max)
{
  const size_t stride_size = 4096;
  const size_t num_of_strides = len / stride_size;
  const size_t remainder_size = len - stride_size * num_of_strides;

  //
  // Calculate min and max of arr1
  //
  const auto minmax = std::minmax_element(arr1, arr1 + len);
  arr1min = *minmax.first;
  arr1max = *minmax.second;

  //
  // In rare cases, the two input arrays are identical.
  //
  auto mism = std::mismatch(arr1, arr1 + len, arr2, arr2 + len);
  if (mism.first == arr1 + len && mism.second == arr2 + len) {
    rmse = 0.0;
    linfty = 0.0;
    psnr = std::numeric_limits<T>::infinity();
    return;
  }

  auto sum_vec = std::vector<T>(num_of_strides + 1);
  auto linfty_vec = std::vector<T>(num_of_strides + 1);

  //
  // Calculate summation and l-infty of each stride
  //
  // (Uncomment the following line to enable OpenMP)
  // #pragma omp parallel for
  for (size_t stride_i = 0; stride_i < num_of_strides; stride_i++) {
    T linfty = 0.0;
    auto buf = std::array<T, stride_size>();
    for (size_t i = 0; i < stride_size; i++) {
      const size_t idx = stride_i * stride_size + i;
      auto diff = std::abs(arr1[idx] - arr2[idx]);
      linfty = std::max(linfty, diff);
      buf[i] = diff * diff;
    }
    sum_vec[stride_i] = std::accumulate(buf.begin(), buf.end(), T{0.0});
    linfty_vec[stride_i] = linfty;
  }

  //
  // Calculate summation and l-infty of the remaining elements
  //
  T last_linfty = 0.0;
  auto last_buf = std::array<T, stride_size>{};  // must be enough for
                                                 // `remainder_size` elements.
  for (size_t i = 0; i < remainder_size; i++) {
    const size_t idx = stride_size * num_of_strides + i;
    auto diff = std::abs(arr1[idx] - arr2[idx]);
    last_linfty = std::max(last_linfty, diff);
    last_buf[i] = diff * diff;
  }
  sum_vec[num_of_strides] =
      std::accumulate(last_buf.begin(), last_buf.begin() + remainder_size, T{0.0});
  linfty_vec[num_of_strides] = last_linfty;

  //
  // Now calculate linfty
  //
  linfty = *(std::max_element(linfty_vec.begin(), linfty_vec.end()));

  //
  // Now calculate rmse and psnr
  // Note: psnr is calculated in dB, and follows the equation described in:
  // http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/VELDHUIZEN/node18.html
  // Also refer to https://www.mathworks.com/help/vision/ref/psnr.html
  //
  const auto msr = std::accumulate(sum_vec.begin(), sum_vec.end(), T{0.0}) / T(len);
  rmse = std::sqrt(msr);
  auto range_sq = *minmax.first - *minmax.second;
  range_sq *= range_sq;
  psnr = std::log10(range_sq / msr) * T{10.0};
}
template void
EP::calc_stats(const float*, const float*, size_t, float&, float&, float&, float&, float&);
template void EP::calc_stats(const double*,
                                const double*,
                                size_t,
                                double&,
                                double&,
                                double&,
                                double&,
                                double&);

auto EP::chunk_volume(const std::array<size_t, 3>& vol_dim,
                         const std::array<size_t, 3>& chunk_dim)
    -> std::vector<std::array<size_t, 6>>
{
  // Step 1: figure out how many segments are there along each axis.
  auto n_segs = std::array<size_t, 3>();
  for (size_t i = 0; i < 3; i++) {
    n_segs[i] = vol_dim[i] / chunk_dim[i];
    if ((vol_dim[i] % chunk_dim[i]) > (chunk_dim[i] / 2))
      n_segs[i]++;
    // In case the volume has one dimension that's too small, let's make it
    // at least one segment anyway.
    if (n_segs[i] == 0)
      n_segs[i] = 1;
  }

  // Step 2: calculate the starting indices of each segment along each axis.
  auto x_tics = std::vector<size_t>(n_segs[0] + 1);
  for (size_t i = 0; i < n_segs[0]; i++)
    x_tics[i] = i * chunk_dim[0];
  x_tics[n_segs[0]] = vol_dim[0];  // the last tic is the length in X

  auto y_tics = std::vector<size_t>(n_segs[1] + 1);
  for (size_t i = 0; i < n_segs[1]; i++)
    y_tics[i] = i * chunk_dim[1];
  y_tics[n_segs[1]] = vol_dim[1];  // the last tic is the length in Y

  auto z_tics = std::vector<size_t>(n_segs[2] + 1);
  for (size_t i = 0; i < n_segs[2]; i++)
    z_tics[i] = i * chunk_dim[2];
  z_tics[n_segs[2]] = vol_dim[2];  // the last tic is the length in Z

  // Step 3: fill in details of each chunk
  auto n_chunks = n_segs[0] * n_segs[1] * n_segs[2];
  auto chunks = std::vector<std::array<size_t, 6>>(n_chunks);
  size_t chunk_idx = 0;
  for (size_t z = 0; z < n_segs[2]; z++)
    for (size_t y = 0; y < n_segs[1]; y++)
      for (size_t x = 0; x < n_segs[0]; x++) {
        chunks[chunk_idx][0] = x_tics[x];                  // X start
        chunks[chunk_idx][1] = x_tics[x + 1] - x_tics[x];  // X length
        chunks[chunk_idx][2] = y_tics[y];                  // Y start
        chunks[chunk_idx][3] = y_tics[y + 1] - y_tics[y];  // Y length
        chunks[chunk_idx][4] = z_tics[z];                  // Z start
        chunks[chunk_idx][5] = z_tics[z + 1] - z_tics[z];  // Z length
        chunk_idx++;
      }

  return chunks;
}

template <typename T>
auto EP::gather_chunk(const T* vol, dims_type vol_dim, const std::array<size_t, 6>& chunk)
    -> std::vector<T>
{
  auto buf = std::vector<T>();
  if (chunk[0] + chunk[1] > vol_dim[0] || chunk[2] + chunk[3] > vol_dim[1] ||
      chunk[4] + chunk[5] > vol_dim[2])
    return buf;

  auto len = chunk[1] * chunk[3] * chunk[5];
  buf.resize(len);

  size_t idx = 0;
  for (size_t z = chunk[4]; z < chunk[4] + chunk[5]; z++) {
    const size_t plane_offset = z * vol_dim[0] * vol_dim[1];
    for (size_t y = chunk[2]; y < chunk[2] + chunk[3]; y++) {
      const size_t col_offset = plane_offset + y * vol_dim[0];
      for (size_t x = chunk[0]; x < chunk[0] + chunk[1]; x++)
        buf[idx++] = vol[col_offset + x];
    }
  }

  // Will be subject to Named Return Value Optimization.
  return buf;
}
template auto EP::gather_chunk(const float*, dims_type, const std::array<size_t, 6>&)
    -> std::vector<float>;
template auto EP::gather_chunk(const double*, dims_type, const std::array<size_t, 6>&)
    -> std::vector<double>;

template <typename T>
void EP::scatter_chunk(std::vector<T>& big_vol,
                          dims_type vol_dim,
                          const std::vector<T>& small_vol,
                          const std::array<size_t, 6>& chunk)
{
  size_t idx = 0;
  for (size_t z = chunk[4]; z < chunk[4] + chunk[5]; z++) {
    const size_t plane_offset = z * vol_dim[0] * vol_dim[1];
    for (size_t y = chunk[2]; y < chunk[2] + chunk[3]; y++) {
      const size_t col_offset = plane_offset + y * vol_dim[0];
      for (size_t x = chunk[0]; x < chunk[0] + chunk[1]; x++)
        big_vol[col_offset + x] = small_vol[idx++];
    }
  }
}
template void EP::scatter_chunk(vecf_type&, dims_type, 
                                const vecf_type&, const std::array<size_t, 6>&);
template void EP::scatter_chunk(vecd_type&, dims_type, 
                                const vecd_type&, const std::array<size_t, 6>&);
