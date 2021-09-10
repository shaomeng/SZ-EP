#include "EP.h"

#include <algorithm>
#include <numeric> // std::accumulate()
#include <cassert>
#include <iostream>
#include <string> // std::stol()

using EP::vecd_type;
using EP::vecf_type;
using EP::vec8_type;
using EP::dims_type;


// Declare SZ functions that will be used
extern "C" {
  int SZ_Init( const char* config_file );
  unsigned char* SZ_compress_args(int dataType, void *data, size_t *outSize, int errBoundMode, 
                                  double absErrBound, double relBoundRatio, double pwrBoundRatio, 
                                  size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
  void *SZ_decompress(int dataType, unsigned char *bytes, size_t byteLength, 
                      size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
  void SZ_Finalize();
};


auto test_sz(  vecf_type&        chunk,
               double            tol,
               dims_type         dims ) -> std::pair<size_t, vecf_type>
{
  // Note that the following definitions are from include/sz/defines.h
  const int SZ_FLOAT = 0; // #define SZ_FLOAT 0
  const int ABS = 0;      // #define ABS 0

  size_t compressed_buf_len = 0;

  auto* compressed_buf = SZ_compress_args( SZ_FLOAT, chunk.data(), &compressed_buf_len,
                                           ABS, tol, tol, tol, 
                                           0, 0, dims[2], dims[1], dims[0] );
  const auto* reconstructed_buf = SZ_decompress( SZ_FLOAT, compressed_buf, compressed_buf_len,
                                                 0, 0, dims[2], dims[1], dims[0] );

  const float* p = reinterpret_cast<const float*>(reconstructed_buf);
  const auto total_len = dims[0] * dims[1] * dims[2];
  auto buf = vecf_type( total_len );
  std::copy( p, p + total_len, buf.begin() );
  
  return {compressed_buf_len, buf};
}


int main(int argc, char* argv[])
{
  if( argc != 6 ) {
    std::cout << "Usage: " << argv[0] << " input_name dim_x dim_y dim_z tolerance" << std::endl;
    return 0;
  }

  // Hard code the desired chunk size, and num of threads
  const auto chunk_dims = dims_type{ 128, 128, 128 };

  const char* in_name = argv[1];
  const auto in_dims = dims_type{ std::stoul(argv[2]), std::stoul(argv[3]), std::stoul(argv[4]) };
  const double tol = std::stod( argv[5] );
  const auto total_len = in_dims[0] * in_dims[1] * in_dims[2];

  // Read in the whole volume in float
  auto in_buf = EP::read_whole_file<float>( in_name );
  assert( in_buf.size() == total_len );
  auto out_buf = std::vector<float>( total_len );

  // Calculate the offset and size of each chunk;
  const auto chunk_def = EP::chunk_volume(in_dims, chunk_dims);
  const auto num_chunks = chunk_def.size();
  auto chunks = std::vector< std::vector<float> >( num_chunks );
  auto comp_len = std::vector<size_t>( num_chunks );
  
  SZ_Init( "./sz.config" ); // Content in this file will be overwritten.

  #pragma omp parallel for num_threads(2)
  for( size_t i = 0; i < num_chunks; i++ ) {

    // Gather a chunk from the big volume
    chunks[i] = EP::gather_chunk( in_buf.data(), in_dims, chunk_def[i] );

    // Test it using SZ
    auto result = test_sz( chunks[i], tol, 
                           {chunk_def[i][1], chunk_def[i][3], chunk_def[i][5]} );

    // Put this chunk to the output big volume
    EP::scatter_chunk( out_buf, in_dims, result.second, chunk_def[i] );

    // Also collect the compressed size
    comp_len[i] = result.first;
  }

  SZ_Finalize();

  // Gather statistics
  float rmse, linfty, psnr, arr1min, arr1max;
  EP::calc_stats(in_buf.data(), out_buf.data(), total_len,
                 rmse, linfty, psnr, arr1min, arr1max);
  auto total_bytes = std::accumulate( comp_len.begin(), comp_len.end(), size_t{0} );
  printf("-> SZ achieved bpp = %.2f\n", double(total_bytes * 8) / double(total_len) );
  printf("|- Original data range = (%.2e, %.2e)\n", arr1min, arr1max);
  printf("|- Reconstructed data RMSE = %.2e, L-Infty = %.2e, PSNR = %.2fdB\n", rmse, linfty, psnr);

}
