#include "EP.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <string> // std::stol()

int main(int argc, char* argv[])
{
  if( argc != 5 ) {
    std::cout << "Usage: " << argv[0] << " input_name dim_x dim_y dim_z" << std::endl;
    return 0;
  }

  // Hard code the desired chunk size, and num of threads
  const auto chunk_dims = EP::dims_type{ 128, 128, 128 };

  const char* in_name = argv[1];
  const auto in_dims = EP::dims_type{ size_t(std::stol(argv[2])), 
                                      size_t(std::stol(argv[3])), 
                                      size_t(std::stol(argv[4])) };
  const auto total_len = in_dims[0] * in_dims[1] * in_dims[2];

  // Read in the whole volume in float
  auto in_buf = EP::read_whole_file<float>( in_name );
  assert( in_buf.size() == total_len );
  auto out_buf = std::vector<float>( total_len );

  // Calculate the offset and size of each chunk;
  const auto chunk_def = EP::chunk_volume(in_dims, chunk_dims);
  const auto num_chunks = chunk_def.size();
  auto chunks = std::vector< std::vector<float> >( num_chunks );
  

  #pragma omp parallel for num_threads(4)
  for( size_t i = 0; i < num_chunks; i++ ) {

    // Gather a chunk from the big volume
    chunks[i] = EP::gather_chunk( in_buf.data(), in_dims, chunk_def[i] );

    // Put this chunk to the big volume
    EP::scatter_chunk( out_buf, in_dims, chunks[i], chunk_def[i] );
  }

  auto is_equal = std::equal( in_buf.begin(), in_buf.end(), out_buf.begin() );
  std::cout << std::boolalpha << is_equal << std::endl;
  
}
