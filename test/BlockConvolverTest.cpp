#define BOOST_TEST_MODULE BlockConvolverTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <map>
#include <fstream>

#include <bbcat-base/misc.h>
#include "../src/BlockConvolver.h"

USE_BBC_AUDIOTOOLBOX

using free_type = void (*)(void*);

template <typename T>
std::unique_ptr<T[], free_type> fftw_malloc_unique(size_t len) {
  std::unique_ptr<T[], free_type> data((T*)fftwf_malloc(len * sizeof(T)), &fftwf_free);
  memset(data.get(), 0, len * sizeof(T));
  return data;
}

template <typename T>
std::unique_ptr<T[], free_type> generate_random(size_t len, size_t num_nonzero, long int seedval=0)
{
  auto data = fftw_malloc_unique<T>(len);
  
  srand48(seedval);
  for (size_t i = 0; i < num_nonzero; i++) {
    data[lrand48() % len] = drand48();
  }
  
  return data;
}

/// Simple convolution implementation.
/// @param len output and input lengths
/// @param out[out] output data
/// @param out[in] input data
/// @param ir_len length of impulse response
/// @param ir impulse response data
void convolve(size_t len, float *out, float *in, size_t ir_len, float *ir)
{
  for (size_t i = 0; i < len; i++) {
    float total = 0.0;
    for (size_t j = 0; j < ir_len; j++) {
      if ((int)i - (int)j < 0) break;
      total += in[i - j] * ir[j];
    }
    out[i] = total;
  }
}

/// Fade up linearly over a block.
void fade_up(float *out, float *in, size_t n) {
  float i_scale = 1.0 / (n-1);
  
  for (size_t i = 0; i < n; i++)
  {
    float a_v = (float)i * i_scale;
    out[i] = a_v * in[i];
  }
}

/// Fade down linearly over a block.
void fade_down(float *out, float *in, size_t n) {
  float i_scale = 1.0 / (n-1);
  
  for (size_t i = 0; i < n; i++)
  {
    float a_v = (float)i * i_scale;
    out[i] = (1 - a_v) * in[i];
  }
}

/// Mix b into a.
void mix_into(float *a, float *b, size_t n)
{
  for (size_t i = 0; i < n; i++)
    a[i] += b[i];
}

/// Does a buffer contain only zeros?
bool all_zeros(float *in, size_t len) {
  for (size_t i = 0; i < len; i++)
    if (in[i] != 0.0f)
      return false;
  
  return true;
}

/// Write a float array to a file, readable by numpy.loadtext
void write_array(const std::string &base, const std::string &name, size_t n, float *x) {
  if (base != "") {
    std::ofstream out(base + name);
    assert(out.is_open());
    for (size_t i = 0; i < n; i++)
      out << x[i] << std::endl;
  }
}

struct ConvolutionTest {
  /// block size for BlockConvolver
  size_t block_size;
  /// number of blocks to process
  size_t num_blocks;
  /// number of samples to process
  size_t len;
  /// impulse response to set before starting the test; don't set (creating a fade up) if negative)
  int initial_ir;
  /// the number of blocks to create the BlockConvolver with; will be set to
  /// the maximum filter length used if it's too low.
  size_t max_num_blocks;
  /// pass NULL for blocks containing all zero samples?
  bool null_for_zeros;
  /// Which constructor to use for BlockConvolver
  enum {
    NO_FILTER,                 ///< standard, without a filter
    WITH_FILTER_NUM_BLOCKS,    ///< specify a filter and a number of blocks
    WITH_FILTER_NO_NUM_BLOCKS  ///< specify just a filter
  } constructor;
  
  typedef std::pair<std::unique_ptr<float[], free_type>, size_t> ir_spec;
  /// available inpulse responses; pairs of a float pointer and a length
  std::vector<ir_spec> irs;
  /// Which ir index to use for each input block; should be num_blocks long.
  /// Negative numbers indicate no filter.
  std::vector<int> ir_for_block;
  
  /// Input data; len long.
  std::unique_ptr<float[], free_type> input;
  
  /// Output data from the test implementation; len long.
  std::unique_ptr<float[], free_type> test_output;
  /// Output data from the real implementation; len long.
  std::unique_ptr<float[], free_type> real_output;
  
  ConvolutionTest(size_t block_size, size_t num_blocks)
    :block_size(block_size)
    ,num_blocks(num_blocks)
    ,len(block_size * num_blocks)
    ,initial_ir(-1)
    ,max_num_blocks(0)
    ,null_for_zeros(false)
    ,constructor(NO_FILTER)
    ,input(fftw_malloc_unique<float>(len))
    ,test_output(fftw_malloc_unique<float>(len))
    ,real_output(fftw_malloc_unique<float>(len))
  {
  }
  
  void run_test_convolve(const std::string &out_dir="") {
    // For each filter, create modified input data with fades applied, then convolve and add to the output.
    for (size_t i = 0; i < irs.size(); i++) {
      auto input_for_ir = fftw_malloc_unique<float>(len);
      
      // Apply fades for each block
      for (size_t block = 0; block < num_blocks; block++) {
        // Was this filter active in the last block , or this block?
        bool last_block = (int)i == (block == 0 ? initial_ir : ir_for_block[block - 1]);
        bool this_block = (int)i == ir_for_block[block];
        
        // appropriate fades for this block
        if (last_block && this_block)
          memcpy(input_for_ir.get() + block * block_size, input.get() + block * block_size, block_size * sizeof(float));
        else if (!last_block && this_block)
          fade_up(input_for_ir.get() + block * block_size, input.get() + block * block_size, block_size);
        else if (last_block && !this_block)
          fade_down(input_for_ir.get() + block * block_size, input.get() + block * block_size, block_size);
      }
      
      write_array(out_dir, "input_" + std::to_string(i) + ".dat", len, input_for_ir.get());
      
      // convolve and mix
      auto output_for_ir = fftw_malloc_unique<float>(len);
      convolve(len, output_for_ir.get(), input_for_ir.get(), irs[i].second, irs[i].first.get());
      mix_into(test_output.get(), output_for_ir.get(), len);
    }
  }
  
  void run_real_convolve(BlockConvolver::Context &ctx) {
    // Create filters from irs, while finding the required number of blocks.
    std::vector<BlockConvolver::Filter> filters;
    for (size_t i = 0; i < irs.size(); i++) {
      BlockConvolver::Filter filter(&ctx, irs[i].second, irs[i].first.get());
      if (filter.num_blocks() > max_num_blocks) max_num_blocks = filter.num_blocks();
      filters.push_back(std::move(filter));
    }
    
    std::unique_ptr<BlockConvolver> convolver;
    
    // Create the filter
    switch (constructor) {
      case NO_FILTER:
        convolver.reset(new BlockConvolver(&ctx, max_num_blocks));
        break;
      case WITH_FILTER_NUM_BLOCKS:
        convolver.reset(new BlockConvolver(&ctx, &filters[initial_ir], max_num_blocks));
        break;
      case WITH_FILTER_NO_NUM_BLOCKS:
        convolver.reset(new BlockConvolver(&ctx, &filters[initial_ir]));
        break;
    }
    
    // don't set the starting filter if initial_ir is invalid
    if (initial_ir >= 0)
      convolver->set_filter(&filters[initial_ir]);
    
    // filter each block
    for (size_t block = 0; block < num_blocks; block++) {
      // crossfade if the filter has changed
      int this_filter = ir_for_block[block];
      int last_filter = block == 0 ? initial_ir : ir_for_block[block - 1];
      if (last_filter != this_filter) {
        if (this_filter >= 0)
          convolver->crossfade_filter(&filters[this_filter]);
        else
          convolver->crossfade_filter(NULL);
      }
      
      // filter the block to the output; possibly passing null if zeros are detected
      if (null_for_zeros && all_zeros(input.get() + block * block_size, block_size)) {
        convolver->filter_block(NULL, real_output.get() + block * block_size);
      } else {
        convolver->filter_block(input.get() + block * block_size, real_output.get() + block * block_size);
      }
    }
  }
  
  /** Run the test, possibly writing the input and output datas to files.
   *  @param out_dir Output dir for data files; should end with a trailing slash.
   */
  void run(BlockConvolver::Context &ctx, const std::string &out_dir="") {
    run_test_convolve(out_dir);
    run_real_convolve(ctx);
    
    write_array(out_dir, "input.dat", len, input.get());
    write_array(out_dir, "test_output.dat", len, test_output.get());
    write_array(out_dir, "real_output.dat", len, real_output.get());
    
    for (size_t i = 0; i < len; i++)
      BOOST_CHECK_SMALL(test_output.get()[i] - real_output.get()[i], 1e-6f);
  }
};


BOOST_AUTO_TEST_CASE( filter_correct_num_blocks )
{
  BlockConvolver::Context ctx(512);
  auto coeff = generate_random<float>(2000, 100);
  
  BOOST_CHECK_EQUAL(BlockConvolver::Filter(&ctx, 1, coeff.get()).num_blocks(), 1);
  BOOST_CHECK_EQUAL(BlockConvolver::Filter(&ctx, 511, coeff.get()).num_blocks(), 1);
  BOOST_CHECK_EQUAL(BlockConvolver::Filter(&ctx, 512, coeff.get()).num_blocks(), 1);
  BOOST_CHECK_EQUAL(BlockConvolver::Filter(&ctx, 513, coeff.get()).num_blocks(), 2);
}

BOOST_AUTO_TEST_CASE( single_block )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 1);
  t.irs.emplace_back(generate_random<float>(100, 10, 1), 100);
  t.initial_ir = 0;
  t.ir_for_block = {0};
  t.input = generate_random<float>(512, 200, 0);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( two_blocks )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 2);
  t.irs.emplace_back(generate_random<float>(512 * 3, 20, 1), 512 * 3);
  t.initial_ir = 0;
  t.ir_for_block = {0, 0};
  t.input = generate_random<float>(512 * 2, 300, 0);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( fade_once )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 3);
  t.irs.emplace_back(generate_random<float>(100, 10, 1), 100);
  t.irs.emplace_back(generate_random<float>(512, 10, 2), 512);
  t.initial_ir = 0;
  t.ir_for_block = {0, 1, 1};
  t.input = generate_random<float>(512 * 3, 300, 0);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( fade_at_start_from_silence )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 2);
  t.irs.emplace_back(generate_random<float>(512, 10, 1), 512);
  t.ir_for_block = {0, 0};
  t.input = generate_random<float>(512 * 2, 300, 0);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( fade_to_silence )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 3);
  t.irs.emplace_back(generate_random<float>(512, 10, 1), 512);
  t.initial_ir = 0;
  t.ir_for_block = {0, -1, -1};
  t.input = generate_random<float>(512 * 3, 300, 0);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( fade_from_silence )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 3);
  t.irs.emplace_back(generate_random<float>(512, 10, 1), 512);
  t.ir_for_block = {-1, 0, 0};
  t.input = generate_random<float>(512 * 3, 300, 0);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( fade_at_start_from_filter )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 2);
  t.irs.emplace_back(generate_random<float>(512, 10, 1), 512);
  t.irs.emplace_back(generate_random<float>(512, 10, 2), 512);
  t.initial_ir = 0;
  t.ir_for_block = {1, 1};
  t.input = generate_random<float>(512 * 2, 300, 0);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( smaller_filter_than_convolver )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 4);
  t.irs.emplace_back(generate_random<float>(512 * 2, 20, 1), 512 * 2);
  t.initial_ir = 0;
  t.ir_for_block = {0, 0, 0, 0};
  t.input = generate_random<float>(512 * 4, 300, 0);
  t.max_num_blocks = 3;
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( different_num_blocks )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 4);
  t.irs.emplace_back(generate_random<float>(512 * 2, 20, 1), 512 * 2);
  t.irs.emplace_back(generate_random<float>(512 * 3, 20, 2), 512 * 3);
  t.initial_ir = 0;
  t.ir_for_block = {0, 1, 1, 1};
  t.input = generate_random<float>(512 * 4, 300, 0);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( zero_input_blocks )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 5);
  t.irs.emplace_back(generate_random<float>(512 * 2, 20, 1), 512 * 2);
  t.initial_ir = 0;
  t.ir_for_block = {0, 0, 0, 0, 0};
  t.input = generate_random<float>(512 * 5, 300, 0);
  memset(t.input.get() + 512, 0, 512 * 3 * sizeof(float));
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( null_input_blocks )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 5);
  t.irs.emplace_back(generate_random<float>(512 * 2, 20, 1), 512 * 2);
  t.initial_ir = 0;
  t.ir_for_block = {0, 0, 0, 0, 0};
  t.input = generate_random<float>(512 * 5, 300, 0);
  memset(t.input.get() + 512, 0, 512 * 3 * sizeof(float));
  t.null_for_zeros = true;
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( lots_of_filters )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 9);
  t.irs.emplace_back(generate_random<float>(512 * 2, 20, 1), 512 * 2);
  t.irs.emplace_back(generate_random<float>(512 * 3, 20, 2), 512 * 3);
  t.irs.emplace_back(generate_random<float>(512 * 1, 20, 3), 512 * 1);
  t.irs.emplace_back(generate_random<float>(512 * 4, 20, 4), 512 * 4);
  t.initial_ir = 0;
  t.ir_for_block = {0, 1, 2, 3, 3, 2, 2, 1, 0};
  t.input = generate_random<float>(512 * 9, 500, 0);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( construct_with_filter )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 3);
  t.irs.emplace_back(generate_random<float>(512 * 2, 20, 1), 512 * 2);
  t.initial_ir = 0;
  t.ir_for_block = {0, 0, 0};
  t.constructor = ConvolutionTest::WITH_FILTER_NUM_BLOCKS;
  t.input = generate_random<float>(512 * 3, 300, 0);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( construct_with_filter_no_blocks )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 3);
  t.irs.emplace_back(generate_random<float>(512 * 2, 20, 1), 512 * 2);
  t.initial_ir = 0;
  t.ir_for_block = {0, 0, 0};
  t.constructor = ConvolutionTest::WITH_FILTER_NO_NUM_BLOCKS;
  t.input = generate_random<float>(512 * 3, 300, 0);
  t.run(ctx);
}
