#define BOOST_TEST_MODULE test_foo
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <bbcat-base/misc.h>
#include "../src/BlockConvolver.h"

USE_BBC_AUDIOTOOLBOX

using free_type = void (*)(void*);

template <typename T>
std::unique_ptr<T, free_type> fftw_malloc_unique(size_t len) {
  std::unique_ptr<T, free_type> data((T*)fftwf_malloc(len * sizeof(T)), &fftwf_free);
  memset(data.get(), 0, len * sizeof(T));
  return data;
}

template <typename T>
std::unique_ptr<T, free_type> generate_random(size_t len, long int seedval=0)
{
  auto data = fftw_malloc_unique<T>(len);
  
  srand48(seedval);
  for (size_t i = 0; i < len; i++) {
    data.get()[i] = drand48();
  }
  
  return data;
}

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

void fade_up(float *out, float *in, size_t n) {
  float i_scale = 1.0 / (n-1);
  
  for (size_t i = 0; i < n; i++)
  {
    float a_v = (float)i * i_scale;
    out[i] = a_v * in[i];
  }
}

void fade_down(float *out, float *in, size_t n) {
  float i_scale = 1.0 / (n-1);
  
  for (size_t i = 0; i < n; i++)
  {
    float a_v = (float)i * i_scale;
    out[i] = (1 - a_v) * in[i];
  }
}

void mix_into(float *a, float *b, size_t n)
{
  for (size_t i = 0; i < n; i++)
    a[i] += b[i];
}

typedef std::pair<std::unique_ptr<float, free_type>, size_t> ir_spec;

struct ConvolutionTest {
  size_t block_size;
  size_t num_blocks;
  size_t len;
  int initial_block;
  
  std::vector<ir_spec> irs;
  std::vector<size_t> ir_for_block;
  
  std::unique_ptr<float, free_type> input;
  
  std::unique_ptr<float, free_type> test_output;
  std::unique_ptr<float, free_type> real_output;
  
  ConvolutionTest(size_t block_size, size_t num_blocks)
    :block_size(block_size)
    ,num_blocks(num_blocks)
    ,len(block_size * num_blocks)
    ,initial_block(-1)
    ,input(fftw_malloc_unique<float>(len))
    ,test_output(fftw_malloc_unique<float>(len))
    ,real_output(fftw_malloc_unique<float>(len))
  {
  }
  
  void run_test_convolve() {
    for (size_t i = 0; i < irs.size(); i++) {
      auto input_for_ir = fftw_malloc_unique<float>(len);
      
      for (size_t block = 0; block < num_blocks; block++) {
        bool this_block = i == ir_for_block[block];
        bool last_block = i == (block == 0 ? initial_block : ir_for_block[block - 1]);
        
        if (last_block && this_block)
          memcpy(input_for_ir.get() + block * block_size, input.get() + block * block_size, block_size * sizeof(float));
        else if (!last_block && this_block)
          fade_up(input_for_ir.get() + block * block_size, input.get() + block * block_size, block_size);
        else if (last_block && !this_block)
          fade_down(input_for_ir.get() + block * block_size, input.get() + block * block_size, block_size);
      }
      
      auto output_for_ir = fftw_malloc_unique<float>(len);
      convolve(len, output_for_ir.get(), input_for_ir.get(), irs[i].second, irs[i].first.get());
      mix_into(test_output.get(), output_for_ir.get(), len);
    }
  }
  
  void run_real_convolve(BlockConvolver::Context &ctx) {
    size_t max_num_blocks = 0;
    std::vector<BlockConvolver::Filter> filters;
    for (size_t i = 0; i < irs.size(); i++) {
      BlockConvolver::Filter filter(&ctx, block_size, irs[i].second, irs[i].first.get());
      if (filter.num_blocks() > max_num_blocks)
        max_num_blocks = filter.num_blocks();
      filters.push_back(std::move(filter));
    }
    
    BlockConvolver convolver(&ctx, block_size, max_num_blocks);
    
    if (initial_block != -1)
      convolver.set_filter(&filters[initial_block]);
    
    for (size_t block = 0; block < num_blocks; block++) {
      bool this_filter = ir_for_block[block];
      bool last_filter = block == 0 ? initial_block : ir_for_block[block - 1];
      
      if (last_filter != this_filter)
        convolver.crossfade_filter(&filters[this_filter]);
      
      convolver.filter_block(input.get() + block * block_size, real_output.get() + block * block_size);
    }
  }
  
  void run(BlockConvolver::Context &ctx) {
    run_test_convolve();
    run_real_convolve(ctx);
    
    for (size_t i = 0; i < len; i++)
      BOOST_CHECK_CLOSE(test_output.get()[i], real_output.get()[i], 0.15);
  }
};


BOOST_AUTO_TEST_CASE( filter_correct_num_blocks )
{
  BlockConvolver::Context ctx(512);
  auto coeff = generate_random<float>(2000);
  
  BOOST_CHECK_EQUAL(BlockConvolver::Filter(&ctx, 512, 1, coeff.get()).num_blocks(), 1);
  BOOST_CHECK_EQUAL(BlockConvolver::Filter(&ctx, 512, 511, coeff.get()).num_blocks(), 1);
  BOOST_CHECK_EQUAL(BlockConvolver::Filter(&ctx, 512, 512, coeff.get()).num_blocks(), 1);
  BOOST_CHECK_EQUAL(BlockConvolver::Filter(&ctx, 512, 513, coeff.get()).num_blocks(), 2);
}

BOOST_AUTO_TEST_CASE( single_block )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 1);
  t.initial_block = 0;
  t.irs.emplace_back(generate_random<float>(100, 0), 100);
  t.ir_for_block.emplace_back(0);
  t.input = generate_random<float>(512, 1);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( two_block )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 2);
  t.initial_block = 0;
  t.irs.emplace_back(generate_random<float>(512 * 3, 1), 512 * 3);
  t.ir_for_block = {0, 0};
  t.input = generate_random<float>(1024, 0);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( fade_once )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 3);
  t.initial_block = 0;
  t.irs.emplace_back(generate_random<float>(100, 1), 100);
  t.irs.emplace_back(generate_random<float>(512, 2), 512);
  t.ir_for_block = {0, 1, 1};
  t.input = generate_random<float>(512 * 3, 0);
  t.run(ctx);
}

BOOST_AUTO_TEST_CASE( fade_at_start_from_filter )
{
  BlockConvolver::Context ctx(512);
  ConvolutionTest t(512, 2);
  t.initial_block = 0;
  t.irs.emplace_back(generate_random<float>(512, 1), 512);
  t.irs.emplace_back(generate_random<float>(512, 2), 512);
  t.ir_for_block = {1, 1};
  t.input = generate_random<float>(512 * 2, 0);
  t.run(ctx);
}
