#include <Convolve.h>
#include <assert.h>
#include <string.h>

#ifdef __SSE3__
#  include <xmmintrin.h>
#  include <pmmintrin.h>
#endif


#define TD_SIZE (block_size * 2)
#define FD_SIZE (block_size + 1)

ConvolverContext::ConvolverContext(size_t block_size)
  :block_size(block_size)
{
  float *td = fftwf_alloc_real(TD_SIZE);
  fftwf_complex *fd = fftwf_alloc_complex(FD_SIZE);
  
  td_to_fd = fftwf_plan_dft_r2c_1d(TD_SIZE, td, fd, FFTW_DESTROY_INPUT | FFTW_MEASURE);
  fd_to_td = fftwf_plan_dft_c2r_1d(TD_SIZE, fd, td, FFTW_DESTROY_INPUT | FFTW_MEASURE);
  
  fftwf_free(td);
  fftwf_free(fd);
}

ConvolverContext::~ConvolverContext()
{
  fftwf_destroy_plan(td_to_fd);
  fftwf_destroy_plan(fd_to_td);
}

Filter::Filter(ConvolverContext *context, size_t block_size, size_t filter_length, float *coefficients)
  :block_size(block_size)
{
  assert(context->block_size == block_size);

  float *td = fftwf_alloc_real(TD_SIZE);
  
  for (size_t offset = 0; offset < filter_length; offset += block_size)
  {
    size_t this_block_size = offset + block_size < filter_length ? block_size : filter_length - offset;
    
    // copy to the first half (or less) of td, and zero the second half
    memcpy(td, coefficients + offset, this_block_size * sizeof(*coefficients));
    memset(td + this_block_size, 0, (TD_SIZE - this_block_size) * sizeof(*coefficients));
    
    // fft into fd
    fftwf_complex *fd = fftwf_alloc_complex(FD_SIZE);
    fftwf_execute_dft_r2c(context->td_to_fd, td, fd);
    
    blocks.emplace_back(
        std::unique_ptr<fftwf_complex, void (*)(void*)>(fd, &fftwf_free)
    );
  }
  
  free(td);
}

#define SPECTRA_IDX(i) ((spectra_ofs + (i)) % num_blocks)
#define FILTER_IDX(i) ((filter_ofs + (i)) % (num_blocks + 1))

BlockConvolver::BlockConvolver(ConvolverContext *context, size_t block_size, size_t num_blocks)
  :context(context)
  ,block_size(block_size)
  ,num_blocks(num_blocks)
  ,filter_queue(num_blocks + 1, NULL)
  ,filter_ofs(0)
  ,spectra_ofs(0)
  ,current_td(fftwf_alloc_real(block_size * 2), &fftwf_free)
  ,multiply_out_a(fftwf_alloc_complex(block_size + 1), &fftwf_free)
  ,multiply_out_b(fftwf_alloc_complex(block_size + 1), &fftwf_free)
  ,out_td_a(fftwf_alloc_real(block_size * 2), &fftwf_free)
  ,out_td_b(fftwf_alloc_real(block_size * 2), &fftwf_free)
{
  assert(context->block_size == block_size);
  
  for (size_t i = 0; i < num_blocks; i++)
  {
    spectra_queue.emplace_back(
        std::unique_ptr<fftwf_complex, void (*)(void*)>(fftwf_alloc_complex(FD_SIZE), &fftwf_free)
    );
    
    memset(spectra_queue[i].get(), 0, FD_SIZE * sizeof(fftwf_complex));
  }
  
  memset(current_td.get(), 0, block_size * 2 * sizeof(float));
}

void BlockConvolver::crossfade_filter(Filter *filter)
{
  assert(filter->block_size == block_size);
  filter_queue[FILTER_IDX(0)] = filter;
}

void BlockConvolver::set_filter(Filter *filter)
{
  assert(filter->block_size == block_size);
  for (size_t i = 0; i < filter_queue.size(); i++)
    filter_queue[i] = filter;
}

void complex_mul_sum_cpp(fftwf_complex *out, fftwf_complex *a, fftwf_complex *b, size_t n)
{
  for (size_t i = 0; i < n; i++) {
    // implement out[i] += a[i] * b[i]:
    out[i][0] += a[i][0] * b[i][0] - a[i][1] * b[i][1];
    out[i][1] += a[i][0] * b[i][1] + a[i][1] * b[i][0];
  }
}

#ifdef __SSE3__
inline __m128 complex_mul_pair_ps(__m128 a, __m128 b)
{
    __m128 a_r = _mm_moveldup_ps(a);
    __m128 a_i = _mm_movehdup_ps(a);
    __m128 tmp1 = _mm_mul_ps(a_r, b);
    __m128 shufd = _mm_shuffle_ps(b, b, _MM_SHUFFLE(2,3,0,1));
    __m128 tmp2 = _mm_mul_ps(a_i, shufd);
    return _mm_addsub_ps(tmp1, tmp2);
}

void complex_mul_sum_sse(fftwf_complex *out, fftwf_complex *a, fftwf_complex *b, size_t n)
{
  assert(2 * sizeof(fftwf_complex) == 4 * sizeof(float));
  
  size_t offset;
  for (offset = 0; offset + 2 <= n; offset += 2) {
    __m128 av = _mm_load_ps((float *)(a + offset));
    __m128 bv = _mm_load_ps((float *)(b + offset));
    __m128 res = complex_mul_pair_ps(av, bv);
    
    __m128 out_old = _mm_load_ps((float *)(out + offset));
    __m128 out_new = _mm_add_ps(out_old, res);
    _mm_store_ps((float *)(out + offset), out_new);
  }
  
  complex_mul_sum_cpp(out + offset, a + offset, b + offset, n - offset);
}
#endif

void complex_mul_sum(fftwf_complex *out, fftwf_complex *a, fftwf_complex *b, size_t n)
{
#ifdef __SSE3__
  complex_mul_sum_sse(out, a, b, n);
#else
  complex_mul_sum_cpp(out, a, b, n);
#endif
}

void fade_output_norm(float *out, float *a, float *b, size_t n)
{
  for (size_t i = 0; i < n; i++)
  {
    float b_v = (float)i / (float)(n-1);
    float a_v = 1.0f - b_v;
    out[i] = (a_v * a[i] + b_v * b[i]) / ((float)n * 2.0f);
  }
}

void output_norm(float *out, float *a, size_t n)
{
  for (size_t i = 0; i < n; i++)
  {
    out[i] = a[i] / ((float)n * 2.0f);
  }
}

void BlockConvolver::filter_block(float *in, float *out)
{
  // copy in to second half of current_td
  memcpy(current_td.get() + block_size, in, block_size * sizeof(*current_td));
  // fft to spectra_queue[0]
  fftwf_execute_dft_r2c(context->td_to_fd, current_td.get(), spectra_queue[SPECTRA_IDX(0)].get());
  // copy in to first half of current_td
  memcpy(current_td.get(), in, block_size * sizeof(*current_td));
  
  // clear multiply_out_a
  memset(multiply_out_a.get(), 0, (block_size + 1) * sizeof(fftwf_complex));
  
  // multiply the spectra of all filters with the corresponding input spectra,
  // summing into multiply_out_a, for blocks where a crossfade is not needed.
  // set need_crossfade if any blocks that need crossfading are encountered.
  bool need_crossfade = false;
  for (size_t i = 0; i < num_blocks; i++)
  {
    Filter *old_filter = filter_queue[FILTER_IDX(i+1)];
    Filter *new_filter = filter_queue[FILTER_IDX(i  )];
    
    if (old_filter == new_filter) {
      complex_mul_sum(
          multiply_out_a.get(),
          new_filter->blocks[i].get(),
          spectra_queue[SPECTRA_IDX(i)].get(),
          block_size + 1);
    } else
      need_crossfade = true;
  }
  
  if (need_crossfade) {
    // copy multiply_out_a (containing all output blocks that don't need
    // crossfading) into multiply_out_b, then sum the old and new outputs into
    // multiply_out_a and multiply_out_b, respectively.
    memcpy(multiply_out_b.get(), multiply_out_a.get(), (block_size + 1) * sizeof(*multiply_out_a));
    for (size_t i = 0; i < num_blocks; i++)
    {
      Filter *old_filter = filter_queue[FILTER_IDX(i+1)];
      Filter *new_filter = filter_queue[FILTER_IDX(i  )];
      if (old_filter != new_filter) {
        complex_mul_sum(
            multiply_out_a.get(),
            old_filter->blocks[i].get(),
            spectra_queue[SPECTRA_IDX(i)].get(),
            block_size + 1);
        complex_mul_sum(
            multiply_out_b.get(),
            new_filter->blocks[i].get(),
            spectra_queue[SPECTRA_IDX(i)].get(),
            block_size + 1);
      }
    }

    // fft multiply_out_a and multiply_out_b into out_td_a and out_td_b, then
    // crossfade the second half of these to the output.
    fftwf_execute_dft_c2r(context->fd_to_td, multiply_out_a.get(), out_td_a.get());
    fftwf_execute_dft_c2r(context->fd_to_td, multiply_out_b.get(), out_td_b.get());
    fade_output_norm(out, out_td_a.get() + block_size, out_td_b.get() + block_size, block_size);
  } else {
    // all blocks were summed into multiply_out_a; fft this into out_td_a, and
    // copy the second half to the output.
    fftwf_execute_dft_c2r(context->fd_to_td, multiply_out_a.get(), out_td_a.get());
    output_norm(out, out_td_a.get() + block_size, block_size);
  }
  
  // older blocks are at higher indices, so move spectra_ofs and filter_ofs left, with wrap-around.
  spectra_ofs--;
  if (spectra_ofs < 0) spectra_ofs += num_blocks;
  filter_ofs--;
  if (filter_ofs < 0) filter_ofs += num_blocks + 1;
  
  // By default the next filter to use is the previous.
  filter_queue[FILTER_IDX(0)] = filter_queue[FILTER_IDX(1)];
}

#undef TD_SIZE
#undef FD_SIZE
#undef SPECTRA_IDX
#undef FILTER_IDX
