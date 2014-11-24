#include "BlockConvolver.h"
#include <assert.h>
#include <string.h>

#ifdef __SSE3__
#  include <xmmintrin.h>
#  include <pmmintrin.h>
#endif

BBC_AUDIOTOOLBOX_START

#define TD_SIZE (block_size * 2)
#define FD_SIZE (block_size + 1)

BlockConvolver::Context::Context(size_t block_size)
  :block_size(block_size)
{
  float *td = fftwf_alloc_real(TD_SIZE);
  fftwf_complex *fd = fftwf_alloc_complex(FD_SIZE);
  
  td_to_fd = fftwf_plan_dft_r2c_1d(TD_SIZE, td, fd, FFTW_DESTROY_INPUT | FFTW_MEASURE);
  fd_to_td = fftwf_plan_dft_c2r_1d(TD_SIZE, fd, td, FFTW_DESTROY_INPUT | FFTW_MEASURE);
  
  fftwf_free(td);
  fftwf_free(fd);
}

BlockConvolver::Context::~Context()
{
  fftwf_destroy_plan(td_to_fd);
  fftwf_destroy_plan(fd_to_td);
}

BlockConvolver::Filter::Filter(Context *context, size_t block_size, size_t filter_length, float *coefficients)
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

template <typename T>
BlockConvolver::Buffer<T>::Buffer(size_t len)
  :data((T*)fftwf_malloc(len * sizeof(T)), &fftwf_free)
  ,len(len)
  ,zero(true)
{
  memset(data.get(), 0, len * sizeof(T));
}

template <typename T>
T *BlockConvolver::Buffer<T>::read_ptr() {
  return data.get();
}

template <typename T>
T *BlockConvolver::Buffer<T>::write_ptr() {
  zero = false;
  return data.get();
}

template <typename T>
void BlockConvolver::Buffer<T>::clear() {
  zero = true;
  memset(data.get(), 0, len * sizeof(T));
}

#define SPECTRA_IDX(i) ((spectra_ofs + (i)) % num_blocks)
#define FILTER_IDX(i) ((filter_ofs + (i)) % (num_blocks + 1))

BlockConvolver::BlockConvolver(Context *context, size_t block_size, size_t num_blocks)
  :context(context)
  ,block_size(block_size)
  ,num_blocks(num_blocks)
  ,filter_queue(num_blocks + 1, NULL)
  ,filter_ofs(0)
  ,spectra_ofs(0)
  ,last_tail(block_size)
  ,current_td_old(TD_SIZE)
  ,current_td_new(TD_SIZE)
  ,multiply_out(FD_SIZE)
  ,out_td(TD_SIZE)
{
  assert(context->block_size == block_size);
  
  for (size_t i = 0; i < num_blocks; i++)
  {
    spectra_queue_old.emplace_back(FD_SIZE);
    spectra_queue_new.emplace_back(FD_SIZE);
  }
}

void BlockConvolver::crossfade_filter(const Filter *filter)
{
  if (filter != NULL)
  {
    assert(filter->block_size == block_size);
    assert(filter->num_blocks() <= num_blocks);
  }
  
  filter_queue[FILTER_IDX(0)] = filter;
}

void BlockConvolver::set_filter(const Filter *filter)
{
  if (filter != NULL)
  {
    assert(filter->block_size == block_size);
    assert(filter->num_blocks() <= num_blocks);
  }
  
  for (size_t i = 0; i < filter_queue.size(); i++)
    filter_queue[i] = filter;
}

/** Complex multiply and sum; C++ version.
 *  implements out[i] += a[i] * b[i] for 0 <= i < n */
void complex_mul_sum_cpp(fftwf_complex *out, fftwf_complex *a, fftwf_complex *b, size_t n)
{
  for (size_t i = 0; i < n; i++) {
    // implement out[i] += a[i] * b[i]:
    out[i][0] += a[i][0] * b[i][0] - a[i][1] * b[i][1];
    out[i][1] += a[i][0] * b[i][1] + a[i][1] * b[i][0];
  }
}

#ifdef __SSE3__
/** Complex multiplication of a and b using SSE */
inline __m128 complex_mul_pair_ps(__m128 a, __m128 b)
{
  __m128 a_r = _mm_moveldup_ps(a);
  __m128 a_i = _mm_movehdup_ps(a);
  __m128 tmp1 = _mm_mul_ps(a_r, b);
  __m128 shufd = _mm_shuffle_ps(b, b, _MM_SHUFFLE(2,3,0,1));
  __m128 tmp2 = _mm_mul_ps(a_i, shufd);
  return _mm_addsub_ps(tmp1, tmp2);
}

/** Complex multiply and sum of two pairs of complex numbers using SSE.
 *  implements:
 *    out[0] += a[0] * b[0];
 *    out[1] += a[1] * b[1];
 */
inline void complex_mul_pair_accumulate(fftwf_complex *out, fftwf_complex *a, fftwf_complex *b)
{
  __m128 av = _mm_load_ps((float *)a);
  __m128 bv = _mm_load_ps((float *)b);
  __m128 res = complex_mul_pair_ps(av, bv);

  __m128 out_old = _mm_load_ps((float *)out);
  __m128 out_new = _mm_add_ps(out_old, res);
  _mm_store_ps((float *)out, out_new);
}

/** Complex multiply and sum; SSE version.
 *  implements out[i] += a[i] * b[i] for 0 <= i < n */
void complex_mul_sum_sse(fftwf_complex *out, fftwf_complex *a, fftwf_complex *b, size_t n)
{
  assert(2 * sizeof(fftwf_complex) == 4 * sizeof(float));
  
  // Manually unroll this 3 times; this is the fastest on my machine, and
  // significantly faster than not unrolling.
  size_t offset;
  for (offset = 0; offset + 6 <= n; offset += 6) {
    complex_mul_pair_accumulate(out + offset, a + offset, b + offset);
    complex_mul_pair_accumulate(out + offset + 2, a + offset + 2, b + offset + 2);
    complex_mul_pair_accumulate(out + offset + 4, a + offset + 4, b + offset + 4);
  }
  
  complex_mul_sum_cpp(out + offset, a + offset, b + offset, n - offset);
}
#endif

/** Complex multiply and sum; this uses SSE if available.
 *  implements out[i] += a[i] * b[i] for 0 <= i < n */
void complex_mul_sum(fftwf_complex *out, fftwf_complex *a, fftwf_complex *b, size_t n)
{
#ifdef __SSE3__
  complex_mul_sum_sse(out, a, b, n);
#else
  complex_mul_sum_cpp(out, a, b, n);
#endif
}

/** Mix b into a.
 *   implements a[i] += b[i] for 0 <= i < n */
void mix_into(float *a, float *b, size_t n)
{
  for (size_t i = 0; i < n; i++)
    a[i] += b[i];
}

/** Produce two versions of the block of n samples in in, one faded down across
 * the block, and the other faded up across the block.
 */
void fade_down_and_up(float *in, float *down, float *up, size_t n)
{
  // pull divide out of the loop
  float i_scale = 1.0 / (n-1);
  
  for (size_t i = 0; i < n; i++)
  {
    float a_v = (float)i * i_scale;
    float b_v = 1.0f - a_v;
    up[i]   = a_v * in[i];
    down[i] = b_v * in[i];
  }
}

/** Normalise the fft output. For a block size of n, this divides n samples by
 * 2n (the size of fft used).
 */
void output_norm(float *out, float *in, size_t n)
{
  float norm = 1.0 / (2 * n);
  for (size_t i = 0; i < n; i++)
  {
    out[i] = in[i] * norm;
  }
}

/** Is a buffer null or contain all zeros? */
bool all_zeros(float *in, size_t len) {
  if (in == NULL)
    return true;
  
  for (size_t i = 0; i < len; i++)
    if (in[i] != 0.0f)
      return false;
  
  return true;
}

void BlockConvolver::filter_block(float *in, float *out)
{
  // Pad and fft in into spectra_queue_old and spectra_queue_new, fading if necessary.
  if (all_zeros(in, block_size)) {
    // zero if input is zero
    spectra_queue_old[SPECTRA_IDX(0)].clear();
    spectra_queue_new[SPECTRA_IDX(0)].clear();
  } else {
    // has the filter changed?
    if (filter_queue[FILTER_IDX(1)] != filter_queue[FILTER_IDX(0)]) {
      // if so, produce a version fading down in current_td_old and up in current_td_new, then fft both
      fade_down_and_up(in, current_td_old.write_ptr(), current_td_new.write_ptr(), block_size);
      fftwf_execute_dft_r2c(context->td_to_fd, current_td_old.read_ptr(), spectra_queue_old[SPECTRA_IDX(0)].write_ptr());
      fftwf_execute_dft_r2c(context->td_to_fd, current_td_new.read_ptr(), spectra_queue_new[SPECTRA_IDX(0)].write_ptr());
    } else {
      // otherwise just fft directly to new spectra
      memcpy(current_td_new.write_ptr(), in, block_size * sizeof(float));
      fftwf_execute_dft_r2c(context->td_to_fd, current_td_new.read_ptr(), spectra_queue_new[SPECTRA_IDX(0)].write_ptr());
      spectra_queue_old[SPECTRA_IDX(0)].clear();
    }
  }
  
  // Multiply the spectra of all filters with the corresponding input spectra,
  // summing into multiply_out. Blocks in spectra_queue_new are multiplied with
  // the current filter block, and spectra_queue_old with the previous filter
  // block.
  multiply_out.clear();
  for (size_t i = 0; i < num_blocks; i++)
  {
    const Filter *old_filter = filter_queue[FILTER_IDX(i+1)];
    const Filter *new_filter = filter_queue[FILTER_IDX(i  )];
    
    if (old_filter != NULL && i < old_filter->blocks.size() && !spectra_queue_old[SPECTRA_IDX(i)].zero)
      complex_mul_sum(
          multiply_out.write_ptr(),
          old_filter->blocks[i].get(),
          spectra_queue_old[SPECTRA_IDX(i)].read_ptr(),
          FD_SIZE);
    if (new_filter != NULL && i < new_filter->blocks.size() && !spectra_queue_new[SPECTRA_IDX(i)].zero)
      complex_mul_sum(
          multiply_out.write_ptr(),
          new_filter->blocks[i].get(),
          spectra_queue_new[SPECTRA_IDX(i)].read_ptr(),
          FD_SIZE);
  }
  
  // Inverse fft, then send the first half plus the last tail to the output,
  // and write the second half to last_tail, to be used in the next block.
  if (!multiply_out.zero) {
    fftwf_execute_dft_c2r(context->fd_to_td, multiply_out.read_ptr(), out_td.write_ptr());
    
    // Mix last_tail into the first half out_td, and replace last_tail with the second half of out_td.
    if (!last_tail.zero)
      mix_into(out_td.write_ptr(), last_tail.read_ptr(), block_size);
    memcpy(last_tail.write_ptr(), out_td.read_ptr() + block_size, block_size * sizeof(float));
    
    output_norm(out, out_td.read_ptr(), block_size);
  } else if (!last_tail.zero) {
    // no spectra, just send the last tail if nonzero
    output_norm(out, last_tail.read_ptr(), block_size);
    last_tail.clear();
  } else {
    // no spectra or last tail, zero the output
    memset(out, 0, block_size * sizeof(float));
  }
  
  // older blocks are at higher indices, so move spectra_ofs and filter_ofs left, with wrap-around
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

BBC_AUDIOTOOLBOX_END
