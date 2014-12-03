#include "BlockConvolver.h"
#include <assert.h>
#include <string.h>

#include <fftw3.h>

#ifdef __SSE3__
#  include <xmmintrin.h>
#  include <pmmintrin.h>
#endif

BBC_AUDIOTOOLBOX_START

struct BlockConvolver::Context::fft_data
{
  fftwf_plan plan;
};

// aliases so we don't have to keep using BlockConvolver:: in function headers
typedef BlockConvolver::real_t real_t;
typedef BlockConvolver::complex_t complex_t;

// Helper to allocate a fftw-compatible block of memory in a unique_ptr
template <typename T>
std::unique_ptr<T[], void (*)(void*)> fftw_malloc_unique(size_t len)
{
  std::unique_ptr<T[], void (*)(void*)> data((T*)fftwf_malloc(len * sizeof(T)), &fftwf_free);
  memset(data.get(), 0, len * sizeof(T));
  return data;
}


BlockConvolver::Context::Context(size_t block_size)
  :block_size(block_size)
  ,td_size(block_size * 2)
  ,fd_size(block_size + 1)
  ,td_to_fd(new fft_data)
  ,fd_to_td(new fft_data)
{
  auto td = fftw_malloc_unique<real_t>(td_size);
  auto fd = fftw_malloc_unique<complex_t>(fd_size);
  
  td_to_fd->plan = fftwf_plan_dft_r2c_1d(td_size, td.get(), fd.get(), FFTW_DESTROY_INPUT | FFTW_MEASURE);
  fd_to_td->plan = fftwf_plan_dft_c2r_1d(td_size, fd.get(), td.get(), FFTW_DESTROY_INPUT | FFTW_MEASURE);
}

BlockConvolver::Context::~Context()
{
  fftwf_destroy_plan(td_to_fd->plan);
  fftwf_destroy_plan(fd_to_td->plan);
}


BlockConvolver::Filter::Filter(Context *ctx, size_t filter_length, real_t *coefficients)
  :ctx(ctx)
{
  auto td = fftw_malloc_unique<real_t>(ctx->td_size);
  
  for (size_t offset = 0; offset < filter_length; offset += ctx->block_size)
  {
    size_t this_block_size = offset + ctx->block_size < filter_length ? ctx->block_size : filter_length - offset;
    
    // copy to the first half (or less) of td, and zero the second half
    memcpy(td.get(), coefficients + offset, this_block_size * sizeof(*coefficients));
    memset(td.get() + this_block_size, 0, (ctx->td_size - this_block_size) * sizeof(*coefficients));
    
    // fft into fd
    auto fd = fftw_malloc_unique<complex_t>(ctx->fd_size);
    fftwf_execute_dft_r2c(ctx->td_to_fd->plan, td.get(), fd.get());
    
    blocks.emplace_back(std::move(fd));
  }
}


template <typename T>
BlockConvolver::Buffer<T>::Buffer(size_t len)
  :data(fftw_malloc_unique<T>(len))
  ,len(len)
  ,zero(true)
{
  memset(data.get(), 0, len * sizeof(T));
}

template <typename T>
const T *BlockConvolver::Buffer<T>::read_ptr()
{
  return data.get();
}

template <typename T>
T *BlockConvolver::Buffer<T>::write_ptr()
{
  zero = false;
  return data.get();
}

template <typename T>
void BlockConvolver::Buffer<T>::clear()
{
  if (!zero)
  {
    zero = true;
    memset(data.get(), 0, len * sizeof(T));
  }
}


BlockConvolver::BlockConvolver(Context *ctx, size_t num_blocks)
  :ctx(ctx)
  ,num_blocks(num_blocks)
  ,filter_queue(num_blocks + 1, NULL)
  ,filter_ofs(0)
  ,spectra_ofs(0)
  ,last_tail(ctx->block_size)
  ,current_td_old(ctx->td_size)
  ,current_td_new(ctx->td_size)
  ,multiply_out(ctx->fd_size)
  ,out_td(ctx->td_size)
{
  for (size_t i = 0; i < num_blocks; i++)
  {
    spectra_queue_old.emplace_back(Buffer<complex_t>(ctx->fd_size));
    spectra_queue_new.emplace_back(Buffer<complex_t>(ctx->fd_size));
  }
}

// num_blocks = 0 by default; if so use num_blocks from filter.
BlockConvolver::BlockConvolver(Context *ctx, const Filter *filter, size_t num_blocks)
  :BlockConvolver(ctx, num_blocks > 0 ? num_blocks : filter->num_blocks())
{
  set_filter(filter);
}

void BlockConvolver::crossfade_filter(const Filter *filter)
{
  if (filter != NULL)
  {
    assert(filter->ctx->block_size == ctx->block_size);
    assert(filter->num_blocks() <= num_blocks);
  }
  
  filters(0) = filter; // set filter through reference
}

void BlockConvolver::set_filter(const Filter *filter)
{
  if (filter != NULL)
  {
    assert(filter->ctx->block_size == ctx->block_size);
    assert(filter->num_blocks() <= num_blocks);
  }
  
  for (size_t i = 0; i < filter_queue.size(); i++)
  {
    filter_queue[i] = filter;
  }
}

const BlockConvolver::Filter *&BlockConvolver::filters(size_t i)
{
  return filter_queue[(filter_ofs + i) % (num_blocks + 1)];
}

BlockConvolver::Buffer<complex_t> &BlockConvolver::spectra_old(size_t i)
{
  return spectra_queue_old[(spectra_ofs + i) % num_blocks];
}

BlockConvolver::Buffer<complex_t> &BlockConvolver::spectra_new(size_t i)
{
  return spectra_queue_new[(spectra_ofs + i) % num_blocks];
}

void BlockConvolver::rotate_queues()
{
  // older blocks are at higher indices, so move spectra_ofs and filter_ofs left, with wrap-around
  spectra_ofs = (spectra_ofs + num_blocks - 1) % num_blocks;
  filter_ofs = (filter_ofs + (num_blocks + 1) - 1) % (num_blocks + 1);
  
  // By default the next filter to use is the previous.
  filters(0) = filters(1);
}

/** Complex multiply and sum; C++ version.
 *  implements out[i] += a[i] * b[i] for 0 <= i < n */
void complex_mul_sum_cpp(complex_t *out, const complex_t *a, const complex_t *b, size_t n)
{
  for (size_t i = 0; i < n; i++)
  {
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
inline void complex_mul_pair_accumulate(complex_t *out, const complex_t *a, const complex_t *b)
{
  __m128 av = _mm_load_ps((real_t *)a);
  __m128 bv = _mm_load_ps((real_t *)b);
  __m128 res = complex_mul_pair_ps(av, bv);

  __m128 out_old = _mm_load_ps((real_t *)out);
  __m128 out_new = _mm_add_ps(out_old, res);
  _mm_store_ps((real_t *)out, out_new);
}

/** Complex multiply and sum; SSE version.
 *  implements out[i] += a[i] * b[i] for 0 <= i < n */
void complex_mul_sum_sse(complex_t *out, const complex_t *a, const complex_t *b, size_t n)
{
  assert(2 * sizeof(complex_t) == 4 * sizeof(real_t));
  
  // Manually unroll this 3 times; this is the fastest on my machine, and
  // significantly faster than not unrolling.
  size_t offset;
  for (offset = 0; offset + 6 <= n; offset += 6)
  {
    complex_mul_pair_accumulate(out + offset, a + offset, b + offset);
    complex_mul_pair_accumulate(out + offset + 2, a + offset + 2, b + offset + 2);
    complex_mul_pair_accumulate(out + offset + 4, a + offset + 4, b + offset + 4);
  }
  
  complex_mul_sum_cpp(out + offset, a + offset, b + offset, n - offset);
}
#endif

/** Complex multiply and sum; this uses SSE if available.
 *  implements out[i] += a[i] * b[i] for 0 <= i < n */
void complex_mul_sum(complex_t *out, const complex_t *a, const complex_t *b, size_t n)
{
#ifdef __SSE3__
  complex_mul_sum_sse(out, a, b, n);
#else
  complex_mul_sum_cpp(out, a, b, n);
#endif
}

/** Mix b into a.
 *   implements a[i] += b[i] for 0 <= i < n */
void mix_into(real_t *a, const real_t *b, size_t n)
{
  for (size_t i = 0; i < n; i++) a[i] += b[i];
}

/** Produce two versions of the block of n samples in in, one faded down across
 * the block, and the other faded up across the block.
 */
void fade_down_and_up(const real_t *in, real_t *down, real_t *up, size_t n)
{
  // pull divide out of the loop
  real_t i_scale = 1.0 / n;
  
  for (size_t i = 0; i < n; i++)
  {
    real_t a_v = (real_t)i * i_scale;
    real_t b_v = 1.0f - a_v;
    up[i]   = a_v * in[i];
    down[i] = b_v * in[i];
  }
}

/** Normalise the fft output. For a block size of n, this divides n samples by
 * 2n (the size of fft used).
 */
void output_norm(real_t *out, const real_t *in, size_t n)
{
  real_t norm = 1.0 / (2 * n);
  for (size_t i = 0; i < n; i++)
  {
    out[i] = in[i] * norm;
  }
}

/** Is a buffer null or contain all zeros? */
bool all_zeros(const real_t *in, size_t len) {
  if (in == NULL) return true;
  
  for (size_t i = 0; i < len; i++)
  {
    if (in[i] != 0.0f) return false;
  }
  
  return true;
}

void BlockConvolver::filter_block(const real_t *in, real_t *out)
{
  // Pad and fft in into spectra_queue_old and spectra_queue_new, fading if necessary.
  if (all_zeros(in, ctx->block_size))
  {
    // zero if input is zero
    spectra_old(0).clear();
    spectra_new(0).clear();
  }
  else
  {
    // has the filter changed?
    if (filters(1) != filters(0))
    {
      // if so, produce a version fading down in current_td_old and up in current_td_new, then fft both
      fade_down_and_up(in, current_td_old.write_ptr(), current_td_new.write_ptr(), ctx->block_size);
      // fft, and clear the second half as this may have been modified by fftw
      fftwf_execute_dft_r2c(ctx->td_to_fd->plan, current_td_old.write_ptr(), spectra_old(0).write_ptr());
      memset(current_td_old.write_ptr() + ctx->block_size, 0, ctx->block_size * sizeof(real_t));
      fftwf_execute_dft_r2c(ctx->td_to_fd->plan, current_td_new.write_ptr(), spectra_new(0).write_ptr());
      memset(current_td_new.write_ptr() + ctx->block_size, 0, ctx->block_size * sizeof(real_t));
    }
    else
    {
      // otherwise just fft directly to new spectra.
      memcpy(current_td_new.write_ptr(), in, ctx->block_size * sizeof(real_t));
      // fft, and clear the second half as this may have been modified by fftw
      fftwf_execute_dft_r2c(ctx->td_to_fd->plan, current_td_new.write_ptr(), spectra_new(0).write_ptr());
      memset(current_td_new.write_ptr() + ctx->block_size, 0, ctx->block_size * sizeof(real_t));
      // all in spectra_new for this block
      spectra_old(0).clear();
    }
  }
  
  // Multiply the spectra of all filters with the corresponding input spectra,
  // summing into multiply_out. Blocks in spectra_queue_new are multiplied with
  // the current filter block, and spectra_queue_old with the previous filter
  // block.
  multiply_out.clear();
  for (size_t i = 0; i < num_blocks; i++)
  {
    const Filter *old_filter = filters(i+1);
    const Filter *new_filter = filters(i);
    
    if (old_filter != NULL && i < old_filter->blocks.size() && !spectra_old(i).zero)
    {
      complex_mul_sum(
          multiply_out.write_ptr(),
          old_filter->blocks[i].get(),
          spectra_old(i).read_ptr(),
          ctx->fd_size);
    }
    if (new_filter != NULL && i < new_filter->blocks.size() && !spectra_new(i).zero)
    {
      complex_mul_sum(
          multiply_out.write_ptr(),
          new_filter->blocks[i].get(),
          spectra_new(i).read_ptr(),
          ctx->fd_size);
    }
  }
  
  // Inverse fft, then send the first half plus the last tail to the output,
  // and write the second half to last_tail, to be used in the next block.
  if (!multiply_out.zero)
  {
    fftwf_execute_dft_c2r(ctx->fd_to_td->plan, multiply_out.write_ptr(), out_td.write_ptr());
    
    // Mix last_tail into the first half out_td, and replace last_tail with the second half of out_td.
    if (!last_tail.zero)
      mix_into(out_td.write_ptr(), last_tail.read_ptr(), ctx->block_size);
    memcpy(last_tail.write_ptr(), out_td.read_ptr() + ctx->block_size, ctx->block_size * sizeof(real_t));
    
    output_norm(out, out_td.read_ptr(), ctx->block_size);
  }
  else if (!last_tail.zero)
  {
    // no spectra, just send the last tail if nonzero
    output_norm(out, last_tail.read_ptr(), ctx->block_size);
    last_tail.clear();
  }
  else
  {
    // no spectra or last tail, zero the output
    memset(out, 0, ctx->block_size * sizeof(real_t));
  }
  
  rotate_queues();
}

BBC_AUDIOTOOLBOX_END
