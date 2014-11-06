#include <Convolve.h>
#include <assert.h>
#include <string.h>

ConvolverContext::ConvolverContext(size_t block_size)
  :block_size(block_size)
{
  int dft_td_size = block_size * 2;
  int dft_fd_size = dft_td_size/2+1; // complex
  
  float *td = (float *)fftwf_malloc(dft_td_size * sizeof(float));
  fftwf_complex *fd = (fftwf_complex *)fftwf_malloc(dft_fd_size * sizeof(fftwf_complex));
  
  td_to_fd = fftwf_plan_dft_r2c_1d(dft_td_size, td, fd, FFTW_DESTROY_INPUT | FFTW_MEASURE);
  fd_to_td = fftwf_plan_dft_c2r_1d(dft_td_size, fd, td, FFTW_DESTROY_INPUT | FFTW_MEASURE);
  
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
  int dft_td_size = block_size * 2;
  int dft_fd_size = dft_td_size/2+1; // complex

  float *td = (float *)fftwf_malloc(dft_td_size * sizeof(*td));
  
  for (size_t offset = 0; offset < filter_length; offset += block_size)
  {
    size_t this_block_size = offset + block_size < filter_length ? block_size : filter_length - offset;
    
    // copy to the first half (or less) of td, and zero the second half
    memcpy(td, coefficients + offset, this_block_size * sizeof(*coefficients));
    memset(td + this_block_size, 0, (dft_td_size - this_block_size) * sizeof(*coefficients));
    
    // fft into fd
    fftwf_complex *fd = (fftwf_complex *)fftwf_malloc(dft_fd_size * sizeof(*fd));
    fftwf_execute_dft_r2c(context->td_to_fd, td, fd);
    
    blocks.emplace_back(
        std::unique_ptr<fftwf_complex, void (*)(void*)>(fd, &fftwf_free)
    );
  }
  
  free(td);
}

BlockConvolver::BlockConvolver(ConvolverContext *context, size_t block_size, size_t num_blocks)
  :context(context)
  ,block_size(block_size)
  ,num_blocks(num_blocks)
  ,filter_queue(num_blocks + 1, NULL)
  ,current_td((float *)fftwf_malloc(block_size * 2 * sizeof(float)), &fftwf_free)
  ,multiply_out_a((fftwf_complex *)fftwf_malloc((block_size + 1) * sizeof(fftwf_complex)), &fftwf_free)
  ,multiply_out_b((fftwf_complex *)fftwf_malloc((block_size + 1) * sizeof(fftwf_complex)), &fftwf_free)
  ,out_td_a((float *)fftwf_malloc((block_size * 2) * sizeof(float)), &fftwf_free)
  ,out_td_b((float *)fftwf_malloc((block_size * 2) * sizeof(float)), &fftwf_free)
  ,queue_ofs(0)
  ,filter_ofs(0)
{
  assert(context->block_size == block_size);
  
  int dft_td_size = block_size * 2;
  int dft_fd_size = dft_td_size/2+1; // complex
  
  for (size_t i = 0; i < num_blocks; i++)
  {
    spectra_queue.emplace_back(
        std::unique_ptr<fftwf_complex, void (*)(void*)>(
          (fftwf_complex *)fftwf_malloc(dft_fd_size * sizeof(fftwf_complex)),
          &fftwf_free
        )
    );
    
    memset(spectra_queue[i].get(), 0, dft_fd_size * sizeof(fftwf_complex));
  }
  
  memset(current_td.get(), 0, block_size * 2 * sizeof(float));
}

void BlockConvolver::crossfade_filter(Filter *filter)
{
  assert(filter->block_size == block_size);
  filter_queue[filter_ofs] = filter;
}

void BlockConvolver::set_filter(Filter *filter)
{
  assert(filter->block_size == block_size);
  for (size_t i = 0; i < filter_queue.size(); i++)
    filter_queue[i] = filter;
}

void complex_mul_sum(fftwf_complex *out, fftwf_complex *a, fftwf_complex *b, size_t n)
{
  for (size_t i = 0; i < n; i++) {
    // implement out[i] += a[i] * b[i]:
    out[i][0] += a[i][0] * b[i][0] - a[i][1] * b[i][1];
    out[i][1] += a[i][0] * b[i][1] + a[i][1] * b[i][0];
  }
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
  // fft to spectra_queue[queue_ofs]
  fftwf_execute_dft_r2c(context->td_to_fd, current_td.get(), spectra_queue[queue_ofs].get());
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
    Filter *old_filter = filter_queue[(filter_ofs + i + 1) % filter_queue.size()];
    Filter *new_filter = filter_queue[(filter_ofs + i    ) % filter_queue.size()];
    
    if (old_filter == new_filter) {
      complex_mul_sum(
          multiply_out_a.get(),
          new_filter->blocks[i].get(),
          spectra_queue[(queue_ofs + i) % num_blocks].get(),
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
      Filter *old_filter = filter_queue[(filter_ofs + i + 1) % filter_queue.size()];
      Filter *new_filter = filter_queue[(filter_ofs + i    ) % filter_queue.size()];
      if (old_filter != new_filter) {
        complex_mul_sum(
            multiply_out_a.get(),
            old_filter->blocks[i].get(),
            spectra_queue[(queue_ofs + i) % num_blocks].get(),
            block_size + 1);
        complex_mul_sum(
            multiply_out_b.get(),
            new_filter->blocks[i].get(),
            spectra_queue[(queue_ofs + i) % num_blocks].get(),
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
  
  // older blocks are at higher indices, so move queue_ofs left, with wrap-around.
  queue_ofs--;
  if (queue_ofs < 0) queue_ofs += num_blocks;
  filter_ofs--;
  if (filter_ofs < 0) filter_ofs += filter_queue.size();
  
  // By default the next filter to use is the previous.
  filter_queue[filter_ofs % filter_queue.size()] = filter_queue[(filter_ofs + 1) % filter_queue.size()];
}
