#ifndef __CONVOLVE__
#define __CONVOLVE__

#include <fftw3.h>
#include <vector>
#include <memory>
#include <bbcat-base/misc.h>

BBC_AUDIOTOOLBOX_START

class ConvolverContext {
public:
  ConvolverContext(size_t block_size);
  ~ConvolverContext();
  size_t block_size;
  fftwf_plan td_to_fd;
  fftwf_plan fd_to_td;
};

class Filter {
public:
  Filter(ConvolverContext *context, size_t block_size, size_t filter_length, float *coefficients);
  
  size_t block_size;
  std::vector<std::unique_ptr<fftwf_complex, void (*)(void*)> > blocks;
};

class BlockConvolver {
public:
  void filter_block(float *in, float *out);
  
  void crossfade_filter(Filter *filter);
  void set_filter(Filter *filter);
  
  BlockConvolver(ConvolverContext *context, size_t block_size, size_t num_blocks);
  
private:
  ConvolverContext *context;
  size_t block_size;
  size_t num_blocks;
  
  // on each frame, the output is a crossfade between:
  // - filter_queue[1:] (old)
  // - filter_queue[:-1] (new)
  // after filter_block, filter_queue has been shifted one along, with filter_queue[0] left at filter_queue[0]
  // set_filter simply writes to filter_queue[0].
  // this should be num_blocks + 1 in length.
  std::vector<Filter *> filter_queue;
  int filter_ofs;
  
  // Each of these corresponds to the matching item in filter_queue.
  // should contain the FFT of last + current (i.e. twice the block width)
  // num_blocks in length.
  std::vector<std::unique_ptr<fftwf_complex, void (*)(void*)> > spectra_queue;
  int spectra_ofs;
  
  // The current time domain input; length 2n.
  // This is used as a temporary during filter_block, and the first half is
  // used to store the last filter block between calls; the second half is undefined.
  std::unique_ptr<float, void (*)(void*)> current_td;
  
  // temporaries used by filter_block to store the multiplied spectrums; size n+1
  std::unique_ptr<fftwf_complex, void (*)(void*)> multiply_out_a;
  std::unique_ptr<fftwf_complex, void (*)(void*)> multiply_out_b;
  
  // temporaries used to store the time domain outputs before cross fading, size 2n
  std::unique_ptr<float, void (*)(void*)> out_td_a;
  std::unique_ptr<float, void (*)(void*)> out_td_b;
};

BBC_AUDIOTOOLBOX_END

#endif
