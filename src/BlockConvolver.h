#ifndef __BLOCK_CONVOLVER__
#define __BLOCK_CONVOLVER__

#include <fftw3.h>
#include <vector>
#include <memory>
#include <bbcat-base/misc.h>

BBC_AUDIOTOOLBOX_START

class BlockConvolver {
  public:
    class Context {
      public:
        Context(size_t block_size);
        ~Context();
      private:
        size_t block_size;
        fftwf_plan td_to_fd;
        fftwf_plan fd_to_td;
      
      friend class BlockConvolver;
      friend class Filter;
    };

    class Filter {
      public:
        Filter(Context *context, size_t block_size, size_t filter_length, float *coefficients);
        size_t num_blocks() const { return blocks.size(); }

      private:
        size_t block_size;
        std::vector<std::unique_ptr<fftwf_complex, void (*)(void*)> > blocks;
      
      friend class BlockConvolver;
    };

    BlockConvolver(Context *context, size_t block_size, size_t num_blocks);
    
    void filter_block(float *in, float *out);
    
    void crossfade_filter(const Filter *filter);
    void set_filter(const Filter *filter);
    
  private:
    template <typename T>
    struct Buffer {
      explicit Buffer(size_t len);
      std::unique_ptr<T, void (*)(void*)> data;
      size_t len;
      bool zero = true;
      
      T *read_ptr();
      T *write_ptr();
      void clear();
    };
    
    Context *context;
    size_t block_size;
    size_t num_blocks;
    
    // on each frame, the output is a crossfade between:
    // - filter_queue[1:] (old)
    // - filter_queue[:-1] (new)
    // after filter_block, filter_queue has been shifted one along, with filter_queue[0] left at filter_queue[0]
    // set_filter simply writes to filter_queue[0].
    // this should be num_blocks + 1 in length.
    std::vector<const Filter *> filter_queue;
    int filter_ofs;
    
    // Each of these corresponds to the matching item in filter_queue.
    // should contain the FFT of last + current (i.e. twice the block width)
    // num_blocks in length.
    std::vector<Buffer<fftwf_complex>> spectra_queue;
    int spectra_ofs;
    
    // The current time domain input; length 2n.
    // This is used as a temporary during filter_block, and the first half is
    // used to store the last filter block between calls; the second half is undefined.
    Buffer<float> current_td;
    
    // temporaries used by filter_block to store the multiplied spectrums; size n+1
    Buffer<fftwf_complex> multiply_out_a;
    Buffer<fftwf_complex> multiply_out_b;
    
    // temporaries used to store the time domain outputs before cross fading, size 2n
    Buffer<float> out_td_a;
    Buffer<float> out_td_b;
};

BBC_AUDIOTOOLBOX_END

#endif
