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
    
    // The spectra of the input, after padding on the right hand side.
    // If the filter is changed before input block i, spectra_queue_old[i]
    // contains the input faded down (to convolve with the old filter) and
    // spectra_queue_new[i] contains the input faded up (to convolve with the
    // new filter).
    // If the filter is not changed, only spectra_queue_new contains data.
    // num_blocks in length, each of size n+1.
    std::vector<Buffer<fftwf_complex>> spectra_queue_old;
    std::vector<Buffer<fftwf_complex>> spectra_queue_new;
    int spectra_ofs;
    
    // The second half of the ifft output for the last block, added to the first half before output; size n.
    Buffer<float> last_tail;
    
    // Temporaries used by filter_block to store the padded and faded inputs; size 2n, the second half should always be zeros.
    Buffer<float> current_td_old;
    Buffer<float> current_td_new;
    // Temporary used by filter_block to store the multiplied spectrum; size n+1.
    Buffer<fftwf_complex> multiply_out;
    // Temporary used to store the time domain output, size 2n.
    Buffer<float> out_td;
};

BBC_AUDIOTOOLBOX_END

#endif
