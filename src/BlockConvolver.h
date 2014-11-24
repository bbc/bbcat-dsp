#ifndef __BLOCK_CONVOLVER__
#define __BLOCK_CONVOLVER__

#include <fftw3.h>
#include <vector>
#include <memory>
#include <bbcat-base/misc.h>

BBC_AUDIOTOOLBOX_START

/** BlockConvolver applies overlap-add partitioned convolution to fixed size
 * blocks of samples. */
class BlockConvolver {
  public:
    /** Static data required to perform convolution of a particular block size;
     * may be shared between any number of BlockConvolver instances. */
    class Context {
      public:
        /** Create a Context with a given block size.
         * @param block_size Block size in samples.
         */
        Context(size_t block_size);
        ~Context();
      
      private:
        size_t block_size;
        fftwf_plan td_to_fd;
        fftwf_plan fd_to_td;
      
      friend class BlockConvolver;
      friend class Filter;
    };

    /** A filter response which may be shared between many BlockConvolver instances.
     * 
     * This stores the pre-transformed filter blocks.
     */
    class Filter {
      public:
        /** Create a new Filter given the block size and coefficients.
         * @param context Context required for transformations; only needed
         *                during construction; must have the same block size.
         * @param block_size Block size in samples.
         * @param filter_length Length of coefficients in samples.
         * @param coefficients Filter coefficients.
         */
        Filter(Context *context, size_t block_size, size_t filter_length, float *coefficients);
        /** The number of blocks in the filter. */
        size_t num_blocks() const { return blocks.size(); }

      private:
        size_t block_size;
        std::vector<std::unique_ptr<fftwf_complex, void (*)(void*)> > blocks;
        
      friend class BlockConvolver;
    };

    /** Create a BlockConvolver given the block size and number of blocks.
     * @param context Context required for transformations. This must be alive
     *                for at least as long as this BlockConvolver, and must
     *                have the same block size.
     * @param block_size Block size in samples.
     * @param num_blocks Maximum number of blocks of any filter used.
     */
    BlockConvolver(Context *context, size_t block_size, size_t num_blocks);
    
    /** Pass a block of audio through the filter.
     * @param in Input samples; block_size samples will be read.
     * @param out Output samples; block_size samples will be written.
     */
    void filter_block(float *in, float *out);
    
    /** Crossfade to a new filter during the next block.
     * 
     * This is equivalent to:
     * - Creating a new convolver.
     * - Passing the next block of samples through the old and new convolvers,
     *   with the input to the old faded down across the block, and the input
     *   to the new faded up across the block. All subsequent blocks are passed
     *   through the new filter.
     * - Mixing the output of the old and new filters for the next num_blocks blocks.
     * 
     * @param filter Filter to crossfade to; should be alive for as long as it is active.
     */
    void crossfade_filter(const Filter *filter);
    
    /** Switch to a different filter at the start of the next block.
     * @param filter Filter to switch to; should be alive for as long as it is active.
     */
    void set_filter(const Filter *filter);
    
  private:
    /** Buffer which knows if it contains only zeros.
     */
    template <typename T>
    struct Buffer {
      explicit Buffer(size_t len);
      std::unique_ptr<T, void (*)(void*)> data;
      size_t len;
      bool zero;
      
      /** Get a pointer to the data for reading. */
      T *read_ptr();
      /** Get a pointer to the data for writing; clears the zero flag. */
      T *write_ptr();
      /** Zero the buffer and set the zero flag. */
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
    std::vector<Buffer<fftwf_complex> > spectra_queue_old;
    std::vector<Buffer<fftwf_complex> > spectra_queue_new;
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
