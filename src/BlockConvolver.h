#ifndef __BLOCK_CONVOLVER__
#define __BLOCK_CONVOLVER__

#include <fftw3.h>
#include <vector>
#include <memory>
#include <bbcat-base/misc.h>

BBC_AUDIOTOOLBOX_START

/** BlockConvolver applies overlap-add partitioned convolution to fixed size
 * blocks of samples. */
class BlockConvolver
{
  public:
    /// Type for real data (float).
    typedef float real_t;
    /// Type for complex data (fftwf_complex).
    typedef fftwf_complex complex_t;
    
    /** Static data required to perform convolution of a particular block size;
     * may be shared between any number of BlockConvolver and
     * BlockConvolver::Filter instances. */
    class Context
    {
      public:
        /** Create a Context with a given block size.
         * @param block_size Block size in samples.
         */
        Context(size_t block_size);
        ~Context();
      
      private:
        // number of samples in a block
        const size_t block_size;
        // time domain size; block_size * 2
        const size_t td_size;
        // frequency domain size; block_size + 1
        const size_t fd_size;
        
        // fft plans between time and frequency domain blocks.
        fftwf_plan td_to_fd;
        fftwf_plan fd_to_td;
      
      friend class BlockConvolver;
      friend class Filter;
    };

    /** A filter response which may be shared between many BlockConvolver instances.
     * 
     * This stores the pre-transformed filter blocks. */
    class Filter
    {
      public:
        /** Create a new Filter given the block size and coefficients.
         * @param ctx Context required for transformations; must have a
         *            lifetime at least as long as the created object.
         * @param filter_length Length of coefficients in samples.
         * @param coefficients Filter coefficients.
         */
        Filter(Context *ctx, size_t filter_length, real_t *coefficients);
        /** The number of blocks in the filter. */
        size_t num_blocks() const { return blocks.size(); }

      private:
        std::vector<std::unique_ptr<complex_t[], void (*)(void*)> > blocks;
        const Context *ctx;
        
      friend class BlockConvolver;
    };

    /** Create a BlockConvolver given the block size and number of blocks.
     * @param ctx Context required for transformations. This must be alive
     *            for at least as long as the created object.
     * @param num_blocks Maximum number of blocks of any filter used.
     */
    BlockConvolver(Context *ctx, size_t num_blocks);
    
    /** Create a BlockConvolver given the block size and number of blocks.
     *  If filter == NULL, num_blocks must be specified.
     * @param ctx Context required for transformations. This must be alive
     *            for at least as long as the created object.
     * @param filter Initial filter to be used, or NULL for no filter.
     * @param num_blocks Maximum number of blocks of any filter used; using 0
     *                   will take the number of blocks from the passed filter.
     */
    BlockConvolver(Context *ctx, const Filter *filter, size_t num_blocks = 0);
    
    /** Pass a block of audio through the filter.
     * @param in Input samples; block_size samples will be read.
     * @param out Output samples; block_size samples will be written.
     */
    void filter_block(const real_t *in, real_t *out);
    
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
     * @param filter Filter to crossfade to; should be alive for as long as it
     *               is active. Pass NULL for no filter.
     */
    void crossfade_filter(const Filter *filter);
    
    /** Switch to a different filter at the start of the next block.
     * @param filter Filter to switch to; should be alive for as long as it is
     *               active. Pass NULL for no filter.
     */
    void set_filter(const Filter *filter);
    
  private:
    /** Buffer which knows if it contains only zeros.
     */
    template <typename T>
    class Buffer
    {
      public:
        explicit Buffer(size_t len);
        std::unique_ptr<T[], void (*)(void*)> data;
        const size_t len;
        bool zero;
        
        /** Get a pointer to the data for reading. */
        const T *read_ptr();
        /** Get a pointer to the data for writing; clears the zero flag. */
        T *write_ptr();
        /** Zero the buffer and set the zero flag. */
        void clear();
    };
    
    const Context *ctx;
    const size_t num_blocks;
    
    // filters(i) accesses a circular buffer of filters, length num_blocks + 1.
    // on each frame, the input is crossfaded up and down, and passed through
    // - filters[1:] (old)
    // - filters[:-1] (new)
    // After filter_block, filters has been shifted one along, with filters(0)
    // left at filters(1); filters(0) is then set to filters(1), such that the
    // next filter is the same as the previous by default.
    // set_filter simply writes to filters(0).
    const Filter *&filters(size_t i);
    // filter_queue and filter_ofs implement the above circular buffer.
    std::vector<const Filter *> filter_queue;
    size_t filter_ofs;
    
    // a queue of spectra of the input, after zero padding on the right hand side.
    // If the filter is changed before the input block i frames ago, spectra_old(i)
    // contains the input faded down (to convolve with the old filter) and
    // spectra_new(i) contains the input faded up (to convolve with the
    // new filter).
    // If the filter is not changed, only spectra_new(i) contains data.
    // num_blocks in length, each of size n+1.
    Buffer<complex_t> &spectra_old(size_t i);
    Buffer<complex_t> &spectra_new(size_t i);
    // implementation the above circular buffer
    std::vector<Buffer<complex_t> > spectra_queue_old;
    std::vector<Buffer<complex_t> > spectra_queue_new;
    size_t spectra_ofs;
    
    // Rotate the filter and spectra queues.
    void rotate_queues();
    
    // The second half of the ifft output for the last block, added to the first half before output; size n.
    Buffer<real_t> last_tail;
    
    // Temporaries used by filter_block to store the padded and faded inputs; size 2n, the second half should always be zeros.
    Buffer<real_t> current_td_old;
    Buffer<real_t> current_td_new;
    // Temporary used by filter_block to store the multiplied spectrum; size n+1.
    Buffer<complex_t> multiply_out;
    // Temporary used to store the time domain output, size 2n.
    Buffer<real_t> out_td;
};

BBC_AUDIOTOOLBOX_END

#endif
