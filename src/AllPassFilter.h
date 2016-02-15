#ifndef __ALL_PASS_FILTER__
#define __ALL_PASS_FILTER__

#include "RingBuffer.h"

BBC_AUDIOTOOLBOX_START

/*--------------------------------------------------------------------------------*/
/** All pass filter implementation
 *
 * Implements y[n] = c.x[n] + x[n - d] - c.y[n - d]
 *
 * For multiple channels
 */
/*--------------------------------------------------------------------------------*/
template<typename TYPE>
class AllPassFilter
{
public:
  /*--------------------------------------------------------------------------------*/
  /** Constructor specifying channels and samples
   *
   * @param _nchannels number of channels
   * @param _delay delay (per-channel)
   */
  /*--------------------------------------------------------------------------------*/
  AllPassFilter(uint_t _nchannels = 0, uint_t _delay = 0) : nchannels(_nchannels),
                                                            delay(_delay),
                                                            buffer(nchannels * delay),
                                                            coeff(TYPE()) {}
  /*--------------------------------------------------------------------------------*/
  /** Copy constructor
   */
  /*--------------------------------------------------------------------------------*/
  AllPassFilter(const AllPassFilter& obj) : nchannels(obj.nchannels),
                                            delay(obj.delay),
                                            buffer(obj.buffer),
                                            coeff(obj.coeff) {}
  ~AllPassFilter() {}

  /*--------------------------------------------------------------------------------*/
  /** Set number of channels
   *
   * @note can be called at any time but will reset the ring buffer
   */
  /*--------------------------------------------------------------------------------*/
  void SetChannels(uint_t n) {nchannels = n; buffer.SetLength(nchannels * delay);}

  /*--------------------------------------------------------------------------------*/
  /** Set delay
   *
   * @note can be called at any time but will reset the ring buffer
   */
  /*--------------------------------------------------------------------------------*/
  void SetDelay(uint_t n) {delay = n; buffer.SetLength(nchannels * delay);}

  /*--------------------------------------------------------------------------------*/
  /** Set coeff
   *
   * @note avoid changing this during usage, will lead to strange effects because of w[n]
   */
  /*--------------------------------------------------------------------------------*/
  void SetCoeff(TYPE c) {coeff = c;}

  /*--------------------------------------------------------------------------------*/
  /** Process a sample of data
   */
  /*--------------------------------------------------------------------------------*/
  TYPE Process(TYPE x)
  {
    TYPE y = coeff * x + buffer.Read();         // y[n] = c.x[n] + w[n - d]
    buffer.Write(x - coeff * y);                // save w[n] = x[n] - c.y[n]
    return y;
  }
  
  /*--------------------------------------------------------------------------------*/
  /** Process audio
   *
   * @param src source buffer
   * @param dst destination buffer
   * @param srcchannel source starting channel
   * @param nsrcchannels total number of source channels
   * @param dstchannel destination starting channel
   * @param ndstchannels total number of destination channels
   * @param nframes number of sample frames to process
   *
   * @note src can == dst IFF src and dst channel parameters are the same 
   */
  /*--------------------------------------------------------------------------------*/
  void Process(const TYPE *src, TYPE *dst, uint_t srcchannel, uint_t nsrcchannels, uint_t dstchannel, uint_t ndstchannels, uint_t nframes = 1)
  {
    uint_t i;

    // offset buffer pointers by starting channels
    src += srcchannel;
    dst += dstchannel;

    if (nchannels == 1)
    {
      // optimized version for single channel
      for (i = 0; i < nframes; i++, src += nsrcchannels, dst += ndstchannels)
      {
        TYPE x = src[0];                            // take copy of src in case dst == src
        dst[0] = coeff * x + buffer.Read();         // y[n] = c.x[n] + w[n - d]
        buffer.Write(x - coeff * dst[0]);           // save w[n] = x[n] - c.y[n]
      }
    }
    else
    {
      uint_t j, n = nchannels;
      
      // calculate number of channels that can be processed
      n = std::min(n, limited::subz(nsrcchannels, srcchannel));
      n = std::min(n, limited::subz(ndstchannels, dstchannel));

      for (i = 0; i < nframes; i++, src += nsrcchannels, dst += ndstchannels)
      {
        for (j = 0; j < n; j++)
        {
          TYPE x = src[j];                          // take copy of src in case dst == src
          dst[j] = coeff * x + buffer.Read();       // y[n] = c.x[n] + w[n - d]
          buffer.Write(x - coeff * dst[j]);         // save w[n] = x[n] - c.y[n]
        }
        // skip over unused channels in ring buffers
        buffer.Advance(nchannels - n);
      }
    }
  }
  
protected:
  uint_t           nchannels;
  uint_t           delay;
  RingBuffer<TYPE> buffer;
  TYPE             coeff;
};

template<typename TYPE>
class AllPassFilterChain
{
public:
  /*--------------------------------------------------------------------------------*/
  /** Constructor for a chain of all-pass filters
   *
   * @param _nchannels number of channels (same across all filters)
   * @param nfilters number of filters
   * @param delays optional array of delays for each filter
   * @param coeffs optional array of coeffs for each filter
   */
  /*--------------------------------------------------------------------------------*/
  AllPassFilterChain(uint_t _nchannels = 0, uint_t nfilters = 0, const uint_t *delays = NULL, const TYPE *coeffs = NULL) : filters(nfilters),
                                                                                                                           nchannels(_nchannels) {
    uint_t i;
    for (i = 0; i < filters.size(); i++)
    {
      filters[i].SetChannels(nchannels);
      if (delays) filters[i].SetDelay(delays[i]);
      if (coeffs) filters[i].SetCoeff(coeffs[i]);
    }
  }
  /*--------------------------------------------------------------------------------*/
  /** Copy constructor
   */
  /*--------------------------------------------------------------------------------*/
  AllPassFilterChain(const AllPassFilterChain& obj) : filters(obj.filters),
                                                      nchannels(obj.nchannels) {}
  ~AllPassFilterChain() {}

  /*--------------------------------------------------------------------------------*/
  /** Set the number of filters and optionally set the delays and coeffs for each filter
   */
  /*--------------------------------------------------------------------------------*/
  void SetFilterCount(uint_t nfilters, const uint_t *delays = NULL, const TYPE *coeffs = NULL)
  {
    uint_t i;

    filters.resize(nfilters);

    for (i = 0; i < filters.size(); i++)
    {
      filters[i].SetChannels(nchannels);
      if (delays) filters[i].SetDelay(delays[i]);
      if (coeffs) filters[i].SetCoeff(coeffs[i]);
    }
  }
  
  /*--------------------------------------------------------------------------------*/
  /** Set number of channels across all filters
   *
   * @note can be called at any time but will reset the filters
   */
  /*--------------------------------------------------------------------------------*/
  void SetChannels(uint_t _nchannels)
  {
    uint_t i;

    nchannels = _nchannels;
    for (i = 0; i < filters.size(); i++)
    {
      filters[i].SetChannels(nchannels);
    }
  }

  /*--------------------------------------------------------------------------------*/
  /** Set delay for a filter
   *
   * @note can be called at any time but will reset the filter
   */
  /*--------------------------------------------------------------------------------*/
  void SetDelay(uint_t n, uint_t delay)
  {
    if (n < filters.size()) filters[n].SetDelay(delay);
  }

  /*--------------------------------------------------------------------------------*/
  /** Set coeff for a filter
   *
   * @note avoid changing this during usage, will lead to strange effects because of w[n]
   */
  /*--------------------------------------------------------------------------------*/
  void SetCoeff(uint_t n, TYPE coeff)
  {
    if (n < filters.size()) filters[n].SetCoeff(coeff);
  }

  /*--------------------------------------------------------------------------------*/
  /** Process audio through chain of filters
   *
   * @param src source buffer
   * @param dst destination buffer
   * @param srcchannel source starting channel
   * @param nsrcchannels total number of source channels
   * @param dstchannel destination starting channel
   * @param ndstchannels total number of destination channels
   * @param nframes number of sample frames to process
   *
   * @note src can == dst IFF src and dst channel parameters are the same 
   */
  /*--------------------------------------------------------------------------------*/
  void Process(const TYPE *src, TYPE *dst, uint_t srcchannel, uint_t nsrcchannels, uint_t dstchannel, uint_t ndstchannels, uint_t nframes = 1)
  {
    uint_t i;

    for (i = 0; i < filters.size(); i++)
    {
      if (i == 1)
      {
        // after first filter, set the source to be the destination, including the channel parameters
        src          = dst;
        srcchannel   = dstchannel;
        nsrcchannels = ndstchannels;
      }

      // process audio
      filters[i].Process(src, dst, srcchannel, nsrcchannels, dstchannel, ndstchannels, nframes);
    }
  }
  
protected:
  std::vector<AllPassFilter<TYPE> > filters;
  uint_t                            nchannels;
};

BBC_AUDIOTOOLBOX_END

#endif
