#ifndef __MULTI_LAYER_BUFFER__
#define __MULTI_LAYER_BUFFER__

#include <stdlib.h>
#include <string.h>

#include <vector>

#include "SoundFormatConversions.h"
#include "SoundMixing.h"

BBC_AUDIOTOOLBOX_START

/*--------------------------------------------------------------------------------*/
/** Templated class to handle buffer containing multiple layers of mixed data
 *
 * It can be used where there are multiple contributions to an output (for example,
 * the output of multiple renderers being mixed together) but when different contributions
 * arrive in different block sizes
 *
 * e.g.:
 *
 * renderer 1 outputs arbitrary block sizes because it operates sample by sample
 * renderer 2 outputs in blocks of 256 samples because it is using convolution and a partition size of 256 samples  
 * renderer 3 outputs in blocks of 512 samples because it is using convolution and a partition size of 512 samples  
 *
 * Normally it would take considerable management to work out how many samples can be output because all 3 renderers
 * have contributed to it
 *
 * This class does that management.  At some point, the buffer may look like this:
 *
 * * = frame written
 * - = frame not yet written
 *
 * |***********************--------------------------------------- layer 0
 * |************************************-------------------------- layer 1
 * |*****--------------------------------------------------------- layer 2
 * |*****************************--------------------------------- layer 3
 * |     ^ minposition                  ^ maxposition
 * |^^^^^ this data can be read and removed because it has been written by all layers
 *
 * @note the buffer expands as necessary with incoming data
 */
/*--------------------------------------------------------------------------------*/
template<typename T>
class MultilayerBuffer
{
public:
  /*--------------------------------------------------------------------------------*/
  /** Constructor
   *
   * @param _channels number of channels in buffer
   * @param _layers number of layers
   */
  /*--------------------------------------------------------------------------------*/
  MultilayerBuffer(uint_t _channels = 0, uint_t _layers = 0) {Setup(_channels, _layers);}

  /*--------------------------------------------------------------------------------*/
  /** Copy constructor
   */
  /*--------------------------------------------------------------------------------*/
  MultilayerBuffer(const MultilayerBuffer& obj) : buffer(obj.buffer),
                                                  positions(obj.positions),
                                                  channels(obj.channels),
                                                  minposition(obj.minposition),
                                                  maxposition(obj.maxposition) {}
  ~MultilayerBuffer() {}

  /*--------------------------------------------------------------------------------*/
  /** Assignment operator
   */
  /*--------------------------------------------------------------------------------*/
  MultilayerBuffer& operator = (const MultilayerBuffer& obj)
  {
    buffer      = obj.buffer;
    positions   = obj.positions;
    channels    = obj.channels;
    minposition = obj.minposition;
    maxposition = obj.maxposition;
  }
  
  /*--------------------------------------------------------------------------------*/
  /** Setup / reset
   *
   * @param _channels number of channels in buffer
   * @param _layers number of layers
   */
  /*--------------------------------------------------------------------------------*/
  void Setup(uint_t _channels, uint_t _layers)
  {
    channels = _channels;
    positions.resize(_layers);
    Reset();
  }

  /*--------------------------------------------------------------------------------*/
  /** Reset buffer and positions
   */
  /*--------------------------------------------------------------------------------*/
  void Reset()
  {
    if (positions.size()) memset(&positions[0], 0, positions.size() * sizeof(positions[0]));
    if (buffer.size())    memset(&buffer[0], 0, buffer.size() * sizeof(buffer[0]));
    minposition = maxposition = 0;
  }

  /*--------------------------------------------------------------------------------*/
  /** Change number of layers
   *
   * @param _layers number of layers
   *
   * @note increasing the number of layers resets the available data back to zero but *no* data is lost
   */
  /*--------------------------------------------------------------------------------*/
  void SetLayers(uint_t _layers)
  {
    // if layers are to be *added*, minimum position has to be reset to 0
    if (_layers > (uint_t)positions.size()) minposition = 0;
    positions.resize(_layers);
  }

  /*--------------------------------------------------------------------------------*/
  /** Add layer at end
   */
  /*--------------------------------------------------------------------------------*/
  uint_t AddLayer() {uint_t layer = GetLayers(); SetLayers(layer + 1); return layer;}

  /*--------------------------------------------------------------------------------*/
  /** Delete layer - positions are moved accordingly
   */
  /*--------------------------------------------------------------------------------*/
  void DeleteLayer(uint_t layer)
  {
    if (layer < positions.size())
    {
      // remove layer from array of positions
      positions.erase(positions.begin() + layer);
    }
  }
   
  /*--------------------------------------------------------------------------------*/
  /** Return number of channels in buffer
   */
  /*--------------------------------------------------------------------------------*/
  uint_t GetChannels() const {return channels;}

  /*--------------------------------------------------------------------------------*/
  /** Return number of layers in buffer
   */
  /*--------------------------------------------------------------------------------*/
  uint_t GetLayers()   const {return (uint_t)positions.size();}
  
  /*--------------------------------------------------------------------------------*/
  /** Ensure there is enough space in buffer for new data for a layer
   *
   * @note it is not necessary to call this function if using WriteLayer() below and is only
   * @note needed when GetWritableLayer() is used
   */
  /*--------------------------------------------------------------------------------*/
  void ReserveSpace(uint_t layer, uint_t nframes)
  {
    if (layer < (uint_t)positions.size())
    {
      uint_t maxframes = (positions[layer] + nframes) * channels;
      if (maxframes > buffer.size()) buffer.resize(maxframes);
    }
  }
  
  /*--------------------------------------------------------------------------------*/
  /** Write data into buffer for the specified layer
   *
   * @param layer layer number
   * @param src source buffer
   * @param srcchannel channel into source buffer to read from
   * @param nsrcchannels number of channels in source buffer (in total)
   * @param ndstchannel channel into destination buffer to write to
   * @param nchannels number of channels to write
   * @param nframes number of frames to write
   * @param inc true to increment layer position after writing
   *
   * @note the buffer may be expanded to ensure it can hold the data
   * @note if inc = false, the caller is responsible to ensuring the layer's position is correctly updated (use LayerWritten() below)
   */
  /*--------------------------------------------------------------------------------*/
  void WriteLayer(uint_t layer, const T *src, uint_t srcchannel, uint_t nsrcchannels, uint_t ndstchannel, uint_t nchannels, uint_t nframes, bool inc = true)
  {
    if (layer < (uint_t)positions.size())
    {
      // ensure there is space in the buffer, expanding it if necessary
      ReserveSpace(layer, nframes);

      // mix data into buffer at correct place
      MixSamples(src, srcchannel, nsrcchannels,
                 &buffer[positions[layer] * channels],
                 ndstchannel, channels,
                 nchannels,
                 nframes);

      // increment layer's position if inc is true
      if (inc) LayerWritten(layer, nframes);
    }
  }

  /*--------------------------------------------------------------------------------*/
  /** Return raw pointer to next writable location for a layer
   *
   * @param layer layer number
   *
   * @return pointer to data location
   *
   * @note IMPORTANT: ensure there is space in the buffer by using ReserveSpace() above *before* calling this function!
   */
  /*--------------------------------------------------------------------------------*/
  T *GetWritableLayer(uint_t layer) {return (layer < (uint_t)positions.size()) ? &buffer[positions[layer] * channels] : NULL;}
  
  /*--------------------------------------------------------------------------------*/
  /** Advance layer position by specified amount
   *
   * @param layer layer number
   * @param nframes number of frames to advance the position by
   *
   * @return number of frames available to be read out
   *
   * @note this function *may* result in data being available to be read out, see GetFramesAvailable() below
   */
  /*--------------------------------------------------------------------------------*/
  uint_t LayerWritten(uint_t layer, uint_t nframes)
  {
    if (layer < (uint_t)positions.size())
    {
      // ensure there is enough space in buffer for advancement
      ReserveSpace(layer, nframes);

      // advance position for layer
      positions[layer] += nframes;

      // find minimum position for layers (i.e. the maximum number of frames written by ALL layers)
      uint_t i;
      for (i = 0; i < (uint_t)positions.size(); i++)
      {
        if (i == 0) minposition = positions[i];
        else        minposition = std::min(minposition, positions[i]);
      }

      // find maximum position for layers (i.e. the furthest point written by ANY layer)
      maxposition = std::max(maxposition, positions[layer]);
    }

    return minposition;
  }
  
  /*--------------------------------------------------------------------------------*/
  /** Return number of frames available for reading (that have been written by all layers)
   */
  /*--------------------------------------------------------------------------------*/
  uint_t GetAvailableFrames() const {return minposition;}

  /*--------------------------------------------------------------------------------*/
  /** Get readable buffer and number of frames available
   */
  /*--------------------------------------------------------------------------------*/
  const T *GetReadableBuffer() const {return &buffer[0];}
  const T *GetReadableBuffer(uint_t& availableframes) const {availableframes = minposition; return &buffer[0];}

  /*--------------------------------------------------------------------------------*/
  /** Read data from buffer
   *
   * @param srcchannel source channel within buffer
   * @param dst destination buffer to write to
   * @param dstchannel channel within destination buffer to write to
   * @param ndstchannels total number of channels in destination
   * @param nframes number of frames to read (will be limited)
   * @param inc true to remove read data from buffer after reading
   * @param overwrite true to overwrite data in destination, false to add data to destination
   *
   * @return number of frames actually read
   *
   * @note the supplied number of frames to read will be limited to the number of frames available so large values are acceptable
   */
  /*--------------------------------------------------------------------------------*/
  uint_t ReadBuffer(uint_t srcchannel, T *dst, uint_t dstchannel, uint_t ndstchannels, uint_t nchannels, uint_t nframes, bool inc = true, bool overwrite = true)
  {
    // limit frames to read to number of frames available
    if ((nframes = std::min(nframes, minposition)) > 0)
    {
      if (overwrite)
      {
        // if overwriting, use TransferSamples()
        TransferSamples(&buffer[0], srcchannel, channels,
                        dst,        dstchannel, ndstchannels,
                        nchannels,
                        nframes);
      }
      else
      {
        // if overwriting, use MixSamples()
        MixSamples(&buffer[0], srcchannel, channels,
                   dst,        dstchannel, ndstchannels,
                   nchannels,
                   nframes);
      }

      // if inc is true, remove read data from buffer
      if (inc) BufferRead(nframes);
    }

    return nframes;
  }

  /*--------------------------------------------------------------------------------*/
  /** Read data from buffer
   *
   * @param srcchannel source channel within buffer
   * @param dst destination buffer to write to
   * @param dstchannel channel within destination buffer to write to
   * @param ndstchannels total number of channels in destination
   * @param nframes number of frames to read (will be limited)
   * @param inc true to remove read data from buffer after reading
   * @param overwrite true to overwrite data in destination, false to add data to destination
   *
   * @return number of frames actually read
   *
   * @note the supplied number of frames to read will be limited to the number of frames available so large values are acceptable
   * @note dst will be expanded to fit frames, if necessary
   */
  /*--------------------------------------------------------------------------------*/
  uint_t ReadBuffer(uint_t srcchannel, std::vector<T>& dst, uint_t dstchannel, uint_t ndstchannels, uint_t nchannels, uint_t nframes, bool inc = true, bool overwrite = true)
  {
    // limit frames to read to number of frames available
    if ((nframes = std::min(nframes, minposition)) > 0)
    {
      uint_t samples = ndstchannels * nframes;

      // expand destination buffer if necessary
      if (samples > dst.size()) dst.resize(samples);

      ReadBuffer(srcchannel, &dst[0], dstchannel, ndstchannels, nchannels, nframes, inc, overwrite);
    }

    return nframes;
  }

  /*--------------------------------------------------------------------------------*/
  /** Read data from buffer
   *
   * @param dst destination buffer to write to
   * @param nframes number of frames to read (will be limited)
   * @param inc true to remove read data from buffer after reading
   * @param overwrite true to overwrite data in destination, false to add data to destination
   *
   * @return number of frames actually read
   *
   * @note the supplied number of frames to read will be limited to the number of frames available so large values are acceptable
   * @note all internal channels are transfered
   * @note dst is *assumed* to be of the same number of channels as this object
   * @note dst will be expanded to fit frames, if necessary
   */
  /*--------------------------------------------------------------------------------*/
  uint_t ReadBuffer(std::vector<T>& dst, uint_t nframes, bool inc = true, bool overwrite = true)
  {
    // limit frames to read to number of frames available
    if ((nframes = std::min(nframes, minposition)) > 0)
    {
      uint_t samples = channels * nframes;

      // expand destination buffer if necessary
      if (samples > dst.size()) dst.resize(samples);

      ReadBuffer(0, &dst[0], 0, channels, channels, nframes, inc, overwrite);
    }

    return nframes;
  }
  
  /*--------------------------------------------------------------------------------*/
  /** Remove data from buffer (after it has been read)
   *
   * @param nframes number of frames of data to remove (limited to the number of frames available)
   *
   * @return number of frames actually removed
   */
  /*--------------------------------------------------------------------------------*/
  uint_t BufferRead(uint_t nframes)
  {
    // limit frames to remove to number of frames available
    if ((nframes = std::min(nframes, minposition)) > 0)
    {
      // move min and max positions back
      minposition -= nframes;
      maxposition -= nframes;

      // move each layer's position back
      uint_t i;
      for (i = 0; i < (uint_t)positions.size(); i++)
      {
        positions[i] -= nframes;
      }

      // shift data to front of buffer
      if (maxposition) memmove(&buffer[0], &buffer[nframes * channels], maxposition * channels * sizeof(buffer[0]));
      // and clear newly made space (very important!)
      memset(&buffer[maxposition * channels], 0, nframes * channels * sizeof(buffer[0]));
    }

    // return number of frames removed
    return nframes;
  }

  /*--------------------------------------------------------------------------------*/
  /** Return internal data buffer
   */
  /*--------------------------------------------------------------------------------*/
  const std::vector<T>& GetBuffer() const {return buffer;}

  /*--------------------------------------------------------------------------------*/
  /** Return internal layer positions buffer
   */
  /*--------------------------------------------------------------------------------*/
  const std::vector<uint_t>& GetLayerPositions() const {return positions;}
  
protected:
  std::vector<T>      buffer;
  std::vector<uint_t> positions;
  uint_t              channels;
  uint_t              minposition;
  uint_t              maxposition;
};

BBC_AUDIOTOOLBOX_END

#endif
