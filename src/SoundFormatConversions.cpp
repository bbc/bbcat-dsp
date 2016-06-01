
#include <math.h>

#define BBCDEBUG_LEVEL 1
#include "SoundFormatConversions.h"
#include "SoundFormatRawConversions.h"

BBC_AUDIOTOOLBOX_START

/*--------------------------------------------------------------------------------*/
/** Array of number of bytes for each format
 */
/*--------------------------------------------------------------------------------*/
const uint8_t SoundFormatBits[SampleFormat_Count] =
{
  1,  // SampleFormat_Unknown
  16, // SampleFormat_16bit
  24, // SampleFormat_24bit
  32, // SampleFormat_32bit
  32, // SampleFormat_Float
  64, // SampleFormat_Double
};

/*--------------------------------------------------------------------------------*/
/** Return number of bits per sample for a given sample format
 */
/*--------------------------------------------------------------------------------*/
uint8_t GetBitsPerSample(SampleFormat_t type)
{
  return SoundFormatBits[type];
}

/*--------------------------------------------------------------------------------*/
/** Return number of bytes per sample for a given sample format
 */
/*--------------------------------------------------------------------------------*/
uint8_t GetBytesPerSample(SampleFormat_t type)
{
  return (SoundFormatBits[type] + 7) >> 3;
}


/*--------------------------------------------------------------------------------*/
/** Perform sanity checks / adjustments for source/destination sample blocks
 *
 * @param src_channel starting channel to read from
 * @param src_channels total number of source channels in buffer
 * @param dst_channel starting channel to write to
 * @param dst_channels total number of destination channels in buffer
 * @param nchannels number of channels to transfer/convert
 * @param nframes number of frames to transfer/convert
 * @param allowsinglechannel true to allow many contiguous frames to be changed to a single frame of many channels
 *
 * @return true if resultant transfer is non-empty and valid
 *
 * @note allowsinglechannel should be FALSE when performing operations on a per-frame basis (like interpolation)
 */
/*--------------------------------------------------------------------------------*/
bool BlockTransferSanityChecks(uint_t& src_channel, uint_t& src_channels,
                               uint_t& dst_channel, uint_t& dst_channels,
                               uint_t& nchannels,
                               uint_t& nframes,
                               bool    allowsinglechannel)
{
  bool sane = false;
  
  // sanity checks
  if (src_channels && dst_channels &&
      nframes      && nchannels)
  {
    // restrict input data to sensible values
    src_channel = std::min(src_channel, src_channels - 1);
    dst_channel = std::min(dst_channel, dst_channels - 1);

    nchannels   = std::min(nchannels,   src_channels - src_channel);
    nchannels   = std::min(nchannels,   dst_channels - dst_channel);

    // final sanity check
    if (nchannels)
    {
      if (allowsinglechannel && (nchannels == src_channels) && (nchannels == dst_channels))
      {
        // optimization: if source and destination samples are all contiguous, reduce the process to a single frame of many channels
        nchannels *= nframes;
        nframes    = 1;
      }

      sane = true;
    }
  }

  return sane;
}

/*--------------------------------------------------------------------------------*/
/** Move/convert samples from one format to another and from one buffer to another
 *
 * @param vsrc pointer to source buffer (void * for easy casting) 
 * @param srctype format type of source
 * @param src_be true if source samples are big-endian
 * @param src_channel starting channel to read from
 * @param src_channels total number of source channels in buffer
 * @param vdst pointer to destination buffer (void * for easy casting)
 * @param dsttype format type of destination
 * @param dst_be true if destination samples are big-endian
 * @param dst_channel starting channel to write to
 * @param dst_channels total number of destination channels in buffer
 * @param nchannels number of channels to transfer/convert
 * @param nframes number of frames to transfer/convert
 * @param ditherer pointer to ditherer object or NULL
 *
 * @note this function provides copying, converting, interleaving and de-interleaving functionality
 * @note it allows converting in-place (i.e. dst == src)
 * @note HOWEVER, src and dst MUST NOT OVERLAP UNLESS dst == src
 *
 * The process can be viewed as the following:
 *
 * src:                          dst:
 * src_channel                   dst_channel
 * |-->|                         |----->|
 *     nchannels                        nchannels
 *     |---->|                          |---->|
 * +---+-----+-----+             +------+-----+----+   -
 * |   |sssss|     |             |      |ddddd|    |   |
 * |   |sssss|     |             |      |ddddd|    |   |
 * |   |sssss|     |             |      |ddddd|    |   |
 * |   |sssss|     |             |      |ddddd|    |   |
 * |   |sssss|     |             |      |ddddd|    |   |
 * |   |sss===============================>ddd|    |   | nframes
 * |   |sssss|     |             |      |ddddd|    |   |
 * |   |sssss|     |             |      |ddddd|    |   |
 * |   |sssss|     |             |      |ddddd|    |   |
 * |   |sssss|     |             |      |ddddd|    |   |
 * |   |sssss|     |             |      |ddddd|    |   v
 * +---+-----+-----+             +------+-----+----+   -
 * |-------------->|             |---------------->|
 *   src_channels                   dst_channels
 *
 * s = source format srctype / src_be
 * d = destination format dsttype / dst_be
 *
 * Notes:
 *   dst and src CAN be the same but otherwise MUST NOT overlap
 *
 *   0 <= src_channel <= (src_channels - nchannels)
 *   0 <= dst_channel <= (dst_channels - nchannels)
 *   0 <  nchannels   <= src_channels
 *   0 <  nchannels   <= dst_channels
 */
/*--------------------------------------------------------------------------------*/
void TransferSamples(const void *vsrc,       SampleFormat_t srctype, bool src_be,
                     uint_t     src_channel, uint_t src_channels,
                     void       *vdst,       SampleFormat_t dsttype, bool dst_be,
                     uint_t     dst_channel, uint_t dst_channels,
                     uint_t     nchannels,
                     uint_t     nframes,
                     Ditherer   *ditherer)
{
  // sanity checks
  if (BlockTransferSanityChecks(src_channel, src_channels,
                                dst_channel, dst_channels,
                                nchannels,
                                nframes) &&
      (srctype != SampleFormat_Unknown) && (srctype < SampleFormat_Count) &&
      (dsttype != SampleFormat_Unknown) && (dsttype < SampleFormat_Count))
  {
    const uint8_t *src    = (const uint8_t *)vsrc;
    uint8_t       *dst    = (uint8_t       *)vdst;
    sint_t        srclen  = GetBytesPerSample(srctype); // (signed so that the direction of operation can be backwards as well as forwards)
    sint_t        dstlen  = GetBytesPerSample(dsttype); // (signed so that the direction of operation can be backwards as well as forwards)
    sint_t        srcflen = src_channels * srclen;      // source frame length (signed so that the direction of operation can be backwards as well as forwards)
    sint_t        dstflen = dst_channels * dstlen;      // dest   frame length (signed so that the direction of operation can be backwards as well as forwards)

    // move to desired offsets (starting channel)
    src += src_channel * srclen;
    dst += dst_channel * dstlen;

    if (dstflen > srcflen)
    {
      // destination rectangle is BIGGER than source rectangle, switch to running backwards from the end of the buffers
      src += (nframes - 1) * srcflen;
      dst += (nframes - 1) * dstflen;
      srcflen = -srcflen;
      dstflen = -dstflen;
    }

    // look up correct function for conversion/copy
    CONVERTSAMPLES fn = SoundFormatConversions[(uint_t)src_be][(uint_t)dst_be][srctype][dsttype];
    if (!fn)
    {
      BBCERROR("Unknown copying routine for (%u/%u/%u/%u)", (uint_t)src_be, (uint_t)dst_be, srctype, dsttype);
      return;
    }

    // finally, run the functions
    (*fn)(src, dst, nchannels, nframes, srcflen, dstflen, ditherer);
  }
}

/*--------------------------------------------------------------------------------*/
/** Simple linear, contiguous transfer, (internal endianness -> internal endianness)
 */
/*--------------------------------------------------------------------------------*/
void TransferSamplesLinear(const void *vsrc, SampleFormat_t srctype,
                           void       *vdst, SampleFormat_t dsttype,
                           uint_t     nsamples,
                           Ditherer   *ditherer)
{
  // look up correct function for conversion/copy
  CONVERTSAMPLES fn = SoundFormatConversions[(uint_t)MACHINE_IS_BIG_ENDIAN][(uint_t)MACHINE_IS_BIG_ENDIAN][srctype][dsttype];
  if (!fn)
  {
    BBCERROR("Unknown copying routine for (%u/%u)", srctype, dsttype);
    return;
  }

  // finally, run the function
  (*fn)((const uint8_t *)vsrc, (uint8_t *)vdst, nsamples, 1, 0, 0, ditherer);
}

BBC_AUDIOTOOLBOX_END
