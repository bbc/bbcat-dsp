
#include "SoundMixing.h"

BBC_AUDIOTOOLBOX_START
                          
/*--------------------------------------------------------------------------------*/
/** Mix source samples to destination samples (like TransferSamples() but adding instead of over-writing) with level changing
 *
 * @param src pointer to source buffer
 * @param src_channel starting channel to read from
 * @param src_channels total number of source channels in buffer
 * @param dst pointer to destination buffer
 * @param dst_channel starting channel to write to
 * @param dst_channels total number of destination channels in buffer
 * @param nchannels number of channels to transfer/convert
 * @param nframes number of frames to transfer/convert
 * @param interp level interpolation control
 * @param inc level interpolation rate
 *
 * This function is like the above except that the scaling factor can change linearly using the interp object and inc parameter
 */
/*--------------------------------------------------------------------------------*/
void MixSamples(const Sample_t *src,
                uint_t src_channel, uint_t src_channels,
                Sample_t *dst,
                uint_t dst_channel, uint_t dst_channels,
                uint_t nchannels,
                uint_t nframes,
                Interpolator& interp, Sample_t inc)
{
  // sanity checks
  if (BlockTransferSanityChecks(src_channel, src_channels,
                                dst_channel, dst_channels,
                                nchannels,
                                nframes,
                                false) &&       // cannot allow optimisation to single channel because of interpolation!
      interp.NonZero())
  {
    // move to desired offsets (starting channel)
    src += src_channel;
    dst += dst_channel;

    Sample_t mul = interp;      // current level
    uint_t i, j;
    for (i = 0; i < nframes; i++, src += src_channels, dst += dst_channels)
    {
      for (j = 0; j < nchannels; j++) dst[j] += mul * src[j];
      interp += inc;            // interpolate level
      mul     = interp;         // get new level
    }
  }
}

BBC_AUDIOTOOLBOX_END
