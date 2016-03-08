
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BBCDEBUG_LEVEL 1
#include <bbcat-base/EnhancedFile.h>

#include "BiQuad.h"

BBC_AUDIOTOOLBOX_START

BiQuadCoeffs::BiQuadCoeffs(double num0, double num1, double num2, double den1, double den2) : mul(0.0),
                                                                                              dec(1.0)
{
  memset(&current, 0, sizeof(current));
  memset(&diffs,   0, sizeof(diffs));

  current.num0 = num0;
  current.num1 = num1;
  current.num2 = num2;
  current.den1 = den1;
  current.den2 = den2;
  targets = current;
}

BiQuadCoeffs::BiQuadCoeffs(Filter_t type, double freq, double fs, double gain, double bandwidth, double interp_time) : mul(0.0),
                                                                                                                       dec(1.0)
{
  memset(&current, 0, sizeof(current));
  memset(&diffs,   0, sizeof(diffs));

  current.num0 = 1.0;

  CalcCoeffs(type, freq, fs, gain, bandwidth, interp_time);
}

BiQuadCoeffs::BiQuadCoeffs(const BiQuadCoeffs& obj) : current(obj.current),
                                                      targets(obj.targets),
                                                      diffs(obj.diffs),
                                                      mul(obj.mul),
                                                      dec(obj.dec)
{
}

BiQuadCoeffs::~BiQuadCoeffs()
{
}

/*--------------------------------------------------------------------------------*/
/** Assignment operator
 */
/*--------------------------------------------------------------------------------*/
BiQuadCoeffs& BiQuadCoeffs::operator = (const BiQuadCoeffs& obj)
{
  current = obj.current;
  targets = obj.targets;
  diffs   = obj.diffs;
  mul     = obj.mul;
  dec     = obj.dec;

  return *this;
}

/*--------------------------------------------------------------------------------*/
/** Set from explicit set of coeffs
 *
 * @param interp_samples time in SAMPLES for interpolation to complete
 *
 * @note note difference between interp_samples and interp_time used in other functions!
 * @note default operation is for coeffs to be changed *immediately* to the new coeffs
 * @note set interp_samples to not do this
 */
/*--------------------------------------------------------------------------------*/
void BiQuadCoeffs::SetCoeffs(double num0, double num1, double num2, double den1, double den2, double interp_samples)
{
  targets.num0 = num0;
  targets.num1 = num1;
  targets.num2 = num2;
  targets.den1 = den1;
  targets.den2 = den2;

  // calculate differences between current and target for interpolation
  diffs.num0 = targets.num0 - current.num0;
  diffs.num1 = targets.num1 - current.num1;
  diffs.num2 = targets.num2 - current.num2;
  diffs.den1 = targets.den1 - current.den1;
  diffs.den2 = targets.den2 - current.den2;

  if (interp_samples > 0.0)
  {
    // set initial multiplier and decrement rate
    mul = 1.0;
    dec = 1.0 / interp_samples;
  }
  else
  {
    // zero time interpolation - just set multiplier at zero and current to targets
    mul = dec = 0.0;
    current = targets;
  }
}

/*--------------------------------------------------------------------------------*/
/** Calculate response at specific frequency
 *
 * @param f frequency
 * @param fs sampling rate
 * @param usetargets true to use targets, false to use current values
 *
 * @return complex response
 */
/*--------------------------------------------------------------------------------*/
BiQuadCoeffs::RESPONSE BiQuadCoeffs::CalcResponseComplex(double f, double fs, bool usetargets) const
{
  static const RESPONSE mul = 2.0 * M_PI * RESPONSE(0, 1.0);    // = 2.pi.j
  static const RESPONSE one(1.0);                               // = 1
  const COEFFS& coeffs = usetargets ? targets : current;
  RESPONSE z1 = exp(RESPONSE(f / fs * mul));                    // z^-1
  RESPONSE z2 = exp(RESPONSE(f / fs * 2.0 * mul));              // z^-2

  // Calculate
  //
  //      num0 + num1.z^-1 + num2.z^-2
  // H = -----------------------------
  //       1 + den1.z^-1 + den2.z^-2

  return ((coeffs.num0 + coeffs.num1 * z1 + coeffs.num2 * z2) /
          (one         + coeffs.den1 * z1 + coeffs.den2 * z2));
}

/*--------------------------------------------------------------------------------*/
/** Calculate magnitude response in dB at specific frequency
 *
 * @param f frequency
 * @param fs sampling rate
 * @param usetargets true to use targets, false to use current values
 *
 * @return response in dB
 */
/*--------------------------------------------------------------------------------*/
double BiQuadCoeffs::CalcResponse(double f, double fs, bool usetargets) const
{
  return 20.0 * log10(abs(CalcResponseComplex(f, fs, usetargets)));
}

/*--------------------------------------------------------------------------------*/
/** Calculate biquad coeffs for specified filter
 *
 * @param coeffs coeffs to be written
 * @param type filter type
 * @param freq filter centre frequency
 * @param fs sampling rate
 * @param gain gain in dB used for certain filters
 * @param bandwidth filter bandwidth
 * @param interp_time time in seconds for interpolation to complete
 *
 *
 * @note default operation is for coeffs to be changed *immediately* to the new coeffs
 * @note set interp_time to not do this
 *
Based on the work:

Cookbook formulae for audio EQ biquad filter coefficients
---------------------------------------------------------
by Robert Bristow-Johnson, pbjrbj@viconet.com  a.k.a. robert@audioheads.com

 * Available on the web at

http://www.smartelectronix.com/musicdsp/text/filters005.txt

 * Enjoy.
 *
 * This work is hereby placed in the public domain for all purposes, whether
 * commercial, free [as in speech] or educational, etc.  Use the code and please
 * give me credit if you wish.
 *
 * Tom St Denis -- http://tomstdenis.home.dhs.org
 */
/*--------------------------------------------------------------------------------*/
void BiQuadCoeffs::CalcCoeffs(Filter_t type, double freq, double fs, double gain, double bandwidth, double interp_time)
{
  double  A, omega, sn, cs, alpha, beta;
  // numerator coeffs
  double& b0 = targets.num0;
  double& b1 = targets.num1;
  double& b2 = targets.num2;
  // denominator coeffs
  double  a0 = 1.0;             // this will be used to normalize the coeffs at the end
  double& a1 = targets.den1;
  double& a2 = targets.den2;

  // setup variables
  A     = pow(10.0, gain / 40.0);
  omega = 2.0 * M_PI * freq / fs;
  sn    = sin(omega);
  cs    = cos(omega);
  alpha = sn * sinh(M_LN2 / 2.0 * bandwidth * omega / sn);
  beta  = sqrt(A + A);

  switch (type)
  {
    default:
    case FLAT:
      // flat
      b0 = 1.0;
      b1 = 0.0;
      b2 = 0.0;
      a0 = 1.0;
      a1 = 0.0;
      a2 = 0.0;
      break;

    case LPF6:
      // lowpass filter (6dB/octave) - NOT from Audio EQ Cookbook
      //          sn
      // = ----------------
      //   (1 + sn) - z(-1)
      b0 = sn;
      b1 = 0;
      b2 = 0;
      a0 = 1 + sn;
      a1 = -1;
      a2 = 0;
      break;

    case LPF12:
      // lowpass filter (12dB/octave) - NOT from Audio EQ Cookbook
      //          sn                 sn                           sn^2
      // = ---------------- . ---------------- = -------------------------------------
      //   (1 + sn) - z(-1)   (1 + sn) - z(-1)   (1 + sn)^2 - 2.(1 + sn).z(-1) + z(-2)
      b0 = sn * sn;
      b1 = 0;
      b2 = 0;
      a0 = (1 + sn) * (1 + sn);
      a1 = -2 * (1 + sn);
      a2 = 1;
      break;

    case HPF6:
      // highpass filter (6dB/octave) - NOT from Audio EQ Cookbook
      //       1 - z(-1)
      // = ------------------
      //   1 - (1 - sn).z(-1)
      b0 = 1;
      b1 = -1;
      b2 = 0;
      a0 = 1;
      a1 = -(1 - sn);
      a2 = 0;
      break;

    case HPF12:
      // highpass filter (12dB/octave) - NOT from Audio EQ Cookbook
      //        1 - z(-1)            1 - z(-1)                 1 - 2.z(-1) + z(-2)
      // = ------------------ . ------------------ = ---------------------------------------
      //   1 - (1 - sn).z(-1)   1 - (1 - sn).z(-1)   1 - 2.(1 - sn).z(-1) + (1 - sn)^2.z(-2)
      b0 = 1;
      b1 = -2;
      b2 = 1;
      a0 = 1;
      a1 = -2 * (1 - sn);
      a2 = (1 - sn) * (1 - sn);
      break;

    case BPF:
      // bandpass filter
      b0 = alpha;
      b1 = 0;
      b2 = -alpha;
      a0 = 1 + alpha;
      a1 = -2 * cs;
      a2 = 1 - alpha;
      break;

    case NOTCH:
      // notch reject filter
      b0 = 1;
      b1 = -2 * cs;
      b2 = 1;
      a0 = 1 + alpha;
      a1 = -2 * cs;
      a2 = 1 - alpha;
      break;

    case PEQ:
      // peaking filter
      b0 = 1 + (alpha * A);
      b1 = -2 * cs;
      b2 = 1 - (alpha * A);
      a0 = 1 + (alpha / A);
      a1 = -2 * cs;
      a2 = 1 - (alpha / A);
      break;

    case LSH:
      // lowshelf filter
      b0 = A * ((A + 1) - (A - 1) * cs + beta * sn);
      b1 = 2 * A * ((A - 1) - (A + 1) * cs);
      b2 = A * ((A + 1) - (A - 1) * cs - beta * sn);
      a0 = (A + 1) + (A - 1) * cs + beta * sn;
      a1 = -2 * ((A - 1) + (A + 1) * cs);
      a2 = (A + 1) + (A - 1) * cs - beta * sn;
      break;

    case HSH:
      // highshelf filter
      b0 = A * ((A + 1) + (A - 1) * cs + beta * sn);
      b1 = -2 * A * ((A - 1) + (A + 1) * cs);
      b2 = A * ((A + 1) + (A - 1) * cs - beta * sn);
      a0 = (A + 1) - (A - 1) * cs + beta * sn;
      a1 = 2 * ((A - 1) - (A + 1) * cs);
      a2 = (A + 1) - (A - 1) * cs - beta * sn;
      break;
  }

  {
    // normalize coeffs
    double normalise = 1.0 / a0;
    b0 *= normalise;
    b1 *= normalise;
    b2 *= normalise;
    a1 *= normalise;
    a2 *= normalise;
  }

  // calculate differences between current and target for interpolation
  diffs.num0 = targets.num0 - current.num0;
  diffs.num1 = targets.num1 - current.num1;
  diffs.num2 = targets.num2 - current.num2;
  diffs.den1 = targets.den1 - current.den1;
  diffs.den2 = targets.den2 - current.den2;

  if (interp_time > 0.0)
  {
    // set initial multiplier and decrement rate
    mul = 1.0;
    dec = 1.0 / (interp_time * fs);
  }
  else
  {
    // zero time interpolation - just set multiplier at zero and current to targets
    mul = dec = 0.0;
    current = targets;
  }

  BBCDEBUG2(("Coeffs: type %u freq %0.1lfHz gain %0.1lfdB bw %0.6lfHz -> {%0.6lf, %0.6lf, %0.6lf, %0.6lf, %0.6lf}",
          (uint_t)type, freq, gain, bandwidth,
          b0, b1, b2, a1, a2));

#if BBCDEBUG_LEVEL >= 3
  {
    // write response out to file
    EnhancedFile fp("coeffs.dat", "w");
    double fs     = 48000.0;
    double f1     = 10.0, f2 = 22000.0;
    double frange = log(f2 / f1);
    uint_t i, steps = 1000;

    for (i = 0; i < steps; i++)
    {
      double p = (double)i / (double)(steps - 1);
      double f = f1 * exp(p * frange);
      
      double gain = CalcResponse(f, fs);

      fp.fprintf("%u %0.1lf %0.4le\n", i, f, gain);
    }
  }
#endif
}

/*--------------------------------------------------------------------------------*/
/** Interpolate coeffs from targets and diffs
 *
 * @param count number of iterations to interpolate
 */
/*--------------------------------------------------------------------------------*/
void BiQuadCoeffs::Interpolate(double count)
{
  // only interpolate if not at targets yet
  if (mul > 0.0)
  {
    // decrement multiplier towards 0
    mul -= dec * count;
    mul  = std::max(mul, 0.0);

    // coeffs reach their targets when mul == 0
    current.num0 = targets.num0 - mul * diffs.num0;
    current.num1 = targets.num1 - mul * diffs.num1;
    current.num2 = targets.num2 - mul * diffs.num2;
    current.den1 = targets.den1 - mul * diffs.den1;
    current.den2 = targets.den2 - mul * diffs.den2;
  }
}

/*----------------------------------------------------------------------------------------------------*/

BiQuad::BiQuad(const BiQuadCoeffs& _coeffs) : coeffs(&_coeffs)
{
  memset(w, 0, sizeof(w));
}

BiQuad::BiQuad(const BiQuad& obj) : coeffs(obj.coeffs)
{
  memcpy(w, obj.w, sizeof(w));
}

BiQuad::~BiQuad()
{
}


/*--------------------------------------------------------------------------------*/
/** Copy audio state (NOT coeffs) from another biquad
 */
/*--------------------------------------------------------------------------------*/
void BiQuad::CopyAudioState(const BiQuad& obj)
{
  memcpy(w, obj.w, sizeof(w));
}

/*--------------------------------------------------------------------------------*/
/** Process an array of samples
 *
 * @param input input sample array
 * @param output output sample array
 * @param n number of samples to process
 *
 * @note this function does not provide ANY interpolation of coeffs during processing
 */
/*--------------------------------------------------------------------------------*/
void BiQuad::Process(const Sample_t *input, Sample_t *output, uint_t n)
{
  uint_t i;
#define unroll 1
#if unroll
  for (i = 0; i < n; i+=4)
  {
    output[i]   = Process(input[i]);
    output[i+1] = Process(input[i+1]);
    output[i+2] = Process(input[i+2]);
    output[i+3] = Process(input[i+3]);
  }
  for (i -= 3; i < n; i++)
  {
    output[i] = Process(input[i]);
  }
#else
  for (i = 0; i < n; i++)
  {
    output[i] = Process(input[i]);
  }
#endif
}

/*--------------------------------------------------------------------------------*/
/** Process a block of samples across multiple channels with coeff interpolation
 *
 * @param filters array of biquad filters (one for each channel)
 * @param src source array
 * @param dst destination array (can be the same as src)
 * @param nchannels number of channels to process (= number of filters)
 * @param nsrcchannels number of source channels
 * @param ndstchannels number of destination channels
 * @param nframes number of frames to process
 *
 * @note array of filters must be *at least* nchannels long!
 * @note the *assumption* is that each filter is using the same coeffs as specified in the arguments but this is not necessary
 * @note (each filter will use its own coeffs)
 */
/*--------------------------------------------------------------------------------*/
void BiQuad::Process(BiQuad *filters, const Sample_t *src, Sample_t *dst, uint_t nchannels, uint_t nsrcchannels, uint_t ndstchannels, uint_t nframes, BiQuadCoeffs& coeffs)
{
  uint_t i, j;

  // sanity checking
  nchannels = std::min(nchannels, nsrcchannels);
  nchannels = std::min(nchannels, ndstchannels);

  // process a block of samples using a list of filters
  for (i = 0; i < nframes; i++, src += nsrcchannels, dst += ndstchannels)
  {
    // process each channel through its own filter
    for (j = 0; j < nchannels; j++)
    {
      BBCDEBUG4(("Biquad[%u][%u]: {%0.8lf, %0.8lf, %0.8lf, %0.8lf, %0.8lf}", i, j, coeffs.current.num0, coeffs.current.num1, coeffs.current.num2, coeffs.current.den1, coeffs.current.den2));
      dst[j] = filters[j].Process(src[j]);
    }

    // interpolate coeffs
    coeffs.Interpolate();
  }
}

/*----------------------------------------------------------------------------------------------------*/

BiQuadFilterBank::BiQuadFilterBank() : nchannels(0)
{
}

BiQuadFilterBank::BiQuadFilterBank(const BiQuadFilterBank& obj) : nchannels(0)
{
  // duplicate passed object include audio state

  SetChannels(obj.nchannels);
  SetFilters((uint_t)obj.filters.size());
  
  uint_t i, j;
  for (i = 0; i < filters.size(); i++)
  {
    const FILTER& srcfilter = *obj.filters[i];
    FILTER&       dstfilter = *filters[i];

    // copy filter coeffs and interpolation state
    dstfilter.coeffs = srcfilter.coeffs;

    for (j = 0; j < srcfilter.channels.size(); j++)
    {
      // copy audio state
      dstfilter.channels[j].CopyAudioState(srcfilter.channels[j]);
    }
  }
}

BiQuadFilterBank::~BiQuadFilterBank()
{
  // delete all filters
  SetFilters(0);
}

/*--------------------------------------------------------------------------------*/
/** Set number of filters per channel
 */
/*--------------------------------------------------------------------------------*/
void BiQuadFilterBank::SetFilters(uint_t n)
{
  // remove filters from end if necessary
  while (filters.size() > n)
  {
    FILTER *filter = filters.back();
    delete filter;
    filters.pop_back();
  }

  // add filters if necessary
  while (filters.size() < n)
  {
    FILTER *filter;

    if ((filter = new FILTER) != NULL)
    {
      filters.push_back(filter);
    }
  }

  // ensure the added filters have the right number of channels
  if (nchannels) SetChannels(nchannels);
}

/*--------------------------------------------------------------------------------*/
/** Add a filter using existing coeffs
 */
/*--------------------------------------------------------------------------------*/
void BiQuadFilterBank::AddFilter(const BiQuadCoeffs& coeffs)
{
  FILTER *filter;

  if ((filter = new FILTER) != NULL)
  {
    filter->coeffs = coeffs;
    filters.push_back(filter);
  }

  // ensure the added filter has the right number of channels
  if (nchannels) SetChannels(nchannels);
}

/*--------------------------------------------------------------------------------*/
/** Set number of channels
 */
/*--------------------------------------------------------------------------------*/
void BiQuadFilterBank::SetChannels(uint_t n)
{
  uint_t i;

  nchannels = n;
  for (i = 0; i < filters.size(); i++)
  {
    FILTER& filter = *filters[i];

    // remove filters if necessary
    while (filter.channels.size() > n)
    {
      filter.channels.pop_back();
    }

    // add filters if necessary
    while (filter.channels.size() < n)
    {
      // create a filter for each channel using the correct set of coeffs
      filter.channels.push_back(BiQuad(filter.coeffs));
    }
  }  
}

/*--------------------------------------------------------------------------------*/
/** Reset internal state
 */
/*--------------------------------------------------------------------------------*/
void BiQuadFilterBank::Reset()
{
  uint_t i;

  for (i = 0; i < filters.size(); i++)
  {
    FILTER& filter = *filters[i];
    uint_t  j;

    // reset state of all filters
    for (j = 0; j < filter.channels.size(); j++)
    {
      filter.channels[j].Reset();
    }
  }
}

/*--------------------------------------------------------------------------------*/
/** Process a block of samples through each filter
 *
 * @param src source buffer
 * @param dst destination buffer
 * @param nchannels number of channels to process
 * @param nsrcchannels number of source channels
 * @param ndstchannels number of destination channels
 * @param nframes number of sample frames to process
 */
/*--------------------------------------------------------------------------------*/
void BiQuadFilterBank::Process(const Sample_t *src, Sample_t *dst, uint_t nchannels, uint_t nsrcchannels, uint_t ndstchannels, uint_t nframes)
{
  uint_t i;

  // sanity checking
  nchannels = std::min(nchannels, nsrcchannels);
  nchannels = std::min(nchannels, ndstchannels);

  for (i = 0; i < filters.size(); i++)
  {
    FILTER& filter = *filters[i];
    uint_t  n      = std::min(nchannels, (uint_t)filter.channels.size());    // number of channels to process

    BBCDEBUG4(("Processing biquad %u/%u %u channels (%u/%u) * %u frames", i, (uint_t)filters.size(), n, nsrcchannels, ndstchannels, nframes));

    BiQuad::Process(&filter.channels[0], src, dst, n, nsrcchannels, ndstchannels, nframes, filter.coeffs);

    // because each filter must process the previous filter's output,
    // set the source and source channel count to that of the destination
    // so that subsequent processing takes place in-place
    src = dst;
    nsrcchannels = ndstchannels;
  }
}

/*--------------------------------------------------------------------------------*/
/** Process a block of samples through each filter, one sample at a time, no coefficient
 * interpolation.
 *
 * @param src source buffer
 * @param dst destination buffer
 * @param nframes number of sample frames to process
 *
 * @note Warning: for speed, this function doesn't check whether each filter has the same
 * number of channels as the
 */
/*--------------------------------------------------------------------------------*/
void BiQuadFilterBank::ProcessCascade(const Sample_t *src, Sample_t *dst, uint_t nframes)
{
#define filter_first 0
#if filter_first
  uint_t i;

  for (i = 0; i < filters.size(); i++)
  {
    BBCDEBUG4(("Processing biquad %u/%u %u frames", i, (uint_t)filters.size(), nframes));
    filters[i]->channels[0].Process(src, dst, nframes);

    // because each filter must process the previous filter's output,
    // set the source and source channel count to that of the destination
    // so that subsequent processing takes place in-place
    src = dst;
  }
#else // sample first
  uint_t i, j;
  for (j = 0; j < nframes; j++, src++, dst++)
  {
    *dst = filters[0]->channels[0].Process(*src);
    for (i = 1; i < filters.size(); i++)
    {
      *dst = filters[i]->channels[0].Process(*dst);
    }
  }
#endif
}

/*--------------------------------------------------------------------------------*/
/** Calculate the response of this filterbank at specific frequency
 *
 * @param f frequency
 * @param fs sampling rate
 * @param usetargets true to use target coefficients, false to use current coefficient values
 *
 * @return complex response
 */
/*--------------------------------------------------------------------------------*/
BiQuadCoeffs::RESPONSE BiQuadFilterBank::CalcResponseComplex(double f, double fs, bool usetargets) const
{
  // loop over the filters, create BiQuadCoeffs for each one and multiply responses
  BiQuadCoeffs::RESPONSE r(1.0);
  for (uint_t i = 0; i < filters.size(); i++)
  {
    r *= filters[i]->coeffs.CalcResponseComplex(f, fs, usetargets);
  }
  return r;
}

/*--------------------------------------------------------------------------------*/
/** Calculate magnitude response in dB of this filterbank at specific frequency
 *
 * @param f frequency
 * @param fs sampling rate
 * @param usetargets true to use target coefficients, false to use current coefficient values
 *
 * @return response in dB
 */
/*--------------------------------------------------------------------------------*/
double BiQuadFilterBank::CalcResponse(double f, double fs, bool usetargets) const
{
  return 20.0 * log10(abs(CalcResponseComplex(f, fs, usetargets)));
}

BBC_AUDIOTOOLBOX_END
