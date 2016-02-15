#ifndef __BI_QUAD__
#define __BI_QUAD__

#include <stdlib.h>
#include <string.h>

#include <vector>
#include <complex>
#ifdef __SSE3__
#  include <xmmintrin.h>
#endif

#include <bbcat-base/misc.h>

BBC_AUDIOTOOLBOX_START

/*--------------------------------------------------------------------------------*/
/** A class to manage coeffs for biquad filters, their calculation and interpolation
 *
 * The coeffs are used in the biquad calculation as show below:
 *
 *         num0 * x[n] + num1 * x[n - 1] + num2 * x[n - 2]
 * y[n] = ------------------------------------------------
 *              1 + den1 * y[n - 1] + den2 * y[n - 2]
 */
/*--------------------------------------------------------------------------------*/
class BiQuadCoeffs
{
public:
  // filter types
  typedef enum {
    FLAT,       // flat response
    LPF6,       // low pass filter (6dB/octave)
    HPF6,       // High pass filter (6dB/octave)
    LPF12,      // low pass filter (12dB/octave)
    HPF12,      // High pass filter (12dB/octave)
    BPF,        // band pass filter
    NOTCH,      // Notch Filter
    PEQ,        // Peaking band EQ filter
    LSH,        // Low shelf filter
    HSH,        // High shelf filter
  } Filter_t;

public:
  /*--------------------------------------------------------------------------------*/
  /** Construction from explicit set of coeffs
   */
  /*--------------------------------------------------------------------------------*/
  BiQuadCoeffs(double num0 = 1.0, double num1 = 0.0, double num2 = 0.0, double den1 = 0.0, double den2 = 0.0);
  /*--------------------------------------------------------------------------------*/
  /** Construct from filter description
   *
   * @param type filter type
   * @param freq filter centre frequency
   * @param fs sampling rate
   * @param gain gain in dB used for certain filters
   * @param bandwidth filter bandwidth
   * @param interp_time time in seconds for interpolation to complete
   *
   * @note default operation is for coeffs to be changed *immediately* to the new coeffs
   * @note set interp_time to not do this
   */
  /*--------------------------------------------------------------------------------*/
  BiQuadCoeffs(Filter_t type, double freq, double fs, double gain = 0.0, double bandwidth = 1.0, double interp_time = 0.0);
  /*--------------------------------------------------------------------------------*/
  /** Copy constructor
   */
  /*--------------------------------------------------------------------------------*/
  BiQuadCoeffs(const BiQuadCoeffs& obj);
  ~BiQuadCoeffs();

  /*--------------------------------------------------------------------------------*/
  /** Assignment operator
   */
  /*--------------------------------------------------------------------------------*/
  BiQuadCoeffs& operator = (const BiQuadCoeffs& obj);

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
  void SetCoeffs(double num0, double num1 = 0.0, double num2 = 0.0, double den1 = 0.0, double den2 = 0.0, double interp_samples = 0.0);

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
  typedef std::complex<double> RESPONSE;
  RESPONSE CalcResponseComplex(double f, double fs, bool usetargets = true) const;

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
  double CalcResponse(double f, double fs, bool usetargets = true) const;

  // coeff calculation routines taken from http://www.musicdsp.org/files/biquad.c
  //
  // Tom St Denis -- http://tomstdenis.home.dhs.org

  /*--------------------------------------------------------------------------------*/
  /** Calculate biquad coeffs for specified filter
   *
   * @param type filter type
   * @param freq filter centre frequency
   * @param fs sampling rate
   * @param gain gain in dB used for certain filters
   * @param bandwidth filter bandwidth
   * @param interp_time time in seconds for interpolation to complete
   *
   * @note default operation is for coeffs to be changed *immediately* to the new coeffs
   * @note set interp_time to not do this
   */
  /*--------------------------------------------------------------------------------*/
  void CalcCoeffs(Filter_t type, double freq, double fs, double gain = 0.0, double bandwidth = 1.0, double interp_time = 0.0);

  /*--------------------------------------------------------------------------------*/
  /** Interpolate coeffs from targets and diffs
   *
   * @param count number of iterations to interpolate
   */
  /*--------------------------------------------------------------------------------*/
  void Interpolate(double count = 1.0);

  typedef struct
  {
    double num0, num1, num2;    // numerator coeffs (a0, a1, a2 below)
    double den1, den2;          // denominator coeffs (b1, b2 below)
  } COEFFS;

  COEFFS current;               // current coeffs

  const COEFFS& GetTargets() const {return targets;}

protected:
  COEFFS targets;               // target coeffs
  COEFFS diffs;                 // initial differences between current and target coeffs
  double mul, dec;              // multiplier and decrement used for interpolation
};

/*--------------------------------------------------------------------------------*/
/** Simple class to implement a biquad filter
 *
 * Performs:
 *
 *         num0 * x[n] + num1 * x[n - 1] + num2 * x[n - 2]
 * y[n] = ------------------------------------------------
 *              1 + den1 * y[n - 1] + den2 * y[n - 2]
 *
 */
/*--------------------------------------------------------------------------------*/
class BiQuad
{
public:
  BiQuad(const BiQuadCoeffs& _coeffs);
  BiQuad(const BiQuad& obj);
  ~BiQuad();

  /*--------------------------------------------------------------------------------*/
  /** Copy audio state (NOT coeffs) from another biquad
   */
  /*--------------------------------------------------------------------------------*/
  void CopyAudioState(const BiQuad& obj);

  /*--------------------------------------------------------------------------------*/
  /** Reset internal state
   */
  /*--------------------------------------------------------------------------------*/
  void Reset() {memset(w, 0, sizeof(w));}

  /*--------------------------------------------------------------------------------*/
  /** Process a single sample through biquad
   *
   * @param x input sample
   * @param coeffs biquad coeffs
   *
   * @note implements biquad using direct-form II
   *
   * @return output sample
   */
  /*--------------------------------------------------------------------------------*/
  inline Sample_t Process(Sample_t x)
  {
    Sample_t y = (Sample_t)(x * coeffs.current.num0 + w[0]);
    w[0] = x * coeffs.current.num1 - y * coeffs.current.den1 + w[1];
    w[1] = x * coeffs.current.num2 - y * coeffs.current.den2;
    return y;
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
  void Process(const Sample_t *input, Sample_t *output, uint_t n);

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
  static void Process(BiQuad *filters, const Sample_t *src, Sample_t *dst, uint_t nchannels, uint_t nsrcchannels, uint_t ndstchannels, uint_t nframes, BiQuadCoeffs& coeffs);

protected:
  const BiQuadCoeffs& coeffs;
  double w[2];
};

/*--------------------------------------------------------------------------------*/
/** Array of filter banks - manage multiple biquads per channel
 */
/*--------------------------------------------------------------------------------*/
class BiQuadFilterBank
{
public:
  BiQuadFilterBank();
  BiQuadFilterBank(const BiQuadFilterBank& obj);
  ~BiQuadFilterBank();

  /*--------------------------------------------------------------------------------*/
  /** Set number of filters per channel
   */
  /*--------------------------------------------------------------------------------*/
  void SetFilters(uint_t n);

  /*--------------------------------------------------------------------------------*/
  /** Add a filter using existing coeffs
   */
  /*--------------------------------------------------------------------------------*/
  void AddFilter(const BiQuadCoeffs& coeffs);

  /*--------------------------------------------------------------------------------*/
  /** Get number of filters
   */
  /*--------------------------------------------------------------------------------*/
  uint_t GetFilters() const {return (uint_t)filters.size();}

  /*--------------------------------------------------------------------------------*/
  /** Reset internal state
   */
  /*--------------------------------------------------------------------------------*/
  void Reset();

  /*--------------------------------------------------------------------------------*/
  /** Set number of channels
   */
  /*--------------------------------------------------------------------------------*/
  void SetChannels(uint_t n);

  /*--------------------------------------------------------------------------------*/
  /** Get number of channels
   */
  /*--------------------------------------------------------------------------------*/
  uint_t GetChannels() const {return nchannels;}

  /*--------------------------------------------------------------------------------*/
  /** Return appropriate set of coeffs for a filter (or NULL)
   */
  /*--------------------------------------------------------------------------------*/
  BiQuadCoeffs *GetFilterCoeffs(uint_t i) {return (i < filters.size()) ? &filters[i]->coeffs : NULL;}

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
  void Process(const Sample_t *src, Sample_t *dst, uint_t nchannels, uint_t nsrcchannels, uint_t ndstchannels, uint_t nframes);

  /*--------------------------------------------------------------------------------*/
  /** Process a block of samples through each filter, one sample at a time, no coefficient
   * interpolation.
   *
   * @param src source buffer
   * @param dst destination buffer
   * @param nframes number of sample frames to process
   */
  /*--------------------------------------------------------------------------------*/
  void ProcessCascade(const Sample_t *src, Sample_t *dst, uint_t nframes);

  /*--------------------------------------------------------------------------------*/
  /** Calculate the current response of this filterbank at specific frequency
   *
   * @param f frequency
   * @param fs sampling rate
   * @param usetargets true to use target coefficients, false to use current coefficient values
   *
   * @return complex response
   */
  /*--------------------------------------------------------------------------------*/
  BiQuadCoeffs::RESPONSE CalcResponseComplex(double f, double fs, bool usetargets = true) const;

  /*--------------------------------------------------------------------------------*/
  /** Calculate current magnitude response in dB of this filterbank at specific frequency
   *
   * @param f frequency
   * @param fs sampling rate
   * @param usetargets true to use target coefficients, false to use current coefficient values
   *
   * @return response in dB
   */
  /*--------------------------------------------------------------------------------*/
  double CalcResponse(double f, double fs, bool usetargets = true) const;

protected:
  typedef struct {
    BiQuadCoeffs        coeffs;
    std::vector<BiQuad> channels;
  } FILTER;

protected:
  std::vector<FILTER *> filters;
  uint_t                nchannels;
};

/*--------------------------------------------------------------------------------*/
/** Simple structure to store the coefficients for a BiQuadCascade
 *
 */
/*--------------------------------------------------------------------------------*/
//typedef struct {
//  double                    gain;
//  std::vector<BiQuadCoeffs> biquads;
//} BiQuadCascadeCoeffs;

typedef struct {
  float gain;
  std::vector<float> b1;
  std::vector<float> b2;
  std::vector<float> a1;
  std::vector<float> a2;
} BiQuadCascadeCoeffs;

/*--------------------------------------------------------------------------------*/
/** Simple class to implement a biquad cascade filter.
 *  More restricted and therefore hopefully faster than the above classes.
 *  Maximum cascade length = maxnumfilters (i.e. (2*maxnumfilters)th order IIR)
 *  Data is organised to allow parallelisation.
 *
 *  @note: No coefficient interpolation is carried out with this filter
 *  @note: All data is single precision currently
 *  TODO: template? need to consider loop unrolling and vectorisation here
 *
 */
/*--------------------------------------------------------------------------------*/
//template <class T>
class BiQuadCascade
{
public:
  BiQuadCascade() : numfilters(1), vectorise(false), unroll(false), g(1.0) {};
  /*--------------------------------------------------------------------------------*/
  /** Construct default pass-through filter cascade
   *
   * @param numfilters Use numfilters biquads in cascade
   */
  /*--------------------------------------------------------------------------------*/
  BiQuadCascade(uint_t _numfilters, bool _vectorise = true, bool _unroll = true) :
    numfilters(_numfilters),
    vectorise(_vectorise),
    unroll(_unroll)
  {
    if (numfilters > maxnumfilters)
    {
      BBCERROR("Too many filters - max number of filters is %u", maxnumfilters);
      numfilters = 0;
    }
    if (vectorise && (_numfilters % 4))
    {
      BBCERROR("numfilters must be a factor of 4 when using vectorisation, disabling vectorisation");
      vectorise = false;
    }
    g = 1.0;
    for (uint_t i = 0; i < numfilters; i++)
    {
      b1[i] = 0.0; b2[i] = 0.0;
      a1[i] = 0.0; a2[i] = 0.0;
    }

  }
  /*--------------------------------------------------------------------------------*/
  /** Construct filter cascade from coefficients
   *
   * @param numfilters Use numfilters biquads in cascade
   * @param coefficients Interleaved biquad coefficient vector (length must be 4*numfilters + 1)
   *        (g,b1[0],b2[0],a1[0],a2[0],b1[1],b2[1],a1[1],a2[1]...)
   */
  /*--------------------------------------------------------------------------------*/
  BiQuadCascade(uint_t _numfilters, const std::vector<float>& coefficients,
                bool _vectorise = true, bool _unroll = true) :
                  numfilters(_numfilters),
                  vectorise(_vectorise),
                  unroll(_unroll)
  {
    CheckNumFilters(_numfilters);
    SetCoefficients(coefficients);
  }
  /*--------------------------------------------------------------------------------*/
  /** Construct filter cascade from coefficients
   *
   * @param numfilters Use numfilters biquads in cascade
   * @param g Filter output gain
   * @param b1 First  numerator   coefficients of filters (length numfilters)
   * @param b2 Second numerator   coefficients of filters (length numfilters)
   * @param a1 First  denominator coefficients of filters (length numfilters)
   * @param a2 Second denominator coefficients of filters (length numfilters)
   */
  /*--------------------------------------------------------------------------------*/
  BiQuadCascade(uint_t _numfilters, float _g,
                const std::vector<float> _b1, const std::vector<float> _b2,
                const std::vector<float> _a1, const std::vector<float> _a2,
                bool _vectorise = true,
                bool _unroll = true) :
                  numfilters(_numfilters),
                  vectorise(_vectorise),
                  unroll(_unroll),
                  g(_g)
  {
    CheckNumFilters(_numfilters);
    SetCoefficients(_g, _b1, _b2, _a1, _a2);
  }
  /*--------------------------------------------------------------------------------*/
  /** Destructor
   *
   */
  /*--------------------------------------------------------------------------------*/
  ~BiQuadCascade() {}

  /*--------------------------------------------------------------------------------*/
  /** Set the number of filters to use
   *
   * @note This resets the filter's internal registers
   */
  /*--------------------------------------------------------------------------------*/
  bool SetNumFilters(uint_t _numfilters, bool _vectorise = true, bool _unroll = true)
  {
    unroll    = _unroll;
    vectorise = _vectorise;
    if (!CheckNumFilters(_numfilters))
      return false;
    // clear coefficients for new filters
    if (_numfilters > numfilters)
    {
      for (uint_t i = numfilters; i < _numfilters; i++)
      {
        b1[i] = 0.0; b2[i] = 0.0;
        a1[i] = 0.0; a2[i] = 0.0;
      }
    }
    numfilters = _numfilters;
    Reset();
    return true;
  }

  /*--------------------------------------------------------------------------------*/
  /** Check that the number of required filters in the cascade is allowed and whether
   * it can be vectorised.
   */
  /*--------------------------------------------------------------------------------*/
  bool CheckNumFilters(uint_t _numfilters)
  {
    if (_numfilters > maxnumfilters)
    {
      BBCERROR("SetNumFilters: Too many filters - max number of filters is %u", maxnumfilters);
      return false;
    }
    if (vectorise && (_numfilters % 4))
    {
      BBCERROR("numfilters must be a factor of 4 when using vectorisation, disabling vectorisation");
      vectorise = false;
    }
    return true;
  }

  /*--------------------------------------------------------------------------------*/
  /** Reset filter cascade's registers
   */
  /*--------------------------------------------------------------------------------*/
  void Reset()
  {
    memset(&x,  0, sizeof(x));
    memset(&y,  0, sizeof(y));
    memset(&w0, 0, sizeof(w0));
    memset(&w1, 0, sizeof(w1));
    lastoutput = 0.0;
  }

  /*--------------------------------------------------------------------------------*/
  /** Set the filter cascade's coefficients from a single interleaved vector.
   *  This resets the filter's internal memory registers.
   */
  /*--------------------------------------------------------------------------------*/
  bool SetCoefficients(const std::vector<float> &coefficients)
  {
    if (coefficients.size() != (4*numfilters + 1))
    {
      BBCERROR("BiQuadCascade: coefficients vector must be 4*numfilters + 1 long");
      return false;
    }
    // unpack coefficients vector (g,b1[0],b2[0],a1[0],a2[0],b1[1],b2[1],a1[1],a2[1]...)
    std::vector<float>::const_iterator cf_it = coefficients.begin();
    g = *cf_it++;
    for (uint_t i = 0; i < numfilters; i++)
    {
      b1[i] = *cf_it++; b2[i] = *cf_it++;
      a1[i] = *cf_it++; a2[i] = *cf_it++;
    }
    if (cf_it != coefficients.end())
    {
      BBCERROR("BiQuadCascade: Something went wrong when loading the coefficients!");
      return false;
    }
    else
    {
      Reset();
      return true;
    }
  }

  /*--------------------------------------------------------------------------------*/
  /** Set the filter cascade's coefficients.
   *  This resets the filter's internal memory registers.
   *
   * @param g Filter output gain
   * @param b1 First  numerator   coefficients of filters (length numfilters)
   * @param b2 Second numerator   coefficients of filters (length numfilters)
   * @param a1 First  denominator coefficients of filters (length numfilters)
   * @param a2 Second denominator coefficients of filters (length numfilters)
   */
  /*--------------------------------------------------------------------------------*/
  bool SetCoefficients(float _g,
                       const std::vector<float> _b1, const std::vector<float> _b2,
                       const std::vector<float> _a1, const std::vector<float> _a2)
  {
    if ((_b1.size() != numfilters) && (_b1.size() != _b2.size())
        && (_b2.size() != _a1.size()) && (_a1.size() != _a2.size()))
    {
      BBCERROR("BiQuadCascade: coefficients vectors must all be numfilters long");
      return false;
    }

    g = _g;
    memcpy(&b1[0], &_b1[0], sizeof(float)*numfilters);
    memcpy(&b2[0], &_b2[0], sizeof(float)*numfilters);
    memcpy(&a1[0], &_a1[0], sizeof(float)*numfilters);
    memcpy(&a2[0], &_a2[0], sizeof(float)*numfilters);

    return true;
  }

#ifdef __SSE3__
  /*--------------------------------------------------------------------------------*/
  /** Process one sample through a set of four biquads in parallel using SSE3 intrinsics.
   *
   * @param i Start index in cascade buffers for this filter
   *
   * @note This method uses filter inputs and outputs.
   * @note Therefore causes a delay of numfilters samples through the cascade.
   * @note Outputs should be copied to inputs after all filters have ticked.
   */
  /*--------------------------------------------------------------------------------*/
  inline void VectorisedBiQuadTick_simd(uint_t i)
  {
    // load data
    __m128 xv      = _mm_load_ps((float *)&x [i]);
    __m128 w0v     = _mm_load_ps((float *)&w0[i]);
    __m128 w1v     = _mm_load_ps((float *)&w1[i]);
    __m128 b1v     = _mm_load_ps((float *)&b1[i]);
    __m128 a1v     = _mm_load_ps((float *)&a1[i]);
    __m128 b2v     = _mm_load_ps((float *)&b2[i]);
    __m128 a2v     = _mm_load_ps((float *)&a2[i]);
    // calculate output registers
    __m128 yv      = _mm_add_ps(xv, w0v);
    ;                _mm_store_ps((float *) &y[i], yv);
    // update first internal register
    __m128 tmp1    = _mm_mul_ps(xv, b1v);
    __m128 tmp2    = _mm_mul_ps(yv, a1v);
    __m128 tmp3    = _mm_sub_ps(tmp1, tmp2);
    __m128 w0v_new = _mm_add_ps(tmp3, w1v);
    ;                _mm_store_ps((float *) &w0[i], w0v_new);
    // update second internal register
    __m128 tmp4    = _mm_mul_ps(xv, b2v);
    __m128 tmp5    = _mm_mul_ps(yv, a2v);
    __m128 w1v_new = _mm_sub_ps(tmp4, tmp5);
    ;                _mm_store_ps((float *) &w1[i], w1v_new);
  }
#endif

  /*--------------------------------------------------------------------------------*/
  /** Process one sample through a set of four biquads in parallel.
   *
   * @param i Start index in cascade buffers for this filter
   *
   * @note This method uses filter inputs and outputs.
   * @note Therefore causes a delay of numfilters samples through the cascade.
   * @note Outputs should be copied to inputs after all filters have ticked.
   */
  /*--------------------------------------------------------------------------------*/
  inline void VectorisedBiQuadTick_cpp(uint_t i)
  {
    // calculate output registers
    y[i  ] = x[i  ] + w0[i  ];
    y[i+1] = x[i+1] + w0[i+1];
    y[i+2] = x[i+2] + w0[i+2];
    y[i+3] = x[i+3] + w0[i+3];
    // update first internal register
    w0[i  ] = x[i  ] * b1[i  ] - y[i  ] * a1[i  ] + w1[i  ];
    w0[i+1] = x[i+1] * b1[i+1] - y[i+1] * a1[i+1] + w1[i+1];
    w0[i+2] = x[i+2] * b1[i+2] - y[i+2] * a1[i+2] + w1[i+2];
    w0[i+3] = x[i+3] * b1[i+3] - y[i+3] * a1[i+3] + w1[i+3];
    // update second internal register
    w1[i  ] = x[i  ] * b2[i  ] - y[i  ] * a2[i  ];
    w1[i+1] = x[i+1] * b2[i+1] - y[i+1] * a2[i+1];
    w1[i+2] = x[i+2] * b2[i+2] - y[i+2] * a2[i+2];
    w1[i+3] = x[i+3] * b2[i+3] - y[i+3] * a2[i+3];
  }

  /*--------------------------------------------------------------------------------*/
  /** Process one sample through a single biquad
   *
   * @param i Start index in cascade buffers for this filter
   *
   * @note This method uses filter inputs and outputs, as with the vectorised methods.
   * @note Outputs should be copied to inputs after all filters have ticked.
   */
  /*--------------------------------------------------------------------------------*/
  inline void SingleBiQuadTick(uint_t i)
  {
    y[i]  = x[i] + w0[i];
    w0[i] = x[i] * b1[i] - y[i] * a1[i] + w1[i];
    w1[i] = x[i] * b2[i] - y[i] * a2[i];
  }

  /*--------------------------------------------------------------------------------*/
  /** Process one sample through the filter cascade
   *
   * @param in Input sample
   * @return Output sample
   */
  /*--------------------------------------------------------------------------------*/
  inline float Tick(float in)
  {
    uint_t i;

    if(vectorise)
    {
      x[0] = in;
      for (i = 0; (i+4) <= numfilters; i += 4)
      {
#ifdef __SSE3__
        VectorisedBiQuadTick_simd(i);
#else
        VectorisedBiQuadTick_cpp(i);
#endif
      }
      // copy filter outputs to filter inputs
      memcpy(&x[1], &y[0],sizeof(float)*(numfilters-1));
      return lastoutput = y[numfilters-1];
    }
    else
    {
      y[0]  = in + w0[0];
      w0[0] = in * b1[0] - y[0] * a1[0] + w1[0];
      w1[0] = in * b2[0] - y[0] * a2[0];
      for (i = 1; i < numfilters; i++)
      {
        // calculate output registers
        y [i] = y[i-1] + w0[i];
        w0[i] = y[i-1] * b1[i] - y[i] * a1[i] + w1[i];
        w1[i] = y[i-1] * b2[i] - y[i] * a2[i];
      }
      return lastoutput = y[numfilters-1];
    }
  }

  void ProcessCascade(const float *input, float *dest, uint_t blocksize)
  {
    if (unroll)
    {
      uint_t j;
      for (j = 0; (j+4) <= blocksize; j+=4)
      {
        dest[j]   = Tick(input[j]);
        dest[j+1] = Tick(input[j+1]);
        dest[j+2] = Tick(input[j+2]);
        dest[j+3] = Tick(input[j+3]);
      }
      for (; j < blocksize; j++)
        dest[j] = Tick(input[j]);
    }
    else
    {
      for (uint_t j = 0; j < blocksize; j++)
        dest[j] = Tick(input[j]);
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
  BiQuadCoeffs::RESPONSE CalcResponseComplex(float f, float fs) const
  {
    // loop over the filters, create BiQuadCoeffs for each one and multiply responses
    BiQuadCoeffs::RESPONSE r(1.0);
    BiQuadCoeffs bqcf;
    for (uint_t i = 0; i < numfilters; i++)
    {
      bqcf = BiQuadCoeffs(1.0, b1[i], b2[i], a1[i], a2[i]);
      r *= bqcf.CalcResponseComplex(f, fs);
    }
    return r;
  }

  /*--------------------------------------------------------------------------------*/
  /** Calculate magnitude response in dB at specific frequency
   *
   * @param f frequency
   * @param fs sampling rate
   *
   * @return response in dB
   */
  /*--------------------------------------------------------------------------------*/
  float CalcResponse(float f, float fs) const
  {
    return (float)(20.0 * log10(abs(CalcResponseComplex(f, fs))));
  }

protected:
  static const uint_t maxnumfilters = 12;
  uint_t   numfilters;
  bool     vectorise;
  bool     unroll;
  // TODO: can I used std::vectors here?
  // will the memory be on the stack and byte-aligned?
  // will it be as efficient?
  MEMALIGNED(16, float   b1[maxnumfilters] );  // first numerator coefficient of each filter
  MEMALIGNED(16, float   b2[maxnumfilters] );  // second numerator coefficient of each filter
  MEMALIGNED(16, float   a1[maxnumfilters] );  // first denominator coefficient of each filter
  MEMALIGNED(16, float   a2[maxnumfilters] );  // second denominator coefficient of each filter
  MEMALIGNED(16, float   x [maxnumfilters] );  // input registers for each filter
  MEMALIGNED(16, float   y [maxnumfilters] );  // output registers for each filter
  MEMALIGNED(16, float   w0[maxnumfilters] );  // first internal register of each filter
  MEMALIGNED(16, float   w1[maxnumfilters] );  // second internal register of each filter
  float   lastoutput;         // last output sample of filter cascade
  float   g;   // output gain of filter cascade
};

BBC_AUDIOTOOLBOX_END

#endif
