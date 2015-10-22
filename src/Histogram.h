#ifndef __HISTOGRAM__
#define __HISTOGRAM__

#include <vector>

#include <bbcat-base/EnhancedFile.h>

BBC_AUDIOTOOLBOX_START

/*--------------------------------------------------------------------------------*/
/** Histogram object, accumulates data (of type ITEMTYPE) at indices (of type INDEXTYPE)
 */
/*--------------------------------------------------------------------------------*/
template<typename INDEXTYPE, typename ITEMTYPE>
class Histogram
{
public:
  /*--------------------------------------------------------------------------------*/
  /** Construct from min value, max value and step
   */
  /*--------------------------------------------------------------------------------*/
  Histogram(INDEXTYPE _min  = INDEXTYPE(0),
            INDEXTYPE _max  = INDEXTYPE(1),
            INDEXTYPE _step = INDEXTYPE(1)) {SetRange(_min, _max, _step);}
  Histogram(const Histogram& obj) : min(obj.min),
                                    range(obj.range),
                                    histogram(obj.histogram),
                                    total(obj.total) {}
  ~Histogram() {}

  /*--------------------------------------------------------------------------------*/
  /** Set range of histogram from min value, max value and step
   */
  /*--------------------------------------------------------------------------------*/
  void SetRange(INDEXTYPE _min, INDEXTYPE _max, INDEXTYPE _step = INDEXTYPE(1))
  {
    min   = _min;
    range = _max - _min;

    uint_t n = (uint_t)ceil(range / _step);
    histogram.resize(n);
    Reset();
  }
  uint_t GetSize() const {return histogram.size();}

  /*--------------------------------------------------------------------------------*/
  /** Each index holds count and sum
   */
  /*--------------------------------------------------------------------------------*/
  typedef struct {
    uint_t   count;
    ITEMTYPE data;
  } ITEM;

  /*--------------------------------------------------------------------------------*/
  /** Reset data
   */
  /*--------------------------------------------------------------------------------*/
  void Reset()
  {
    ITEM dummy = {0, ITEMTYPE(0)};
    histogram.assign(histogram.size(), dummy);
    total = dummy;
  }

  /*--------------------------------------------------------------------------------*/
  /** Return min and range
   */
  /*--------------------------------------------------------------------------------*/
  INDEXTYPE GetMin()   const {return min;}
  INDEXTYPE GetRange() const {return range;}
  INDEXTYPE GetMax()   const {return min + range;}

  /*--------------------------------------------------------------------------------*/
  /** Add data to histogram
   */
  /*--------------------------------------------------------------------------------*/
  void Add(INDEXTYPE index, ITEMTYPE val = ITEMTYPE(0))
  {
    ITEM& item = histogram[CalcIndex(index)];
    item.count++;
    item.data += val;
    total.count++;
    total.data += val;
  }

  /*--------------------------------------------------------------------------------*/
  /** Return data
   */
  /*--------------------------------------------------------------------------------*/
  const std::vector<ITEM>& GetData() const {return histogram;}

  /*--------------------------------------------------------------------------------*/
  /** Return total count and sum
   */
  /*--------------------------------------------------------------------------------*/
  const ITEM& GetTotal() const {return total;}

  /*--------------------------------------------------------------------------------*/
  /** Return index into data for given index
   */
  /*--------------------------------------------------------------------------------*/
  uint_t CalcIndex(INDEXTYPE index) const
  {
    index = (INDEXTYPE(histogram.size()) * (index - min)) / range;
    return (uint_t)LIMIT((sint_t)index, 0, (sint_t)(histogram.size() - 1));
  }

  /*--------------------------------------------------------------------------------*/
  /** Return index into data for given index
   */
  /*--------------------------------------------------------------------------------*/
  INDEXTYPE CalcReversedIndex(uint_t index) const
  {
    return min + ((INDEXTYPE(2) * range * (INDEXTYPE)index + INDEXTYPE(1)) / (INDEXTYPE)(2 * histogram.size()));
  }

  /*--------------------------------------------------------------------------------*/
  /** Calculate mean index (traditional histogram mean)
   */
  /*--------------------------------------------------------------------------------*/
  INDEXTYPE CalculateMeanIndex(INDEXTYPE start, INDEXTYPE end) const
  {
    INDEXTYPE sum = INDEXTYPE(0);
    INDEXTYPE div = INDEXTYPE(0);
    uint_t    i, i1 = CalcIndex(start), i2 = CalcIndex(end);

    for (i = i1; i <= i2; i++)
    {
      INDEXTYPE val = INDEXTYPE(histogram[i].count);
      // scale index by a factor of two to use the centre of the bin
      sum += INDEXTYPE(2 * i + 1) * val;
      div += val;
    }

    // divide the result by two to compensate for the scaling above
    return (div != INDEXTYPE(0)) ? sum / (INDEXTYPE(2) * div) : INDEXTYPE(0);
  }
  INDEXTYPE CalculateMeanIndex() const {return CalculateMeanIndex(GetMin(), GetMax());}

  /*--------------------------------------------------------------------------------*/
  /** Calculate mean of data, ignoring indices
   */
  /*--------------------------------------------------------------------------------*/
  ITEMTYPE CalculateMeanData(INDEXTYPE start, INDEXTYPE end) const
  {
    ITEMTYPE sum   = ITEMTYPE(0);
    uint_t   count = 0;
    uint_t   i, i1 = CalcIndex(start), i2 = CalcIndex(end);

    for (i = i1; i <= i2; i++)
    {
      const ITEM& item = histogram[i];
      sum   += item.data;
      count += item.count;
    }

    return count ? sum / ITEMTYPE(count) : ITEMTYPE(0);
  }
  ITEMTYPE CalculateMeanData() const {return CalculateMeanData(GetMin(), GetMax());}

  /*--------------------------------------------------------------------------------*/
  /** Calculate percentile values of indices
   *
   * For each index, returns the percentage of density stored in that index and below
   */
  /*--------------------------------------------------------------------------------*/
  void CalculateIndexPercentiles(std::vector<float>& list) const
  {
    list.resize(histogram.size());
    if (total.count)
    {
      const double mul = 100.0 / (double)total.count;
      double sum = 0.0;
      uint_t i;

      for (i = 0; i < histogram.size(); i++)
      {
        sum    += (double)histogram[i].count;
        list[i] = (float)(mul * sum);
      }
    }
    else list.assign(list.size(), 0.0);
  }

  /*--------------------------------------------------------------------------------*/
  /** Calculate percentile values of data
   *
   * For each index, returns the percentage of data stored in that index and below
   */
  /*--------------------------------------------------------------------------------*/
  void CalculateDataPercentiles(std::vector<float>& list) const
  {
    list.resize(histogram.size());
    if (total.data > (ITEMTYPE)0.0)
    {
      const double mul = 100.0 / (double)total.data;
      double sum = 0.0;
      uint_t i;

      for (i = 0; i < histogram.size(); i++)
      {
        sum    += (double)histogram[i].data;
        list[i] = (float)(mul * sum);
      }
    }
    else list.assign(list.size(), 0.0);
  }

  /*--------------------------------------------------------------------------------*/
  /** Write histogram data to file (for debugging purposes)
   */
  /*--------------------------------------------------------------------------------*/
  bool WriteToFile(const char *filename) const
  {
    EnhancedFile file;
    bool         success = false;
    if (file.fopen(filename, "w"))
    {
      std::vector<float> percentiles;
      ITEM   sum = {0, 0};
      uint_t i;
      CalculateIndexPercentiles(percentiles);
      for (i = 0; i < histogram.size(); i++)
      {
        const ITEM& item = histogram[i];
        sum.count += item.count;
        sum.data  += item.data;
        file.fprintf("%u %0.3lf %0.3f %u %0.14le %u %0.14le %u %0.14le\n",
                     i, (double)CalcReversedIndex(i),
                     percentiles[i],
                     item.count, (double)item.data,
                     sum.count, (double)sum.data,
                     total.count, (double)total.data);
      }
      file.fclose();
      success = true;
    }
    return success;
  }

protected:
  INDEXTYPE         min, range;
  std::vector<ITEM> histogram;
  ITEM              total;
};

BBC_AUDIOTOOLBOX_END

#endif
