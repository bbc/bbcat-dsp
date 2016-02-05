#ifndef __RUNNING_AVERAGE__
#define __RUNNING_AVERAGE__

#include <vector>

#include <bbcat-base/misc.h>

BBC_AUDIOTOOLBOX_START

/*--------------------------------------------------------------------------------*/
/** Simple templated arbitrary type running average
 *
 * @param ITEMTYPE type of each item (MUST be a simple type)
 * @param SUMTYPE type for summation value (can be bigger than ITEMTYPE)
 */
/*--------------------------------------------------------------------------------*/
template<typename ITEMTYPE,typename SUMTYPE>
class RunningAverage
{
public:
  RunningAverage(uint_t n = 0) {SetLength(n);}
  RunningAverage(const RunningAverage& obj) : buffer(obj.buffer),
                                              pos(obj.pos),
                                              sum(obj.sum),
                                              wrapped(obj.wrapped) {}
  ~RunningAverage() {}

  /*--------------------------------------------------------------------------------*/
  /** Set running average length
   */
  /*--------------------------------------------------------------------------------*/
  void SetLength(uint_t n)
  {
    buffer.resize(n);
    Reset();
  }

  /*--------------------------------------------------------------------------------*/
  /** Return running average length
   */
  /*--------------------------------------------------------------------------------*/
  uint_t GetLength() const {return (uint_t)buffer.size();}

  /*--------------------------------------------------------------------------------*/
  /** Reset running average
   */
  /*--------------------------------------------------------------------------------*/
  void Reset()
  {
    // reset contents of buffer
    if (buffer.size()) buffer.assign(buffer.size(), (ITEMTYPE)0);
    pos     = 0;
    sum     = (SUMTYPE)0;
    wrapped = false;
  }

  /*--------------------------------------------------------------------------------*/
  /** Write item into buffer and update running average
   */
  /*--------------------------------------------------------------------------------*/
  void Write(ITEMTYPE val)
  {
    if (buffer.size())
    {
      sum        -= buffer[pos];
      buffer[pos] = val;
      sum        += val;
      if ((++pos) >= buffer.size()) {pos = 0; wrapped = true;}
    }
  }

  /*--------------------------------------------------------------------------------*/
  /** Write item into buffer and update running average
   */
  /*--------------------------------------------------------------------------------*/
  void Write(const ITEMTYPE *src, uint_t count, uint_t stride = 1)
  {
    if (buffer.size())
    {
      while (count)
      {
        uint_t i, n = std::min(count, buffer.size() - pos);

        for (i = 0; i < n; i++, src += stride)
        {
          sum -= buffer[pos];
          buffer[pos] = src[0];
          sum += buffer[pos++];
        }

        count -= n;
        if (pos >= buffer.size()) {pos = 0; wrapped = true;}
      }
    }
  }

  /*--------------------------------------------------------------------------------*/
  /** Update separate running average which re-uses this buffer of items
   *
   * @param sum separate sum to be updated
   * @param len length of running average (*must* be less than the size of this buffer)
   *
   * @return updated running average
   *
   * @note it is assumed that this is called *after* Write() is called
   */
  /*--------------------------------------------------------------------------------*/
  ITEMTYPE AltAverage(SUMTYPE& altsum, uint_t altlen) const
  {
    if (altlen < buffer.size())
    {
      uint_t lastpos = (pos + buffer.size() - 1) % buffer.size();       // use last position - last valid item
      altsum -= buffer[(lastpos + buffer.size() - altlen) % buffer.size()];
      altsum += buffer[lastpos];
      // if not enough data has been written yet, adjust altlen (divider)
      if (!wrapped && (pos < altlen)) altlen = pos; 
    }
    return (ITEMTYPE)(altsum / (SUMTYPE)altlen);
  }

  /*--------------------------------------------------------------------------------*/
  /** Return current running average
   */
  /*--------------------------------------------------------------------------------*/
  ITEMTYPE GetAverage() const {return buffer.size() ? (ITEMTYPE)(sum / (SUMTYPE)(wrapped ? buffer.size() : pos)) : 0;}

  /*--------------------------------------------------------------------------------*/
  /** Return current sum
   */
  /*--------------------------------------------------------------------------------*/
  SUMTYPE GetSum() const {return sum;}

protected:
  std::vector<ITEMTYPE> buffer;
  uint_t                pos;
  SUMTYPE               sum;
  bool                  wrapped;
};

BBC_AUDIOTOOLBOX_END

#endif
