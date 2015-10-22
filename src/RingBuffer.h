#ifndef __RING_BUFFER__
#define __RING_BUFFER__

#include <vector>

#include <bbcat-base/misc.h>

BBC_AUDIOTOOLBOX_START

template<typename ITEMTYPE>
class RingBuffer
{
public:
  RingBuffer(uint_t n = 0) {SetLength(n);}
  RingBuffer(const RingBuffer& obj) : buffer(obj.buffer),
                                      pos(obj.pos) {}
  ~RingBuffer() {}

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
  uint_t GetLength() const {return buffer.size();}

  /*--------------------------------------------------------------------------------*/
  /** Return current position
   */
  /*--------------------------------------------------------------------------------*/
  uint_t GetPosition() const {return pos;}

  /*--------------------------------------------------------------------------------*/
  /** Return delayed position
   */
  /*--------------------------------------------------------------------------------*/
  uint_t GetDelayedPosition(uint_t delay) const {return (pos + buffer.size() - delay) % buffer.size();}

  /*--------------------------------------------------------------------------------*/
  /** Reset position
   */
  /*--------------------------------------------------------------------------------*/
  void ResetPosition() {pos = 0;}

  /*--------------------------------------------------------------------------------*/
  /** Reset contents of buffer
   */
  /*--------------------------------------------------------------------------------*/
  void Reset()
  {
    // reset contents of buffer
    if (buffer.size()) buffer.assign(buffer.size(), (ITEMTYPE)0);
    pos = 0;
  }

  /*--------------------------------------------------------------------------------*/
  /** Write item into buffer and update running average
   */
  /*--------------------------------------------------------------------------------*/
  void Write(ITEMTYPE val)
  {
    if (buffer.size())
    {
      buffer[pos] = val;
      if ((++pos) >= buffer.size()) pos = 0;
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
        uint_t i, n = MIN(count, buffer.size() - pos);

        if (stride == 1)
        {
          memcpy(&buffer[pos], src, n * sizeof(*src));
          src += n;
          pos += n;
        }
        else
        {
          for (i = 0; i < n; i++, src += stride)
          {
            buffer[pos++] = src[0];
          }
        }

        count -= n;
        if (pos >= buffer.size()) pos = 0;
      }
    }
  }

  /*--------------------------------------------------------------------------------*/
  /** Return buffer for specific position and the maximum number of items that can be read
   */
  /*--------------------------------------------------------------------------------*/
  const ITEMTYPE *GetBuffer(uint_t rpos = 0, uint_t *maxitems = NULL) const
  {
    if (buffer.size())
    {
      if (maxitems) *maxitems = buffer.size() - rpos;
      return &buffer[rpos];
    }
    return NULL;
  } 

  /*--------------------------------------------------------------------------------*/
  /** Return delayed buffer ptr and the maximum number of items that can be read
   */
  /*--------------------------------------------------------------------------------*/
  const ITEMTYPE *GetDelayedBuffer(uint_t delay = 0, uint_t *maxitems = NULL) const
  {
    return GetBuffer(GetDelayedPosition(delay), maxitems);
  } 

protected:
  std::vector<ITEMTYPE> buffer;
  uint_t                pos;
};

BBC_AUDIOTOOLBOX_END

#endif
