#ifndef __AUDIO_OBJECT_PARAMETERS__
#define __AUDIO_OBJECT_PARAMETERS__

#include <map>

#include <bbcat-base/3DPosition.h>
#include <bbcat-base/ParameterSet.h>

BBC_AUDIOTOOLBOX_START

/*--------------------------------------------------------------------------------*/
/** A class containing the parameters for rendering audio objects
 *
 * Each input channel to the renderer requires one of these objects which may change
 * over time as parameters for that channel change
 *
 * Each parameter (except the position and othervalues parameters) are either boolean,
 * integer or floating point.  The format is ENTIRELY dictated by the 4 interface
 * functions:
 *
 * double GetGain()            const {return Get_f(Parameter_gain);}
 * bool   GetGain(double& val) const {return Get_f(Parameter_gain, val);}
 * void   SetGain(double val)        {Set_f(Parameter_gain, val);}
 * void   UnsetGain()                {Set_f(Parameter_gain, 1.0, false);}
 *
 * All the Get_x() and Set_x() functions have the same type:
 * _f -> floating point (double)
 * _i -> integer (int)
 * _b -> boolean (bool)
 *
 * The first GetGain() function just returns the gain whether it was set or not
 * (in which case the default value is returned) - this function can be used when
 * just the value is required
 *
 * The second GetGain() function returns whether the gain has been set and puts
 * the gain into the referenced variable - this function can be used to get the
 * gain value and to decide whether to express the parameter in any file or
 * protocol (there is no need to output parameters that have not been set, this
 * function offers a way to determine that)
 *
 * The UnsetGain() function resets the value back to the default and clears the
 * 'set' flag so that the parameter does not get expressed in any file or
 * protocol
 *
 * Therefore, it is *extremely* important that all types and the default values
 * match up!
 */
/*--------------------------------------------------------------------------------*/
class AudioObjectParameters
{
public:
  AudioObjectParameters();
  AudioObjectParameters(const AudioObjectParameters& obj);
  virtual ~AudioObjectParameters();

  /*--------------------------------------------------------------------------------*/
  /** Assignment operator
   */
  /*--------------------------------------------------------------------------------*/
  virtual AudioObjectParameters& operator = (const AudioObjectParameters& obj);

  /*--------------------------------------------------------------------------------*/
  /** Comparison operator
   */
  /*--------------------------------------------------------------------------------*/
  virtual bool operator == (const AudioObjectParameters& obj) const;
  virtual bool operator != (const AudioObjectParameters& obj) const {return !operator == (obj);}

  /*--------------------------------------------------------------------------------*/
  /** Transform this object's position and return new copy
   */
  /*--------------------------------------------------------------------------------*/
  friend AudioObjectParameters operator * (const AudioObjectParameters& obj, const PositionTransform& transform);

  /*--------------------------------------------------------------------------------*/
  /** Transform this object's position
   */
  /*--------------------------------------------------------------------------------*/
  AudioObjectParameters& operator *= (const PositionTransform& transform);

  /*--------------------------------------------------------------------------------*/
  /** Get/Set physical position of this object
   *
   * @note position information is required for every channel
   */
  /*--------------------------------------------------------------------------------*/
  const Position& GetPosition() const {return position;}
  virtual void SetPosition(const Position& position) {this->position = position;}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set gain
   */
  /*--------------------------------------------------------------------------------*/
  double GetGain()            const {return Get_f(Parameter_gain);}
  bool   GetGain(double& val) const {return Get_f(Parameter_gain, val);}
  void   SetGain(double val)        {Set_f(Parameter_gain, val);}
  void   UnsetGain()                {Set_f(Parameter_gain, 1.0, false);}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set width
   */
  /*--------------------------------------------------------------------------------*/
  double GetWidth()            const {return Get_f(Parameter_width);}
  bool   GetWidth(double& val) const {return Get_f(Parameter_width, val);}
  void   SetWidth(double val)        {Set_f(Parameter_width, MAX(val, 0.0));}
  void   UnsetWidth()                {Set_f(Parameter_width, 0.0, false);}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set depth
   */
  /*--------------------------------------------------------------------------------*/
  double GetDepth()            const {return Get_f(Parameter_depth);}
  bool   GetDepth(double& val) const {return Get_f(Parameter_depth, val);}
  void   SetDepth(double val)        {Set_f(Parameter_depth, MAX(val, 0.0));}
  void   UnsetDepth()                {Set_f(Parameter_depth, 0.0, false);}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set height
   */
  /*--------------------------------------------------------------------------------*/
  double GetHeight()            const {return Get_f(Parameter_height);}
  bool   GetHeight(double& val) const {return Get_f(Parameter_height, val);}
  void   SetHeight(double val)        {Set_f(Parameter_height, MAX(val, 0.0));}
  void   UnsetHeight()                {Set_f(Parameter_height, 0.0, false);}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set diffuseness
   */
  /*--------------------------------------------------------------------------------*/
  double GetDiffuseness()            const {return Get_f(Parameter_diffuseness);}
  bool   GetDiffuseness(double& val) const {return Get_f(Parameter_diffuseness, val);}
  void   SetDiffuseness(double val)        {Set_f(Parameter_diffuseness, LIMIT(val, 0.0, 1.0));}
  void   UnsetDiffuseness()                {Set_f(Parameter_diffuseness, 0.0, false);}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set delay
   */
  /*--------------------------------------------------------------------------------*/
  double GetDelay()            const {return Get_f(Parameter_delay);}
  bool   GetDelay(double& val) const {return Get_f(Parameter_delay, val);}
  void   SetDelay(double val)        {Set_f(Parameter_delay, MAX(val, 0.0));}
  void   UnsetDelay()                {Set_f(Parameter_delay, 0.0, false);}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set importance
   */
  /*--------------------------------------------------------------------------------*/
  int    GetImportance()         const {return Get_i(Parameter_importance);}
  bool   GetImportance(int& val) const {return Get_i(Parameter_importance, val);}
  void   SetImportance(int val)        {Set_i(Parameter_importance, val);}
  void   UnsetImportance()             {Set_i(Parameter_importance, 0, false);}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set dialogue
   */
  /*--------------------------------------------------------------------------------*/
  int    GetDialogue()         const {return Get_i(Parameter_dialogue);}
  bool   GetDialogue(int& val) const {return Get_i(Parameter_dialogue, val);}
  void   SetDialogue(int val)        {Set_i(Parameter_dialogue, LIMIT(val, 0, 2));}
  void   UnsetDialogue()             {Set_i(Parameter_dialogue, 0, false);}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set channellock
   */
  /*--------------------------------------------------------------------------------*/
  bool   GetChannelLock()          const {return Get_b(Parameter_channellock);}
  bool   GetChannelLock(bool& val) const {return Get_b(Parameter_channellock, val);}
  void   SetChannelLock(bool val)        {Set_b(Parameter_channellock, val);}
  void   UnsetChannelLock()              {Set_b(Parameter_channellock, false, false);}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set interact
   */
  /*--------------------------------------------------------------------------------*/
  bool   GetInteract()          const {return Get_b(Parameter_interact);}
  bool   GetInteract(bool& val) const {return Get_b(Parameter_interact, val);}
  void   SetInteract(bool val)        {Set_b(Parameter_interact, val);}
  void   UnsetInteract()              {Set_b(Parameter_interact, false, false);}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set interpolate
   */
  /*--------------------------------------------------------------------------------*/
  bool   GetInterpolate()          const {return Get_b(Parameter_interpolate);}
  bool   GetInterpolate(bool& val) const {return Get_b(Parameter_interpolate, val);}
  void   SetInterpolate(bool val)        {Set_b(Parameter_interpolate, val);}
  void   UnsetInterpolate()              {Set_b(Parameter_interpolate, false, false);}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set onscreen
   */
  /*--------------------------------------------------------------------------------*/
  bool   GetOnScreen()          const {return Get_b(Parameter_onscreen);}
  bool   GetOnScreen(bool& val) const {return Get_b(Parameter_onscreen, val);}
  void   SetOnScreen(bool val)        {Set_b(Parameter_onscreen, val);}
  void   UnsetOnScreen()              {Set_b(Parameter_onscreen, false, false);}

  /*--------------------------------------------------------------------------------*/
  /** Get/Set supplementary data
   */
  /*--------------------------------------------------------------------------------*/
  const ParameterSet& GetOtherValues() const {return othervalues;}
  virtual void SetOtherValues(const ParameterSet& othervalues) {this->othervalues = othervalues;}

  /*--------------------------------------------------------------------------------*/
  /** Convert all parameters into text and store them in a ParameterSet object 
   *
   * @param parameters ParameterSet object to receive parameters
   */
  /*--------------------------------------------------------------------------------*/
  virtual void GetAll(ParameterSet& parameters) const;

  /*--------------------------------------------------------------------------------*/
  /** Set parameters from a ParameterSet object
   *
   * @param parameters ParameterSet object holding parameters
   */
  /*--------------------------------------------------------------------------------*/
  virtual void SetAll(const ParameterSet& parameters);

  /*--------------------------------------------------------------------------------*/
  /** Convert parameters to a string
   */
  /*--------------------------------------------------------------------------------*/
  std::string ToString(bool pretty = false) const;

protected:
  /*--------------------------------------------------------------------------------*/
  /** Initialise all parameters to defaults
   */
  /*--------------------------------------------------------------------------------*/
  virtual void InitialiseToDefaults();

  /*--------------------------------------------------------------------------------*/
  /** list of parameters, order are irrelevant
   */
  /*--------------------------------------------------------------------------------*/
  typedef enum {
    Parameter_gain = 0,

    Parameter_width,
    Parameter_depth,
    Parameter_height,

    Parameter_diffuseness,
    Parameter_delay,

    Parameter_importance,
    Parameter_dialogue,

    Parameter_channellock,
    Parameter_interact,
    Parameter_interpolate,
    Parameter_onscreen,

    Parameter_count,
  } Parameter_t;

  /*--------------------------------------------------------------------------------*/
  /** Get functions - return value from values or default value if not valid
   *
   * @param parameter Parameter_xxx parameter index
   *
   * @return value of parameter
   */
  /*--------------------------------------------------------------------------------*/
  bool   Get_b(Parameter_t parameter) const {return values[parameter].value.b;}
  int    Get_i(Parameter_t parameter) const {return values[parameter].value.i;}
  double Get_f(Parameter_t parameter) const {return values[parameter].value.f;}

  /*--------------------------------------------------------------------------------*/
  /** Get functions - set val from values
   *
   * @param parameter Parameter_xxx parameter index
   * @param val variable to receive value
   *
   * @return true if parameter has been set (i.e. default value not used)
   */
  /*--------------------------------------------------------------------------------*/
  bool   Get_b(Parameter_t parameter, bool&   val) const {val = values[parameter].value.b; return values[parameter].isset;}
  bool   Get_i(Parameter_t parameter, int&    val) const {val = values[parameter].value.i; return values[parameter].isset;}
  bool   Get_f(Parameter_t parameter, double& val) const {val = values[parameter].value.f; return values[parameter].isset;}

  /*--------------------------------------------------------------------------------*/
  /** Set functions, set value or default in values
   *
   * @param parameter Parameter_xxx parameter index
   * @param val value to set
   */
  /*--------------------------------------------------------------------------------*/
  void   Set_b(Parameter_t parameter, bool   val, bool set = true) {values[parameter].value.b = val; values[parameter].isset = set;}
  void   Set_i(Parameter_t parameter, int    val, bool set = true) {values[parameter].value.i = val; values[parameter].isset = set;}
  void   Set_f(Parameter_t parameter, double val, bool set = true) {values[parameter].value.f = val; values[parameter].isset = set;}

  typedef struct {
    bool isset;         ///< true if parameter has been explicitly set
    union {
      bool   b;         ///< boolean value
      int    i;         ///< integer value
      double f;         ///< floating-point value
    } value;
  } VALUE;

protected:
  Position     position;                      // audio channel position
  VALUE        values[Parameter_count];       // numerical values
  ParameterSet othervalues;                   // additional, arbitrary parameters
};

BBC_AUDIOTOOLBOX_END

#endif
