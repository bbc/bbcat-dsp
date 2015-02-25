
#include <stdlib.h>
#include <string.h>

#define DEBUG_LEVEL 1
#include "AudioObjectParameters.h"

BBC_AUDIOTOOLBOX_START

AudioObjectParameters::AudioObjectParameters()
{
  InitialiseToDefaults();
}

AudioObjectParameters::AudioObjectParameters(const AudioObjectParameters& obj)
{
  InitialiseToDefaults();
  operator = (obj);
}

AudioObjectParameters::~AudioObjectParameters()
{
}

/*--------------------------------------------------------------------------------*/
/** Initialise all parameters to defaults
 */
/*--------------------------------------------------------------------------------*/
void AudioObjectParameters::InitialiseToDefaults()
{
  memset(values, 0, sizeof(values));
  // set values to their default and 'unset'
  // only need to do this with non-zero parameters
  UnsetGain();

  position    = Position();
  othervalues = ParameterSet();
}

/*--------------------------------------------------------------------------------*/
/** Assignment operator
 */
/*--------------------------------------------------------------------------------*/
AudioObjectParameters& AudioObjectParameters::operator = (const AudioObjectParameters& obj)
{
  position    = obj.position;
  memcpy(values, obj.values, sizeof(values));
  othervalues = obj.othervalues;
  return *this;
}

/*--------------------------------------------------------------------------------*/
/** Comparison operator
 */
/*--------------------------------------------------------------------------------*/
bool AudioObjectParameters::operator == (const AudioObjectParameters& obj) const
{
  return ((position    == obj.position) &&
          (memcmp(values, obj.values, sizeof(values)) == 0) &&
          (othervalues == obj.othervalues));
}

/*--------------------------------------------------------------------------------*/
/** Transform this object's position and return new copy
 */
/*--------------------------------------------------------------------------------*/
AudioObjectParameters operator * (const AudioObjectParameters& obj, const PositionTransform& transform)
{
  AudioObjectParameters res = obj;
  res.SetPosition(obj.GetPosition() * transform);
  return res;
}

/*--------------------------------------------------------------------------------*/
/** Transform this object's position
 */
/*--------------------------------------------------------------------------------*/
AudioObjectParameters& AudioObjectParameters::operator *= (const PositionTransform& transform)
{
  SetPosition(GetPosition() * transform);
  return *this;
}

/*--------------------------------------------------------------------------------*/
/** Convert all parameters into text and store them in a ParameterSet object 
 *
 * @param parameters ParameterSet object to receive parameters
 */
/*--------------------------------------------------------------------------------*/
void AudioObjectParameters::GetAll(ParameterSet& parameters) const
{
  double fval;
  int    ival;
  bool   bval;

  position.SetParameters(parameters, "position");
  if (GetGain(fval))        parameters.Set("gain", fval);
  if (GetWidth(fval))       parameters.Set("width", fval);
  if (GetDepth(fval))       parameters.Set("depth", fval);
  if (GetHeight(fval))      parameters.Set("height", fval);
  if (GetDiffuseness(fval)) parameters.Set("diffuseness", fval);
  if (GetDelay(fval))       parameters.Set("delay", fval);
  if (GetImportance(ival))  parameters.Set("importance", ival);
  if (GetDialogue(ival))    parameters.Set("dialogue", ival);
  if (GetChannelLock(bval)) parameters.Set("channellock", bval);
  if (GetInteract(bval))    parameters.Set("interact", bval);
  if (GetInterpolate(bval)) parameters.Set("interpolate", bval);
  if (GetOnScreen(bval))    parameters.Set("onscreen", bval);

  parameters.Set("othervalues", othervalues);
}

/*--------------------------------------------------------------------------------*/
/** Set parameters from a ParameterSet object
 *
 * @param parameters ParameterSet object holding parameters
 */
/*--------------------------------------------------------------------------------*/
void AudioObjectParameters::SetAll(const ParameterSet& parameters)
{
  double fval;
  int    ival;
  bool   bval;

  position.GetFromParameters(parameters, "position");
  if (parameters.Get("gain", fval))        SetGain(fval);
  if (parameters.Get("width", fval))       SetWidth(fval);
  if (parameters.Get("depth", fval))       SetDepth(fval);
  if (parameters.Get("height", fval))      SetHeight(fval);
  if (parameters.Get("diffuseness", fval)) SetDiffuseness(fval);
  if (parameters.Get("delay", fval))       SetDelay(fval);
  if (parameters.Get("importance", ival))  SetImportance(ival);
  if (parameters.Get("dialogue", ival))    SetDialogue(ival);
  if (parameters.Get("channellock", bval)) SetChannelLock(bval);
  if (parameters.Get("interact", bval))    SetInteract(bval);
  if (parameters.Get("interpolate", bval)) SetInterpolate(bval);
  if (parameters.Get("onscreen", bval))    SetOnScreen(bval);
  
  parameters.Get("othervalues", othervalues);
}

/*--------------------------------------------------------------------------------*/
/** Convert parameters to a string
 */
/*--------------------------------------------------------------------------------*/
std::string AudioObjectParameters::ToString(bool pretty) const
{
  ParameterSet params;

  GetAll(params);

  return params.ToString(pretty);
}

BBC_AUDIOTOOLBOX_END
