/* Auto-generated: DO NOT EDIT! */
#include <bbcat-base/LoadedVersions.h>
BBC_AUDIOTOOLBOX_START
// list of libraries this library is dependant on
extern bool bbcat_register_bbcat_base();

// list of this library's component registration functions

// registration function
bool bbcat_register_bbcat_dsp()
{
  static bool registered = false;
  // prevent registration functions being called more than once
  if (!registered)
  {
    registered = true;
    // register other libraries
	bbcat_register_bbcat_base();

    // register this library's version number
    LoadedVersions::Get().Register("bbcat-dsp", "0.1.2.2-master");
    // register this library's components

  }
  return registered;
}
// automatically call registration functions
volatile const bool bbcat_dsp_registered = bbcat_register_bbcat_dsp();
BBC_AUDIOTOOLBOX_END
