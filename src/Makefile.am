pkgincludedir = $(includedir)/bbcat-dsp-@BBCAT_DSP_MAJORMINOR@/bbcat-dsp

lib_LTLIBRARIES = libbbcat-dsp-@BBCAT_DSP_MAJORMINOR@.la

libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_LDFLAGS =	\
	-no-undefined									\
	-version-info $(BBCAT_DSP_LIBVERSION)			\
    $(BBCAT_DSP_LDFLAGS)

libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_LIBADD =	\
	$(BBCAT_GLOBAL_LIBS)						\
	$(PTHREAD_LIBS)								\
	$(BBCAT_BASE_LIBS)

libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_CPPFLAGS =	\
	$(BBCAT_GLOBAL_CFLAGS)							\
	$(BBCAT_DSP_CFLAGS)								\
	$(PTHREAD_CFLAGS)								\
	$(BBCAT_BASE_CFLAGS)

libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_SOURCES =	\
	FractionalSample.cpp							\
	SoundDelayBuffer.cpp							\
	SoundFormatConversions.cpp						\
	SoundFormatRawConversions.cpp					\
	SoundMixing.cpp

pkginclude_HEADERS =							\
	FractionalSample.h							\
	Interpolator.h								\
	SoundDelayBuffer.h							\
	SoundFormatConversions.h					\
	SoundFormatRawConversions.h					\
	SoundMixing.h

if ENABLE_CONVOLVER
libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_CPPFLAGS += $(APF_CFLAGS) $(SNDFILE_CFLAGS) $(FFTW_CFLAGS)
libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_LIBADD   += $(APF_LIBS) $(SNDFILE_LIBS) $(FFTW_LIBS)
libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_SOURCES  += Convolver.cpp BlockConvolver.cpp
pkginclude_HEADERS                              += Convolver.h BlockConvolver.h
endif

if ENABLE_GPL
if ENABLE_SOFA
libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_CPPFLAGS += $(NETCDF_CFLAGS) $(NETCDFCXX_CFLAGS)
libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_LIBADD   += $(NETCDF_LIBS) $(NETCDFCXX_LIBS)
libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_SOURCES  += SOFA.cpp
pkginclude_HEADERS                              += SOFA.h
endif
endif

noinst_HEADERS =
