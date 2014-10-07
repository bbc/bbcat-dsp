pkgincludedir = $(includedir)/bbcat-dsp-@BBCAT_DSP_MAJORMINOR@/bbcat-dsp

lib_LTLIBRARIES = libap-dsp-@BBCAT_DSP_MAJORMINOR@.la

libap_dsp_@BBCAT_DSP_MAJORMINOR@_la_LDFLAGS =	\
	-no-undefined								\
	-version-info $(BBCAT_DSP_LIBVERSION)		\
    $(BBCAT_DSP_LDFLAGS)

libap_dsp_@BBCAT_DSP_MAJORMINOR@_la_LIBADD =	\
	$(BBCAT_GLOBAL_LIBS)						\
	$(PTHREAD_LIBS)								\
	$(BBCAT_BASE_LIBS)

libap_dsp_@BBCAT_DSP_MAJORMINOR@_la_CPPFLAGS =	\
	$(BBCAT_GLOBAL_CFLAGS)						\
	$(BBCAT_DSP_CFLAGS)							\
	$(PTHREAD_CFLAGS)							\
	$(BBCAT_BASE_CFLAGS)

libap_dsp_@BBCAT_DSP_MAJORMINOR@_la_SOURCES =	\
	FractionalSample.cpp						\
	SoundDelayBuffer.cpp						\
	SoundFormatConversions.cpp					\
	SoundFormatRawConversions.cpp				\
	SoundMixing.cpp

pkginclude_HEADERS =							\
	FractionalSample.h							\
	Interpolator.h								\
	SoundDelayBuffer.h							\
	SoundFormatConversions.h					\
	SoundFormatRawConversions.h					\
	SoundMixing.h

if ENABLE_CONVOLVER
libap_dsp_@BBCAT_DSP_MAJORMINOR@_la_CPPFLAGS += $(APF_CFLAGS) $(SNDFILE_CFLAGS) $(FFTW_CFLAGS)
libap_dsp_@BBCAT_DSP_MAJORMINOR@_la_LIBADD   += $(APF_LIBS) $(SNDFILE_LIBS) $(FFTW_LIBS)
libap_dsp_@BBCAT_DSP_MAJORMINOR@_la_SOURCES  += Convolver.cpp
pkginclude_HEADERS                           += Convolver.h
endif

if ENABLE_GPL
if ENABLE_SOFA
libap_dsp_@BBCAT_DSP_MAJORMINOR@_la_CPPFLAGS += $(NETCDF_CFLAGS) $(NETCDFCXX_CFLAGS)
libap_dsp_@BBCAT_DSP_MAJORMINOR@_la_LIBADD   += $(NETCDF_LIBS) $(NETCDFCXX_LIBS)
libap_dsp_@BBCAT_DSP_MAJORMINOR@_la_SOURCES  += SOFA.cpp
pkginclude_HEADERS                           += SOFA.h
endif
endif

noinst_HEADERS =
