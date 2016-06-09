AUTOMAKE_OPTIONS = subdir-objects
pkgincludedir = $(includedir)/bbcat-dsp-@BBCAT_DSP_MAJORMINOR@/bbcat-dsp

lib_LTLIBRARIES = libbbcat-dsp-@BBCAT_DSP_MAJORMINOR@.la

libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_LDFLAGS =	\
	-no-undefined									\
	-version-info $(BBCAT_DSP_LIBVERSION)			\
    $(BBCAT_DSP_LDFLAGS)

libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_LIBADD =	\
	$(BBCAT_BASE_LIBS)							\
	$(BBCAT_GLOBAL_DSP_LIBS)

libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_CPPFLAGS =	\
	$(BBCAT_BASE_CFLAGS)							\
	$(BBCAT_DSP_CFLAGS)								\
	$(BBCAT_GLOBAL_DSP_CFLAGS)

libbbcat_dsp_sources =							\
	BiQuad.cpp									\
	FractionalSample.cpp						\
	SoundDelayBuffer.cpp						\
	SoundFormatConversions.cpp					\
	SoundFormatRawConversions.cpp				\
	SoundMixing.cpp

pkginclude_HEADERS =							\
	AllPassFilter.h								\
	BiQuad.h									\
	FractionalSample.h							\
	Histogram.h									\
	Interpolator.h								\
	MultilayerBuffer.h							\
	RingBuffer.h								\
	RunningAverage.h							\
	SoundDelayBuffer.h							\
	SoundFormatConversions.h					\
	SoundFormatRawConversions.h					\
	SoundMixing.h								\
	register.h

noinst_HEADERS =

nodist_pkginclude_HEADERS =

distcleanfiles = 

#CLEANFILES = register.cpp

# register.cpp is included in repo but will be updated by the script below 
register.cpp: Makefile ../config.h ../gitbranch.txt $(libbbcat_dsp_sources)
	sh $(prefix)/share/bbcat-common/findregisterfunctions.sh "bbcat-dsp" "@VERSION@-`cat ../gitbranch.txt`" "$(libbbcat_dsp_sources)" "bbcat_base" >$@

libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_SOURCES = $(libbbcat_dsp_sources) register.cpp
