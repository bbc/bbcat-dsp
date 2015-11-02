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

libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la_SOURCES =	\
	FractionalSample.cpp							\
	SoundDelayBuffer.cpp							\
	SoundFormatConversions.cpp						\
	SoundFormatRawConversions.cpp					\
	SoundMixing.cpp									\
	version.cpp

pkginclude_HEADERS =							\
	FractionalSample.h							\
	Histogram.h									\
	Interpolator.h								\
	RingBuffer.h								\
	RunningAverage.h							\
	SoundDelayBuffer.h							\
	SoundFormatConversions.h					\
	SoundFormatRawConversions.h					\
	SoundMixing.h								\
	version.h

noinst_HEADERS =

nodist_pkginclude_HEADERS =

distcleanfiles = 

CLEANFILES = version.h version.cpp

version.h: Makefile ../config.h ../gitbranch.txt
	echo "{\"@PACKAGE@\", \"@VERSION@-`cat ../gitbranch.txt`\"}," >version.h

version.cpp: Makefile ../config.h ../gitbranch.txt
	echo "/* Auto-generated: DO NOT EDIT! */" >version.cpp
	echo "#include <bbcat-base/LoadedVersions.h>" >>version.cpp
	echo "BBC_AUDIOTOOLBOX_START" >>version.cpp
	echo "volatile uint_t __keep_@PACKAGE@_version = bbcat::LoadedVersions::Get().Register(\"@PACKAGE@\", \"@VERSION@-`cat ../gitbranch.txt`\");" >>version.cpp
	echo "BBC_AUDIOTOOLBOX_END" >>version.cpp
