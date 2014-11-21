LDADD = $(BBCAT_GLOBAL_LIBS) \
        $(PTHREAD_LIBS) \
        $(BBCAT_BASE_LIBS) \
        $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        $(SNDFILE_LIBS) \
        $(FFTW_LIBS)

AM_CPPFLAGS = $(BBCAT_GLOBAL_CFLAGS) \
              $(BBCAT_DSP_CFLAGS) \
              $(PTHREAD_CFLAGS) \
              $(BBCAT_BASE_CFLAGS)

check_PROGRAMS =
TESTS =

BlockConvolverTest_SOURCES = BlockConvolverTest.cpp
BlockConvolverTest_LDADD = $(LDADD) $(top_srcdir)/src/libbbcat_dsp_@BBCAT_DSP_MAJORMINOR@_la-BlockConvolver.lo
check_PROGRAMS += BlockConvolverTest
TESTS += BlockConvolverTest
