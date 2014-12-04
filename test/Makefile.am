LDADD = $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        ../src/libbbcat-dsp-@BBCAT_DSP_MAJORMINOR@.la

LDADD = $(BBCAT_GLOBAL_LIBS) \
        $(PTHREAD_LIBS) \
        $(BBCAT_BASE_LIBS) \
        $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        $(SNDFILE_LIBS) \
        $(FFTW_LIBS) \
        ../src/libbbcat-dsp-@BBCAT_DSP_MAJORMINOR@.la

AM_CPPFLAGS = $(BBCAT_GLOBAL_CFLAGS) \
              $(BBCAT_DSP_CFLAGS) \
              $(PTHREAD_CFLAGS) \
              $(BBCAT_BASE_CFLAGS)

check_PROGRAMS =
TESTS =

BlockConvolverTest_SOURCES = BlockConvolverTest.cpp

check_PROGRAMS += BlockConvolverTest
TESTS += BlockConvolverTest
