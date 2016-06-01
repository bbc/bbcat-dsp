#sources
set(_sources
    BiQuad.cpp
	SoundDelayBuffer.cpp
	SoundFormatConversions.cpp
	SoundFormatRawConversions.cpp
	SoundMixing.cpp
)

# public headers
set(_headers
	AllPassFilter.h
    BiQuad.h
	FractionalSample.h
	Histogram.h
	Interpolator.h
	MultilayerBuffer.h
	RingBuffer.h
	RunningAverage.h
	SoundDelayBuffer.h
	SoundFormatConversions.h
	SoundFormatRawConversions.h
	SoundMixing.h
)


# auto-generate register.cpp
bbcat_generate_register_cpp(DEPENDENCIES "bbcat_base")

#include all the parts that are consistent across all projects
include(${BBCAT_COMMON_DIR}/CMakeLibrary.txt)

TARGET_LINK_LIBRARIES(bbcat-dsp bbcat-base)
