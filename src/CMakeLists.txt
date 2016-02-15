#sources
set(_sources
	BiQuad.cpp
	FractionalSample.cpp
	SoundDelayBuffer.cpp
	SoundFormatConversions.cpp
	SoundFormatRawConversions.cpp
	SoundMixing.cpp
	${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)

# public headers
set(_headers
	AllPassFilter.h
	BiQuad.h
	FractionalSample.h
	Histogram.h
	Interpolator.h
	RingBuffer.h
	RunningAverage.h
	SoundDelayBuffer.h
	SoundFormatConversions.h
	SoundFormatRawConversions.h
	SoundMixing.h
	${CMAKE_CURRENT_BINARY_DIR}/version.h
)

#include all the parts that are consistent across all libraries
include(CMakeLists-src.txt)

TARGET_LINK_LIBRARIES(bbcat-dsp bbcat-base)
