# Arbitrarily set minimum version.
# This specific version is used to avoid complaints about the 
# cygwin platform not defining WIN32 anymore.
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.12.2)

# The name of our top-level project.
PROJECT( bbcat-dsp C CXX )

set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 1)
set(PROJECT_VERSION_RELEASE 1)
set(PROJECT_VERSION_BUILD 1)

#include all the parts that are consistent across all libraries
include(CMakeLists-bbcat.txt)

find_package(bbcat-base REQUIRED)
include_directories(${BBCAT_BASE_INCLUDE_DIRS})

set(GLOBAL_INCLUDE_DIRS
	${BBCAT_BASE_INCLUDE_DIRS}
	${INSTALL_INCLUDE_DIR})

set(GLOBAL_FLAGS
	${BBCAT_BASE_FLAGS})

ADD_DEFINITIONS(${GLOBAL_FLAGS})

################################################################################
# sources are contained in the src/ directory
ADD_SUBDIRECTORY( src )

################################################################################
# install files for 'share'
install(DIRECTORY "share/"
		DESTINATION "${INSTALL_SHARE_DIR}"
		PATTERN "share/*")
