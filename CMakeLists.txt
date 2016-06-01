# Arbitrarily set minimum version.
# This specific version is used to avoid complaints about the
# cygwin platform not defining WIN32 anymore.
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.12.2)

# The name of our top-level project.
PROJECT( bbcat-dsp C CXX )

set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 1)
set(PROJECT_VERSION_RELEASE 2)
set(PROJECT_VERSION_BUILD 2)

set(PROJECT_DESCRIPTION "Digital Signal Processing library")

# use bbcat-common as installed by bbcat-base
set(BBCAT_COMMON_DIR ${CMAKE_INSTALL_PREFIX}/share/bbcat-common)

#include all the parts that are consistent across all projects
include(${BBCAT_COMMON_DIR}/CMakeProject.txt)

find_package(bbcat-base REQUIRED)
bbcat_add_library(BBCAT_BASE)

# find ENABLE_GPL option in flags
is_feature_enabled(ENABLE_GPL "GPL support")

# set flags for compiling
add_definitions(${GLOBAL_FLAGS})

# set include directories for compiling
include_directories(${GLOBAL_INCLUDE_DIRS})

# set library directories for linking
link_directories(${GLOBAL_LIB_DIRS})

ADD_DEFINITIONS(${GLOBAL_FLAGS})

################################################################################
# sources are contained in the src/ directory
ADD_SUBDIRECTORY( src )

################################################################################
# install files for 'share'
install(DIRECTORY "share/"
		DESTINATION "${INSTALL_SHARE_DIR}"
		PATTERN "share/*")

# generate package config file
bbcat_generate_pkg_config()
