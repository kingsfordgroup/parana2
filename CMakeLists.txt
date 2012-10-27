CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

# Use GCC
SET( CMAKE_CXX_COMPILER g++ )

# When Clang is ready
#SET( CMAKE_CXX_COMPILER /opt/local/bin/clang++)

MESSAGE( STATUS "Building in ${BUILDENV} Environment" )

Project( parana )

SET(CMAKE_VERBOSE_MAKEFILE ON)

SET( CMAKE_BUILD_TYPE Release)
SET( CMAKE_CXX_FLAGS "-O3 -pg -g -std=c++11 -Wl,-no_pie -Wno-deprecated ${CMAKE_CXX_FLAGS}")

# When Clang is ready
#SET( CMAKE_CXX_FLAGS "-O3 -g -std=c++11 -stdlib=libstdc++ ${CMAKE_CXX_FLAGS}")


FIND_PACKAGE( Boost 1.42 COMPONENTS regex program_options timer system )

SET( PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR} )
SET( CMAKE_INSTALL_PREFIX ${PROJECT_ROOT} )
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_ROOT}/build/bin )

SUBDIRS( src )
