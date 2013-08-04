CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Use GCC
# SET( CMAKE_CXX_COMPILER g++-mp-4.7 )

# When Clang is ready
#SET( CMAKE_CXX_COMPILER /opt/local/bin/clang++)

MESSAGE( STATUS "Building in ${BUILDENV} Environment" )

Project( Parana2 )

IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   MESSAGE ("-- Detected OSX")   
   SET (OperatingSystem "Mac OS X")
   SET (SHARED_LIB_EXT "dylib")
ELSEIF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   MESSAGE ("-- Detected Linux")
   SET (OperatingSystem "Linux")  
   SET (SHARED_LIB_EXT "so")
ELSE ()
   MESSAGE ("WHAT's THE OS? ${CMAKE_SYSTEM_NAME}")
ENDIF()

SET(CMAKE_VERBOSE_MAKEFILE ON)

SET( CMAKE_BUILD_TYPE Release)
SET( CMAKE_CXX_FLAGS "-O3 -DBOOST_HAS_HASH -DHAVE_ANSI_TERM -DNDEBUG -arch x86_64 -std=c++11 -Wl,-no_pie -Wno-deprecated ${CMAKE_CXX_FLAGS}")

# When Clang is ready
#SET( CMAKE_CXX_FLAGS "-O3 -g -std=c++11 -stdlib=libstdc++ ${CMAKE_CXX_FLAGS}")

FIND_PACKAGE( Boost 1.54 COMPONENTS regex program_options filesystem timer system )

SET( PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR} )
SET( CMAKE_INSTALL_PREFIX "." )
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_ROOT}/build/bin )

## Look for GMP ##
FIND_PACKAGE (GMP)

IF (GMP_FOUND)
  MESSAGE ("-- found GMP.")
ELSE()
  SET (FETCH_GMP TRUE)
  MESSAGE ("-- couldn't find GMP; fetching from internet")
ENDIF()

## Look for MPFR ##
FIND_PACKAGE (MPFR)

IF (MPFR_FOUND)
  MESSAGE ("-- found MPFR.")
ELSE()
  SET (FETCH_MPFR TRUE)
  MESSAGE ("-- couldn't find MPFR; fetching from internet")
ENDIF()

## Look for CLN ##
FIND_PACKAGE (CLN)

IF (CLN_FOUND)
  MESSAGE ("-- found CLN.")
ELSE()
  SET (FETCH_CLN TRUE)
  MESSAGE ("-- couldn't find CLN; fetching from internet")
ENDIF()

## Look for PUGIXML ##
FIND_LIBRARY ( PUGIXML_LIBRARY NAMES pugixml HINTS ${LINK_DIRECTORIES})
FIND_PATH ( PUGIXML_INCLUDE_DIR pugixml.hpp HINTS ${INCLUDE_DIRECTORIES} PATH_SUFFIXES pugixml )


## Look for Bio++ ##

FIND_LIBRARY( BPP_CORE_LIBRARY NAMES bpp-core HINTS ${LINK_DIRECTORIES})
FIND_PATH (BPP_CORE_INCLUDE_DIR BppString.h HINTS ${INCLUDE_DIRECTORIES} PATH_SUFFIXES Bpp)
FIND_LIBRARY( BPP_PHYL_LIBRARY NAMES bpp-phyl HINTS ${LINK_DIRECTORIES})
FIND_PATH (BPP_PHYL_INCLUDE_DIR Tree.h HINTS ${INCLUDE_DIRECTORIES} PATH_SUFFIXES Bpp/Phyl)
FIND_LIBRARY( BPP_SEQ_LIBRARY NAMES bpp-seq HINTS ${LINK_DIRECTORIES})
FIND_PATH (BPP_SEQ_INCLUDE_DIR Sequence.h HINTS ${INCLUDE_DIRECTORIES} PATH_SUFFIXES Bpp/Seq)

IF (BPP_CORE_LIBRARY-NOTFOUND OR BPP_CORE_INCLUDE_DIR-NOTFOUND OR
    BPP_PHYL_LIBRARY-NOTFOUND OR BPP_PHYL_INCLUDE_DIR-NOTFOUND OR
    BPP_SEQ_LIBRARY-NOTFOUND OR BPP_SEQ_INCLUDE_DIR-NOTFOUND   )
  MESSAGE ("-- couldn't find Bio++; fetching from the internet")
  SET (FETCH_BPP TRUE)
ELSE()
  MESSAGE ("-- found Bio++.")
ENDIF()

########
##
#   Fetching of external dependencies
##
########
INCLUDE (ExternalProject)

IF (PUGIXML_LIBRARY-NOTFOUND OR PUGIXML_INCLUDE_DIR-NOTFOUND )
  MESSAGE("-- Did not find pugixml library and heades")
  SET(FETCH_PUGIXML TRUE)
ELSE()
  MESSAGE("-- Found pugiXML library and headers (lib = ${PUGIXML_LIBRARY})")
ENDIF()

IF (FETCH_PUGIXML)
  EXTERNALPROJECT_ADD(
    pugixml_proj
    URL http://pugixml.googlecode.com/files/pugixml-1.2.zip
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/pugixml
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/pugixml/build    
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    PATCH_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build && cp <SOURCE_DIR>/scripts/CMakeLists.txt <SOURCE_DIR>"
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=g++ -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    BUILD_COMMAND sh -c "make"
    INSTALL_COMMAND sh -c "mkdir -p <INSTALL_DIR>/include && cp <SOURCE_DIR>/src/*.hpp <INSTALL_DIR>/include"
    BUILD_IN_SOURCE 0
  )
ENDIF()


IF (FETCH_GMP)
  EXTERNALPROJECT_ADD(
    gmp_proj
    URL ftp://ftp.gnu.org/gnu/gmp/gmp-5.1.2.tar.bz2
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/gmp
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/gmp
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND sh -c "CC=g++ CXX=g++ ./configure --prefix=<INSTALL_DIR>"
    BUILD_COMMAND sh -c "make -j4"
    INSTALL_COMMAND sh -c "make install"
  )
ENDIF()


IF (FETCH_MPFR)
  ## If we fetched GMP, use the one we fetched
  IF (FETCH_GMP)
    SET (GMP_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    SET (MPFR_CONFIGURE_STRING "CC=g++ CXX=g++ ./configure --prefix=${GMP_INSTALL_DIR} --with-gmp=${GMP_INSTALL_DIR}")
  ELSE() # Otherwise use the system GMP
    SET (MPFR_CONFIGURE_STRING "CC=g++ CXX=g++ ./configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/external/install")    
  ENDIF()

  EXTERNALPROJECT_ADD(
    mpfr_proj
    URL http://www.mpfr.org/mpfr-current/mpfr-3.1.2.tar.bz2
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/mpfr
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/mpfr
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND sh -c ${MPFR_CONFIGURE_STRING}
    BUILD_COMMAND sh -c "make -j4"
    INSTALL_COMMAND sh -c "make install"
  )
ENDIF()

IF (FETCH_CLN)
  ## If we fetched GMP, use the one we fetched
  IF (FETCH_GMP)
    SET (GMP_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    SET (CLN_CONFIGURE_STRING "CC=g++ CXX=g++ ./configure --prefix=${GMP_INSTALL_DIR} --with-gmp=${GMP_INSTALL_DIR}")
  ELSE() # Otherwise use the system GMP
    SET (CLN_CONFIGURE_STRING "CC=g++ CXX=g++ ./configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/external/install")    
  ENDIF()


  EXTERNALPROJECT_ADD(
    cln_proj
    URL http://www.ginac.de/CLN/cln-1.3.3.tar.bz2
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cln
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cln
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND sh -c ${MPFR_CONFIGURE_STRING}
    BUILD_COMMAND sh -c "make -j4"
    INSTALL_COMMAND sh -c "make install"
  )
ENDIF()


IF (FETCH_BPP)

  EXTERNALPROJECT_ADD(
    bpp_core_proj
    URL http://biopp.univ-montp2.fr/repos/sources/bpp-core-2.1.0.tar.gz
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bpp-core
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bpp-core/build    
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    PATCH_COMMAND sh -c "mkdir -p <BINARY_DIR>"
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    BUILD_COMMAND sh -c "make -j4"
    INSTALL_COMMAND sh -c "make install"
    BUILD_IN_SOURCE 0
  )

  EXTERNALPROJECT_ADD(
    bpp_seq_proj
    DEPENDS bpp_core_proj
    URL http://biopp.univ-montp2.fr/repos/sources/bpp-seq-2.1.0.tar.gz
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bpp-seq
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bpp-seq/build    
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    PATCH_COMMAND sh -c "mkdir -p <BINARY_DIR>"
    CMAKE_ARGS -Dbpp-core_INCLUDE_DIR=<INSTALL_DIR>/include -Dbpp-core_LIBRARY=<INSTALL_DIR>/lib/libbpp-core.${SHARED_LIB_EXT}
               -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    BUILD_COMMAND sh -c "make -j4"
    INSTALL_COMMAND sh -c "make install"
    BUILD_IN_SOURCE 0
  )

  EXTERNALPROJECT_ADD(
    bpp_phyl_proj
    DEPENDS bpp_seq_proj
    URL http://biopp.univ-montp2.fr/repos/sources/bpp-phyl-2.1.0.tar.gz
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bpp-phyl
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bpp-phyl/build    
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    PATCH_COMMAND sh -c "mkdir -p <BINARY_DIR>"
    CMAKE_ARGS -Dbpp-core_INCLUDE_DIR=<INSTALL_DIR>/include -Dbpp-core_LIBRARY=<INSTALL_DIR>/lib/libbpp-core.${SHARED_LIB_EXT}
               -Dbpp-seq_INCLUDE_DIR=<INSTALL_DIR>/include -Dbpp-seq_LIBRARY=<INSTALL_DIR>/lib/libbpp-seq.${SHARED_LIB_EXT}
               -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    BUILD_COMMAND sh -c "make -j4"
    INSTALL_COMMAND sh -c "make install"
    BUILD_IN_SOURCE 0
  )

ENDIF()


########
##
#   END: Fetching of external dependencies
##
########


SUBDIRS( src )

set(CPACK_SET_DESTDIR ON)
SET(CPACK_PACKAGING_INSTALL_PREFIX ".")
SET(CPACK_PACKAGE_DEFAULT_LOCATION "${CPACK_PACKAGE_NAME}")
set(CPACK_PACKAGE_EXECUTABLES "parana2" "Parana2")
set(CPACK_PACKAGE_NAME "Parana2")
set(CPACK_PACKAGE_VENDOR "Rob Patro")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Ensemble method for network history inference")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "9")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_SOURCE_GENERATOR "TBZ2")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-src")
set(CPACK_SOURCE_IGNORE_FILES
  "/doc/;/.git/;.gitignore;/build/;/.bzr/;~$;${CPACK_SOURCE_IGNORE_FILES}")
INCLUDE( CPack )
