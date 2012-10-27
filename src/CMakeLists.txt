INCLUDE(CheckLibraryExists)

SET( PARANA_SRCS FlipKey.cpp TreeUtils.cpp MultiOpt.cpp main.cpp arg_parser.cpp)

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )

INCLUDE_DIRECTORIES(
  ${BOOST_INCLUDEDIR}
  /usr/local/include
  /opt/local/include
  /opt/local/include/gcc47/c++
  /opt/local/include/gcc47/c++/x86_64-apple-darwin11
  ${PROJECT_ROOT}/include/
  ${CMAKE_SOURCE_DIR}
  ./
)


LINK_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/include
  ${BOOST_LIBRARYDIR}
  /usr/lib
  /opt/local/lib
  /usr/local/lib
  /opt/local/lib/gcc47
  ${PROJECT_ROOT}/lib
)

LINK_LIBRARIES( m pthread ${Boost_LIBRARIES} bpp-core bpp-phyl bpp-seq
gmp mpfr cln MPFRPP profiler )

ADD_EXECUTABLE( parana++ ${PARANA_SRCS} )
