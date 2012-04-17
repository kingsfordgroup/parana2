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
  /usr/local/lib
  /opt/local/lib
  ${PROJECT_ROOT}/lib
)

LINK_LIBRARIES( m pthread bfd ${Boost_LIBRARIES} profiler bpp-core bpp-phyl bpp-seq
gmp mpfr cln MPFRPP)

ADD_EXECUTABLE( parana++ ${PARANA_SRCS} )
