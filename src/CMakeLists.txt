INCLUDE(CheckLibraryExists)

SET( PARANA_SRCS main.cpp arg_parser.cpp )

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )

INCLUDE_DIRECTORIES(
  ${BOOST_INCLUDEDIR}
  /usr/local/include
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

LINK_LIBRARIES( m pthread ${Boost_LIBRARIES} profiler bpp-core
bpp-phyl bpp-seq cln)

ADD_EXECUTABLE( parana++ ${PARANA_SRCS} )
