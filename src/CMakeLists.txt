INCLUDE(CheckLibraryExists)

SET( PARANA_SRCS CrossValidationParser.cpp FlipKey.cpp TreeUtils.cpp PhyloXMLParser.cpp MultiOpt.cpp main.cpp arg_parser.cpp)

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/external/install/include      
  ${Boost_INCLUDE_DIR}
  /usr/local/include
  /opt/local/include
  ${PROJECT_ROOT}/include/
  ${CMAKE_SOURCE_DIR}
  ./
)

MESSAGE(STATUS "BOOST_LIBRARYDIR = ${BOOST_LIBRARYDIR}")
MESSAGE(STATUS "Boost Libraries = ${Boost_LIBRARIES}")

LINK_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/external/install/lib  
  ${BOOST_LIBRARYDIR}
  /usr/lib
  /opt/local/lib
  /usr/local/lib
)

FIND_LIBRARY( BPP_CORE_LIBRARY NAMES bpp-core HINTS ${LINK_DIRECTORIES})
FIND_LIBRARY( BPP_PHYL_LIBRARY NAMES bpp-phyl HINTS ${LINK_DIRECTORIES})
FIND_LIBRARY( BPP_SEQ_LIBRARY NAMES bpp-seq HINTS ${LINK_DIRECTORIES})

FIND_LIBRARY( PUGIXML_LIBRARY NAMES pugixml HINTS ${LINK_DIRECTORIES})
FIND_LIBRARY( GMP_LIBRARY NAMES gmp HINTS ${LINK_DIRECTORIES})
FIND_LIBRARY( MPFR_LIBRARY NAMES mpfr HINTS ${LINK_DIRECTORIES})
FIND_LIBRARY( CLN_LIBRARY NAMES cln HINTS ${LINK_DIRECTORIES})

SET( BPP_LIBRARIES ${BPP_CORE_LIBRARY} ${BPP_PHYL_LIBRARY} ${BPP_SEQ_LIBRARY} )
SET( DEPENDENT_LIBS ${Boost_LIBRARIES} ${BPP_LIBRARIES} ${PUGIXML_LIBRARY} 
                    ${GMP_LIBRARY} ${MPFR_LIBRARY} ${CLN_LIBRARY} )

ADD_EXECUTABLE( parana2 ${PARANA_SRCS} )
TARGET_LINK_LIBRARIES( parana2 m pthread ${DEPENDENT_LIBS} )

## scripts
FILE( GLOB SCRIPTS ${PROJECT_ROOT}/scripts/*.py ${PROJECT_ROOT}/scripts/*.sh )
FILE( GLOB CVFILES ${PROJECT_ROOT}/data/cvfiles/*.xml )

INSTALL( PROGRAMS ${EXECUTABLE_OUTPUT_PATH}/parana2 DESTINATION bin )
INSTALL( PROGRAMS ${PROJECT_ROOT}/RunParana2.sh DESTINATION bin )

# The below code replaces the following command:
# INSTALL( FILES ${DEPENDENT_LIBS} DESTINATION lib )
#
# It chases down and copies the appropriate symlinks

set (_resolvedFiles "")
foreach (_file ${DEPENDENT_LIBS})
    get_filename_component(_resolvedFile "${_file}" REALPATH)
    list (APPEND _resolvedFiles "${_resolvedFile}")
endforeach()
install(FILES ${_resolvedFiles} DESTINATION lib)


#INSTALL( FILES ${SCRIPTS} DESTINATION scripts )
INSTALL( FILES ${PROJECT_ROOT}/data/SimpleUndirectedModel.xml
               ${PROJECT_ROOT}/data/UndirectedParsimonyModel.xml 
               DESTINATION data )
INSTALL( FILES ${CVFILES} DESTINATION data/cvfiles )
