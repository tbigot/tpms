cmake_minimum_required(VERSION 2.6)
PROJECT(tpms)

add_definitions(-std=c++11)

set(EXECUTABLE_OUTPUT_PATH build/${CMAKE_BUILD_TYPE})

IF(NOT BPP_LIBRARY)
    FIND_PATH(BPP_INCLUDE_DIR "Bpp/Phyl/TreeTemplate.h")
    SET(BPP_NAMES bpp-phyl bpp-phyllib bpp-phyldll)
    FIND_LIBRARY(BPP_LIBRARY NAMES ${BPP_NAMES})
ENDIF(NOT BPP_LIBRARY)

IF(BPP_LIBRARY)
    MESSAGE("-- Library BPP found here:")
    MESSAGE("   includes : ${BPP_INCLUDE_DIR}")
    MESSAGE("   libraries: ${BPP_LIBRARY}")
ELSE(BPP_LIBRARY)
    MESSAGE(FATAL_ERROR "bpp-phyl required but not found.")
ENDIF(BPP_LIBRARY)


link_directories(${BPP_LIBRARY})
include_directories(${BPP_INCLUDE_DIR})

# Autres packages à trouver : -lbpp-utils -lbpp-numcalc -lbpp-seq

#boost
# search for Boost version 1.40
# Components :
 #filesystem, iostreams, programoptions, python, regex, serialization, signals
#system, thread, wave

#SET(Boost_ADDITIONAL_VERSIONS "1.42" "1.42.0" "1.43.0")
find_package(Boost COMPONENTS filesystem system thread REQUIRED)

link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} )



file(
  GLOB
  classes
  src/classes/*
)


add_executable(tpms_query src/tpms_query.cpp ${classes})
target_link_libraries(tpms_query bpp-phyl bpp-core bpp-seq ${Boost_LIBRARIES})

#STATIC
add_executable(tpms_query_static src/tpms_query.cpp ${classes})
target_link_libraries(tpms_query_static bpp-phyl bpp-core bpp-seq -static)

# add_executable(tpms_jump src/tpms_jump.cpp ${classes})
# target_link_libraries(tpms_jump bpp-phyl bpp-utils bpp-numcalc bpp-seq)

add_executable(tpms_mkdb src/tpms_mkdb.cpp ${classes})
target_link_libraries(tpms_mkdb bpp-phyl bpp-core bpp-seq bpp-raa ${Boost_LIBRARIES})

add_executable(tpms_computations src/tpms_computations.cpp ${classes})
target_link_libraries(tpms_computations bpp-phyl bpp-core bpp-seq bpp-raa ${Boost_LIBRARIES})


#add_executable(tpms_mkdb_static src/tpms_mkdb.cpp ${classes})
#target_link_libraries(STATIC tpms_mkdb_static bpp-phyl bpp-utils bpp-numcalc bpp-seq bpp-raa ${Boost_LIBRARIES} z)
