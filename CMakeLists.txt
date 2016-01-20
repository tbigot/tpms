cmake_minimum_required(VERSION 2.8)
PROJECT(tpms)

add_definitions(-std=gnu++0x)

set(EXECUTABLE_OUTPUT_PATH build/${CMAKE_BUILD_TYPE})

#here is a useful function:
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
  IF(${BUILD_STATIC})
    SET(${lib_name}_NAMES lib${lib_name}.a)
  ELSE(${BUILD_STATIC})
    SET(${lib_name}_NAMES ${lib_name} ${lib_name}lib ${lib_name}dll)
  ENDIF(${BUILD_STATIC})
  FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})

  IF(${lib_name}_LIBRARY)
    MESSAGE("-- Library ${lib_name} ${lib_type} found here:")
    MESSAGE("   includes : ${${lib_name}_INCLUDE_DIR}")
    MESSAGE("   library: ${${lib_name}_LIBRARY}")
  ELSE(${lib_name}_LIBRARY)
    MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
  ENDIF(${lib_name}_LIBRARY)
  
  #add the dependency:
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(IMPROVED_FIND_LIBRARY)

IMPROVED_FIND_LIBRARY(BPP_LIBRARIES bpp-phyl Bpp/Phyl/Tree.h)
IMPROVED_FIND_LIBRARY(BPP_LIBRARIES bpp-seq Bpp/Seq/Alphabet/Alphabet.h)
IMPROVED_FIND_LIBRARY(BPP_LIBRARIES bpp-core Bpp/Clonable.h)



link_directories(${BPP_LIBRARIES})

# Autres packages à trouver : -lbpp-utils -lbpp-numcalc -lbpp-seq

#boost
# search for Boost version 1.40
# Components :
 #filesystem, iostreams, programoptions, python, regex, serialization, signals
#system, thread, wave

#SET(Boost_ADDITIONAL_VERSIONS "1.42" "1.42.0" "1.43.0")
find_package(Boost COMPONENTS filesystem system thread REQUIRED)




file(
  GLOB
  classes
  src/classes/*
)


add_executable(tpms_query src/tpms_query.cpp ${classes})
target_link_libraries(tpms_query bpp-phyl bpp-core bpp-seq ${Boost_LIBRARIES})

add_executable(tpms_mkdb src/tpms_mkdb.cpp ${classes})
target_link_libraries(tpms_mkdb bpp-phyl bpp-core bpp-seq bpp-raa ${Boost_LIBRARIES})

add_executable(tpms_computations src/tpms_computations.cpp ${classes})
target_link_libraries(tpms_computations bpp-phyl bpp-core bpp-seq bpp-raa ${Boost_LIBRARIES})
