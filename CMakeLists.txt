cmake_minimum_required(VERSION 2.8)
PROJECT(tpms CXX)

SET(CMAKE_CXX_FLAGS "-std=c++11")

SET(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -g -Wall -Weffc++ -Wshadow")
SET(CMAKE_C_FLAGS_DEBUG "-g -Wall -Weffc++ -Wshadow")
SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O3 -Wall -Weffc++ -Wshadow -DNDEBUG")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -Wall -Weffc++ -Wshadow -DNDEBUG")
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-std=c++11 -O3 -w -DNDEBUG")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3 -w -DNDEBUG")

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

IMPROVED_FIND_LIBRARY(BPP_LIBRARIES bpp-raa Bpp/Raa/RAA.h)
IMPROVED_FIND_LIBRARY(BPP_LIBRARIES bpp-phyl Bpp/Phyl/Tree.h)
IMPROVED_FIND_LIBRARY(BPP_LIBRARIES bpp-seq Bpp/Seq/Alphabet/Alphabet.h)
IMPROVED_FIND_LIBRARY(BPP_LIBRARIES bpp-core Bpp/Clonable.h)

link_directories(${BPP_LIBRARIES})
message(${BPP_LIBRARIES})


#boost
# search for Boost version 1.40
# Components :
 #filesystem, iostreams, programoptions, python, regex, serialization, signals
#system, thread, wave

#SET(Boost_ADDITIONAL_VERSIONS "1.42" "1.42.0" "1.43.0")
SET(Boost_USE_STATIC_LIBS ${BUILD_STATIC})
find_package(Boost COMPONENTS filesystem system thread REQUIRED)


file(
  GLOB
  classes
  src/classes/*
)

add_executable(tpms_query src/tpms_query.cpp ${classes})
target_link_libraries(tpms_query ${BPP_LIBRARIES} ${Boost_LIBRARIES})

add_executable(tpms_mkdb src/tpms_mkdb.cpp ${classes})
target_link_libraries(tpms_mkdb ${BPP_LIBRARIES} ${Boost_LIBRARIES} -lz)

add_executable(tpms_computations src/tpms_computations.cpp ${classes})
target_link_libraries(tpms_computations ${BPP_LIBRARIES} ${Boost_LIBRARIES})

install(TARGETS tpms_query tpms_mkdb tpms_computations DESTINATION bin)

