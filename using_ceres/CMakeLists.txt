CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(simple_snapping)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(EIGEN_INCLUDE_DIR "/usr/local/include/eigen3")
MESSAGE(${EIGEN_INCLUDE_DIR})

FIND_PACKAGE(Ceres REQUIRED)
INCLUDE_DIRECTORIES(${CERES_INCLUDE_DIRS})


FIND_PACKAGE(GEOS REQUIRED)
INCLUDE_DIRECTORIES(${GEOS_INCLUDE_DIR})

# simple_snapping

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -Wno-long-long -pedantic -O3") #-O3

#SET(CMAKE_BUILD_TYPE "Debug")
#SET(CMAKE_BUILD_TYPE "Release")

MESSAGE(${CERES_LIBRARIES})
MESSAGE(${GEOS_LIBRARY})
MESSAGE(${GEOS_INCLUDE_DIR})

#adding a flag to allow use of unordered maps
ADD_EXECUTABLE(simple_snapping.bin
    Constraints.cc
    geometry_function.cc
    enum_functions.cc
    Parameters.cc
    Data.cc
    simple_snapping.cc
    )
TARGET_LINK_LIBRARIES(simple_snapping.bin ${CERES_LIBRARIES} ${GEOS_LIBRARY})
