CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(simple_snapping)

FIND_PACKAGE(Ceres REQUIRED)
INCLUDE_DIRECTORIES(${CERES_INCLUDE_DIRS})

# simple_snapping

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#adding a flag to allow use of unordered maps
ADD_EXECUTABLE(simple_snapping.bin 
    simple_snapping.cc
    Data.cc)
TARGET_LINK_LIBRARIES(simple_snapping.bin ${CERES_LIBRARIES})
