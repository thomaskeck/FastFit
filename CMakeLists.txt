cmake_minimum_required (VERSION 2.91.0)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

project (FastFit)
set (FastFit_VERSION_MAJOR 3)
set (FastFit_VERSION_MINOR 4)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3 -std=c++11 -Wall -Wextra -g -msse2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -ggdb3 -g -std=c++11 -Wall -Wextra")

configure_file (
  "${PROJECT_SOURCE_DIR}/include/FastFit.h.in"
  "${PROJECT_BINARY_DIR}/include/FastFit.h"
)

set(EIGEN3_INCLUDE_DIR "/local/scratch/ssd/tkeck/externals/development/include/")

find_package(GTest)
find_package( Eigen3 REQUIRED )
include_directories("${EIGEN3_INCLUDE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/include/")
include_directories(BEFORE ${GTEST_INCLUDE_DIRS})

set(FastFit_SOURCES
  "${PROJECT_SOURCE_DIR}/src/FastFit.cc"
)

set(FastFit_TESTS
  "${PROJECT_SOURCE_DIR}/src/test_FastFit.cc"
)

set(FastFit_HEADERS
  "${PROJECT_BINARY_DIR}/include/FastFit.h"
)

add_library(FastFit_static STATIC ${FastFit_SOURCES} ${FastFit_HEADERS})
target_link_libraries(FastFit_static)
add_library(FastFit_shared SHARED ${FastFit_SOURCES} ${FastFit_HEADERS})
target_link_libraries(FastFit_shared)

install(TARGETS FastFit_static FastFit_shared
     LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib
      RUNTIME DESTINATION bin   
   )

if(GTEST_FOUND)
  add_executable(unittests ${FastFit_TESTS} ${FastFit_HEADERS})
  target_link_libraries(unittests ${GTEST_BOTH_LIBRARIES} FastFit_static pthread)
  message(STATUS  ${GTEST_INCLUDE_DIRS})
  target_include_directories(unittests PUBLIC BEFORE ${GTEST_INCLUDE_DIRS})
else()
  message(STATUS "Could not find gtest installation, skip building unittests.")
endif()

