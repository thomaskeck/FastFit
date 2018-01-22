cmake_minimum_required (VERSION 2.91.0)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

project (FastFit)
set (FastFit_VERSION_MAJOR 1)
set (FastFit_VERSION_MINOR 1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3 -std=c++11 -Wall -g -msse2")

configure_file (
  "${PROJECT_SOURCE_DIR}/include/FastFit.h.in"
  "${PROJECT_BINARY_DIR}/include/FastFit.h"
)

find_program(BASF2 "basf2")

if (BASF2)
    set(CMAKE_INSTALL_PREFIX "$ENV{BELLE2_EXTERNALS_DIR}/$ENV{BELLE2_EXTERNALS_SUBDIR}")
    set(EIGEN3_INCLUDE_DIR "$ENV{BELLE2_EXTERNALS_DIR}/include/")
endif()

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

set(FastFit_CINTERFACE
    "${PROJECT_SOURCE_DIR}/src/FastFit_C_API.cc"
    "${PROJECT_SOURCE_DIR}/include/FastFit_C_API.h" 
)

set(FastFit_Python
  "${PROJECT_BINARY_DIR}/PyFastFit/__init__.py"
  "${PROJECT_BINARY_DIR}/PyFastFit/FastFit.py"
)

add_library(FastFit_static STATIC ${FastFit_SOURCES} ${FastFit_HEADERS})
target_link_libraries(FastFit_static)

add_library(FastFit_shared SHARED ${FastFit_SOURCES} ${FastFit_HEADERS})
target_link_libraries(FastFit_shared)

add_library(FastFit_CInterface SHARED ${FastFit_CINTERFACE} ${FastFit_SOURCES} ${FastFit_HEADERS})
target_link_libraries(FastFit_CInterface)

install(TARGETS FastFit_static FastFit_shared FastFit_CInterface 
     LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib
      RUNTIME DESTINATION bin   
   )

 install(FILES ${FastFit_HEADERS} DESTINATION include)

if(BASF2)
  execute_process(COMMAND patch -p1 -d $ENV{BELLE2_LOCAL_DIR} -i ${PROJECT_SOURCE_DIR}/files/Belle2.patch --dry-run RESULT_VARIABLE canpatch)
  if(${canpatch} STREQUAL 0)
    install(CODE "execute_process(COMMAND patch -p1 -d $ENV{BELLE2_LOCAL_DIR} -i ${PROJECT_SOURCE_DIR}/files/Belle2.patch)")
    message(STATUS "Patching basf2 particle vertex fitter will work -> applying patch when installing.")
  else()
    message(STATUS "Cannot patch basf2 particle vertex fitter -> won't do it.")
  endif()
endif()

if(GTEST_FOUND)
  add_executable(unittests ${FastFit_TESTS} ${FastFit_HEADERS})
  target_link_libraries(unittests ${GTEST_BOTH_LIBRARIES} FastFit_static pthread)
  message(STATUS  ${GTEST_INCLUDE_DIRS})
  target_include_directories(unittests PUBLIC BEFORE ${GTEST_INCLUDE_DIRS})
else()
  message(STATUS "Could not find gtest installation, skip building unittests.")
endif()

find_program(PYTHON "python3")

if (PYTHON)
    configure_file(
        "${PROJECT_SOURCE_DIR}/setup.py.in"
        "${PROJECT_SOURCE_DIR}/setup.py"
    )

    add_custom_target(PyFastFit ALL DEPENDS ${FastFit_Python} FastFit_shared FastFit_CInterface)

    add_custom_command(TARGET PyFastFit PRE_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_SOURCE_DIR}/PyFastFit" "${CMAKE_BINARY_DIR}/PyFastFit"
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:FastFit_shared> "${PROJECT_BINARY_DIR}/PyFastFit/"
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:FastFit_CInterface> "${PROJECT_BINARY_DIR}/PyFastFit/"
                       COMMAND ${PYTHON} "${CMAKE_BINARY_DIR}/setup.py" build
    )

    install(CODE "execute_process(COMMAND ${PYTHON} ${PROJECT_BINARY_DIR}/setup.py install --prefix=${CMAKE_INSTALL_PREFIX})")

endif()

set(CPACK_PACKAGE_VERSION "${FastFit_VERSION_MAJOR}.${FastFit_VERSION_MINOR}")
set(CPACK_GENERATOR "RPM;DEB;TGZ")
set(CPACK_PACKAGE_NAME "FastFit")
set(CPACK_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_CONTACT "thomas.keck2@kit.edu")
set(CPACK_PACKAGE_VENDOR "Private")
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")

SET(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
SET(CPACK_DEBIAN_PACKAGE_SECTION "libs")
SET(CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})

include(CPack)
