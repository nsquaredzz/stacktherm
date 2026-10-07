cmake_minimum_required(VERSION 3.20)
project(stacktherm VERSION 0.2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(STACKTHERM_PYTHON "Build the Python bindings used by the studio and the kALDo bridge" ON)
option(STACKTHERM_TESTS "Build the tests" ON)

# ---- dependencies: AMGCL (header-only multigrid) and yaml-cpp, fetched unless
# already present in third_party/
include(FetchContent)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
if(EXISTS ${CMAKE_SOURCE_DIR}/third_party/amgcl/amgcl/amg.hpp)
  set(AMGCL_INCLUDE ${CMAKE_SOURCE_DIR}/third_party/amgcl)
else()
  FetchContent_Declare(amgcl GIT_REPOSITORY https://github.com/ddemidov/amgcl.git GIT_TAG master
                       GIT_SHALLOW ON SOURCE_SUBDIR this-directory-does-not-exist)
  FetchContent_MakeAvailable(amgcl)
  set(AMGCL_INCLUDE ${amgcl_SOURCE_DIR})
endif()
if(EXISTS ${CMAKE_SOURCE_DIR}/third_party/yaml-cpp/CMakeLists.txt)
  add_subdirectory(third_party/yaml-cpp EXCLUDE_FROM_ALL)
else()
  FetchContent_Declare(yaml-cpp GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git GIT_TAG 0.8.0
                       GIT_SHALLOW ON EXCLUDE_FROM_ALL)
  FetchContent_MakeAvailable(yaml-cpp)
endif()
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# ---- core library
file(GLOB STACKTHERM_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/cpp/src/*.cpp)
add_library(stacktherm_core STATIC ${STACKTHERM_SOURCES})
target_include_directories(stacktherm_core PUBLIC ${CMAKE_SOURCE_DIR}/cpp/include)
target_include_directories(stacktherm_core SYSTEM PRIVATE ${AMGCL_INCLUDE})
target_compile_definitions(stacktherm_core PRIVATE AMGCL_NO_BOOST
                           STACKTHERM_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
target_link_libraries(stacktherm_core PRIVATE yaml-cpp::yaml-cpp ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
if(APPLE)
  target_compile_definitions(stacktherm_core PRIVATE ACCELERATE_NEW_LAPACK)
endif()
target_compile_options(stacktherm_core PRIVATE -Wall -Wextra -Wno-unused-parameter)

# ---- command line
if(EXISTS ${CMAKE_SOURCE_DIR}/cpp/cli/main.cpp)
  add_executable(stacktherm cpp/cli/main.cpp)
  target_link_libraries(stacktherm PRIVATE stacktherm_core)
endif()

# ---- tests
if(STACKTHERM_TESTS)
  enable_testing()
  file(GLOB STACKTHERM_TESTS_SRC CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/cpp/tests/test_*.cpp)
  foreach(src ${STACKTHERM_TESTS_SRC})
    get_filename_component(name ${src} NAME_WE)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE stacktherm_core)
    target_compile_definitions(${name} PRIVATE STACKTHERM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    add_test(NAME ${name} COMMAND ${name})
  endforeach()
endif()

# ---- Python bindings
if(STACKTHERM_PYTHON AND EXISTS ${CMAKE_SOURCE_DIR}/cpp/python/bindings.cpp)
  set(PYBIND11_FINDPYTHON ON)
  find_package(pybind11 CONFIG QUIET)
  if(pybind11_FOUND)
    pybind11_add_module(_core cpp/python/bindings.cpp)
    target_link_libraries(_core PRIVATE stacktherm_core)
    set_target_properties(_core PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/stacktherm)
  else()
    message(STATUS "pybind11 not found: Python bindings are skipped")
  endif()
endif()
