
cmake_minimum_required(VERSION 3.10)

project(quant3d VERSION 1.0)

### Default to C++11, CLI and tests need C++17, however
set (CMAKE_CXX_STANDARD 17)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build this in a source directory (or any directory with a CMakeLists.txt file). "
                        "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

### Build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
  "Debug" "Release" "RelWithDebInfo")
endif()

### Cmake package finder files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(macrOS)

if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(quant3d_STANDALONE TRUE)
else()
    set(quant3d_STANDALONE FALSE)
endif()

# ### Binary output
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Build/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Build/lib/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Build/bin/)

### Options
option(quant3d_ENABLE_LOGGING    "Enable logging for quant3d" ${quant3d_STANDALONE})
option(quant3d_BUILD_CLI         "Build CLI app for quant3d"  ${quant3d_STANDALONE})
option(quant3d_SUBMODULES_MANUAL "Skip automatic submodule download" OFF)
set(quant3d_SUBMODULES_MANUAL ON)

### Load external libs
add_subdirectory(extern)

### Add library
set(quant3d_LIB_NAME quant3d)
add_subdirectory(src)

### Add CLI
if (quant3d_BUILD_CLI)
    add_subdirectory(cli)
endif()
