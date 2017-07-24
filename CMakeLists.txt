# Copyright (c) 2015-2017, EPFL/Blue Brain Project
# All rights reserved. Do not distribute without permission.
# Responsible Author: Cyrille Favreau <cyrille.favreau@epfl.ch>
#
# This file is part of SIMDVoxelizer <https://github.com/favreau/SIMDVoxelizer>

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
set(SIMDVoxelizer_VERSION_ABI 1)
project(SIMDVoxelizer VERSION 0.1.0)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake
                              ${CMAKE_SOURCE_DIR}/CMake/common)

set(SIMDVOXELIZER_DESCRIPTION "SIMD Voxelizer")

# Check for C++11 support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(SIMDVOXELIZER_BUILD_ISA "ALL" CACHE STRING "Target ISA (SSE,AVX,AVX2,or ALL)")
set_property(CACHE SIMDVOXELIZER_BUILD_ISA PROPERTY STRINGS ALL SSE AVX AVX2)

# compiler selection for the non-MIC part
if(APPLE)
  set(SIMDVOXELIZER_COMPILER "CLANG" CACHE STRING "compiler to use (ICC,GCC, or CLANG)")
  set_property(CACHE SIMDVOXELIZER_COMPILER PROPERTY STRINGS ICC GCC CLANG)
else(APPLE)
  set(SIMDVOXELIZER_COMPILER "GCC" CACHE STRING "compiler to use (ICC or GCC)")
  set_property(CACHE SIMDVOXELIZER_COMPILER PROPERTY STRINGS ICC GCC)
endif()

# the arch we're targeting for the non-MIC/non-xeon phi part of ospray
set(SIMDVOXELIZER_BUILD_ISA "ALL" CACHE STRING "Target ISA (SSE,AVX,AVX2,or ALL)")
set_property(CACHE SIMDVOXELIZER_BUILD_ISA PROPERTY STRINGS ALL SSE AVX AVX2)

add_subdirectory(simdvoxelizer)
add_subdirectory(apps/Voxelizer)
#add_subdirectory(apps/Mesher)

