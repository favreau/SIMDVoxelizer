# Copyright (c) 2015-2017, EPFL/Blue Brain Project
# All rights reserved. Do not distribute without permission.
# Responsible Author: Cyrille Favreau <cyrille.favreau@epfl.ch>
#
# This file is part of SIMDVoxelizer <https://github.com/favreau/SIMDVoxelizer>


# Configure the output directories. To allow IMPI to do its magic we
# will put *executables* into the (same) build directory, but tag
# mic-executables with ".mic". *libraries* cannot use the
# ".mic"-suffix trick, so we'll put libraries into separate
# directories (names 'intel64' and 'mic', respectively)
macro(CONFIGURE_ISPC)
  set(SIMDVOXELIZER_EXE_SUFFIX "")
  set(SIMDVOXELIZER_LIB_SUFFIX "")
  set(SIMDVOXELIZER_ISPC_SUFFIX ".o")
  set(THIS_IS_MIC OFF)
  set(__XEON__ ON)
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    include(${PROJECT_SOURCE_DIR}/CMake/icc.cmake)
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    include(${PROJECT_SOURCE_DIR}/CMake/gcc.cmake)
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
    include(${PROJECT_SOURCE_DIR}/CMake/clang.cmake)
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    include(${PROJECT_SOURCE_DIR}/CMake/clang.cmake)
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
    include(${PROJECT_SOURCE_DIR}/CMake/msvc.cmake)
  else()
    message(FATAL_ERROR "Unsupported compiler specified: '${CMAKE_CXX_COMPILER_ID}'")
  endif()

  if(SIMDVOXELIZER_BUILD_ISA STREQUAL "ALL")
    set(SIMDVOXELIZER_ISPC_TARGET_LIST sse4 avx avx2)
    if(ISPC_BUILD_ENABLE_KNL)
      set(SIMDVOXELIZER_ISPC_TARGET_LIST sse4 avx avx2 avx512knl-i32x16)
    endif()
  elseif(ISPC_BUILD_ISA STREQUAL "AVX512")
    set(SIMDVOXELIZER_ISPC_TARGET_LIST avx512knl-i32x16)
  elseif(ISPC_BUILD_ISA STREQUAL "AVX2")
    set(SIMDVOXELIZER_ISPC_TARGET_LIST avx2)
  elseif(ISPC_BUILD_ISA STREQUAL "AVX")
    set(SIMDVOXELIZER_ISPC_TARGET_LIST avx)
  elseif(ISPC_BUILD_ISA STREQUAL "SSE")
    set(SIMDVOXELIZER_ISPC_TARGET_LIST sse4)
  else ()
    message(ERROR "Invalid ISPC_BUILD_ISA value. Please select one of SSE, AVX, AVX2, or ALL.")
  endif()

  if(NOT ISPC_COMPILER_SUPPORTS_AVX)
    if(NOT ISPC_WARNED_MISSING_AVX)
      message("Warning: Need at least version ${GCC_VERSION_REQUIRED_AVX} of gcc for AVX. Disabling AVX.\nTo compile for AVX, please switch to either 'ICC'-compiler, or upgrade your gcc version.")
      set(ISPC_WARNED_MISSING_AVX ON CACHE INTERNAL "Warned about missing AVX support.")
    endif()
  endif()

  if(NOT ISPC_COMPILER_SUPPORTS_AVX2)
    if(NOT ISPC_WARNED_MISSING_AVX2)
      message("Warning: Need at least version ${GCC_VERSION_REQUIRED_AVX2} of gcc for AVX2. Disabling AVX2.\nTo compile for AVX2, please switch to either 'ICC'-compiler, or upgrade your gcc version.")
      set(ISPC_WARNED_MISSING_AVX2 ON CACHE INTERNAL "Warned about missing AVX2 support.")
    endif()
  endif()

  if(NOT ISPC_COMPILER_SUPPORTS_AVX512)
    if(NOT ISPC_WARNED_MISSING_AVX2)
      message("Warning: Need at least version ${GCC_VERSION_REQUIRED_AVX512} of gcc for AVX512. Disabling AVX512.\nTo compile for AVX512, please switch to either 'ICC'-compiler, or upgrade your gcc version.")
      set(ISPC_WARNED_MISSING_AVX512 ON CACHE INTERNAL "Warned about missing AVX512 support.")
    endif()
  endif()

  if(THIS_IS_MIC)
    # whether to build in MIC/xeon phi support
    set(ISPC_BUILD_COI_DEVICE OFF CACHE BOOL "Build COI Device for MIC support?")
  endif()

  include(${PROJECT_SOURCE_DIR}/CMake/ispc.cmake)

  include_directories(${PROJECT_SOURCE_DIR})

  INCLUDE_DIRECTORIES_ISPC(${PROJECT_SOURCE_DIR})

  # for auto-generated cmakeconfig etc
  include_directories(${PROJECT_BINARY_DIR})
  INCLUDE_DIRECTORIES_ISPC(${PROJECT_BINARY_DIR})
endmacro()

set(SIMDVOXELIZER_SOURCES
  ispc/tasksys.cpp
  OctreeNode.cpp
  Octree.cpp
)

set(SIMDVOXELIZER_PUBLIC_HEADERS
  OctreeNode.h
  Octree.h
)

set(SIMDVOXELIZER_ISPC_SOURCES
  ispc/SIMDSparseVoxelizer.ispc
)

find_package(glm REQUIRED SYSTEM)

include(ispc)
CONFIGURE_ISPC()
INCLUDE_DIRECTORIES_ISPC(${SIMDVOXELIZER_INCLUDE_DIRS})
ISPC_COMPILE(${SIMDVOXELIZER_ISPC_SOURCES})
list(APPEND SIMDVOXELIZER_SOURCES ${ISPC_OBJECTS})

include_directories(${GLM_INCLUDE_DIR})

ADD_LIBRARY(
    simdvoxelizer
    ${SIMDVOXELIZER_SOURCES})

TARGET_LINK_LIBRARIES(
    simdvoxelizer
    ${GLM_LIBRARIES} pthread stdc++)
