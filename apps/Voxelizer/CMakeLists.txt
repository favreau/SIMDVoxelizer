# Copyright (c) 2015-2017, EPFL/Blue Brain Project
# All rights reserved. Do not distribute without permission.
# Responsible Author: Cyrille Favreau <cyrille.favreau@epfl.ch>
#
# This file is part of SIMDVoxelizer <https://github.com/favreau/SIMDVoxelizer>
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License version 3.0 as published
# by the Free Software Foundation.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

set(PROJECT_NAME Voxelizer)
set(VOXELIZER_SOURCES main.cpp)
set(VOXELIZER_LINK_LIBRARIES PUBLIC simdvoxelizer)

include_directories(${PROJECT_SOURCE_DIR} ${CMAKE_BINARY_DIR})

ADD_EXECUTABLE(
    ${PROJECT_NAME}
    ${VOXELIZER_SOURCES}
    )

TARGET_LINK_LIBRARIES(
    ${PROJECT_NAME}
    ${VOXELIZER_LINK_LIBRARIES}
)

# ------------------------------------------------------------
INSTALL(TARGETS ${PROJECT_NAME} DESTINATION bin)
# ------------------------------------------------------------
