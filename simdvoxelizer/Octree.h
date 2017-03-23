/* Copyright (c) 2015-2017, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Cyrille Favreau <cyrille.favreau@epfl.ch>
 *                     Grigori Chevtchenko <grigori.chevtchenko@epfl.ch>
 *
 * This file is part of SIMDVoxelizer <https://github.com/favreau/SIMDVoxelizer>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef _Octree_h_
#define _Octree_h_

#include <vector>
#include <map>
#include <glm/glm.hpp>

#include "OctreeNode.h"

typedef std::map< uint64_t, OctreeNode > OctreeLevelMap;

class Octree
{
public:
    Octree( const std::vector<Event>& events, float voxelSize, const glm::vec3& min, const glm::vec3& max );
    ~Octree();

    const OctreeNode* getRoot() const;
    uint32_t getOctreeSize() const;

private:
    void _flattenChildren( const OctreeNode* node, uint32_t level );

    inline int32_t _pow2roundup( int32_t x )
    {
        if (x < 0)
            return 0;
        --x;
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        return x+1;
    }

    uint32_t _octreeSize;
    uint32_t _octreeDepth;

    std::vector< OctreeLevelMap > _octree;
};

#endif //_Octree_h_
