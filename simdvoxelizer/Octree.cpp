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

#include <iostream>
#include <string.h>

#include "Octree.h"

Octree::Octree( const std::vector< Event >& events, float leafSize,
                const glm::vec3& min, const glm::vec3& max )
{
    const glm::uvec3 octreeSize( _pow2roundup( std::ceil((max.x - min.x) / leafSize )),
                                 _pow2roundup( std::ceil((max.y - min.y) / leafSize )),
                                 _pow2roundup( std::ceil((max.z - min.z) / leafSize )));

    // This octree is always isotropic
    _octreeSize = std::max( std::max( octreeSize.x,
                                      octreeSize.y),
                                      octreeSize.z);

    _octreeDepth = std::log2( _octreeSize ) + 1u;
    _octree.resize( _octreeDepth );

    std::cout<<"Octree depth: " << _octreeDepth << std::endl;

    const int32_t leafLevel = _octreeDepth - 1u;
    for( const Event& event : events )
    {
        const uint64_t xpos = std::floor(( event.position.x - min.x ) / leafSize );
        const uint64_t ypos = std::floor(( event.position.y - min.y ) / leafSize );
        const uint64_t zpos = std::floor(( event.position.z - min.z ) / leafSize );

        const uint64_t index = xpos + ypos * (uint64_t)_octreeSize +
                               zpos * (uint64_t)_octreeSize * (uint64_t)_octreeSize;

        auto it = _octree[ leafLevel ].find( index );
        if( it == _octree[ leafLevel ].end( ))
        {
            OctreeNode* child = nullptr;
            for( int32_t level = leafLevel; level >= 0; --level )
            {
                bool newNode = false;
                const uint64_t leafFactor = std::pow( 2, leafLevel - level );
                const glm::vec3 center( leafSize * leafFactor * ( xpos / leafFactor + 0.5f ) + min.x,
                                        leafSize * leafFactor * ( ypos / leafFactor + 0.5f ) + min.y,
                                        leafSize * leafFactor * ( zpos / leafFactor + 0.5f ) + min.z );


                const uint64_t nBlock = _octreeSize / leafFactor;
                const uint64_t indexAtLevel = std::floor( xpos / leafFactor ) +
                                              nBlock * std::floor( ypos / leafFactor ) +
                                              nBlock * nBlock * std::floor( zpos / leafFactor );

                const float size = leafSize * leafFactor;

                auto itr = _octree[ level ].find( indexAtLevel );
                if( itr == _octree[ level ].end( ))
                {
                    _octree[ level ].insert(
                                OctreeLevelMap::value_type(
                                    indexAtLevel, OctreeNode( center, size )));
                    newNode = true;
                }

                if( level == leafLevel )
                {
                    _octree[ level ].at( indexAtLevel ).addEvent( event );
                }
                else
                    _octree[ level ].at( indexAtLevel ).addMaxValue( event.value );

                if(( level != leafLevel ) && ( child != nullptr ))
                    _octree[level].at( indexAtLevel ).setChild( child );

                if( newNode )
                    child = &(_octree[ level ].at( indexAtLevel ));
                else
                    child = nullptr;
            }
        }
        else
        {
            for( int32_t level = leafLevel; level >= 0; --level )
            {
                const uint64_t leafFactor = std::pow( 2, leafLevel - level );
                const uint64_t nBlock = _octreeSize / leafFactor;
                const uint64_t indexAtLevel = std::floor( xpos / leafFactor ) +
                                              nBlock * std::floor( ypos / leafFactor ) +
                                              nBlock * nBlock * std::floor( zpos / leafFactor );
                if( level == leafLevel )
                {
                    _octree[ level ].at( indexAtLevel ).addEvent( event );
                }
                else
                    _octree[ level ].at( indexAtLevel ).addMaxValue( event.value );
            }
        }
    }

    std::cout<< "Octree number of leaves: ["<<leafLevel<<"]: " << _octree[leafLevel].size() << std::endl;
}

Octree::~Octree()
{}

const OctreeNode* Octree::getRoot() const
{
    return &(_octree[ 0 ].at( 0 ));
}

uint32_t Octree::getOctreeSize() const
{
    return _octreeSize;
}
