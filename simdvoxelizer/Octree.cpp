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
    : _volumeDim( glm::uvec3( 0u, 0u, 0u ))
    , _volumeSize( 0u )
{
    glm::uvec3 octreeSize( _pow2roundup( std::ceil((max.x - min.x) / leafSize )),
                           _pow2roundup( std::ceil((max.y - min.y) / leafSize )),
                           _pow2roundup( std::ceil((max.z - min.z) / leafSize )));

    // This octree is always isotropic
    _octreeSize = std::max( std::max( octreeSize.x,
                                      octreeSize.y),
                                      octreeSize.z);

    _octreeDepth = std::log2( _octreeSize ) + 1u;
    _octree.resize( _octreeDepth );

    std::cout<<"DEPTH: " << _octreeDepth << " " << _octree.size() << std::endl;

    for( const Event& event : events )
    {
        const uint64_t xpos = std::floor(( event.position.x - min.x ) / leafSize );
        const uint64_t ypos = std::floor(( event.position.y - min.y ) / leafSize );
        const uint64_t zpos = std::floor(( event.position.z - min.z ) / leafSize );

        const uint64_t indexX = xpos;
        const uint64_t indexY = ypos * (uint64_t)_octreeSize;
        const uint64_t indexZ = zpos * (uint64_t)_octreeSize * (uint64_t)_octreeSize;

        auto it = _octree[ 0 ].find( indexX + indexY + indexZ );
        if( it == _octree[ 0 ].end( ))
        {
            OctreeNode* child = nullptr;
            for( uint32_t level = 0; level < _octreeDepth; ++level )
            {
                bool newNode = false;
                const uint64_t divisor = std::pow( 2, level );
                const glm::vec3 center( leafSize * divisor * ( xpos / divisor + 0.5f ) + min.x,
                                        leafSize * divisor * ( ypos / divisor + 0.5f ) + min.y,
                                        leafSize * divisor * ( zpos / divisor + 0.5f ) + min.z );


                const uint64_t nBlock = _octreeSize / divisor;
                const uint64_t index = std::floor( xpos / divisor ) +
                                       nBlock * std::floor( ypos / divisor ) +
                                       nBlock * nBlock * std::floor( zpos / divisor );

                const float size = leafSize * std::pow( 2, level );

                auto itr = _octree[ level ].find( index );
                if( itr == _octree[ level ].end( ))
                {
                    _octree[ level ].insert(
                                OctreeLevelMap::value_type(
                                    index, OctreeNode( center, size )));
                    newNode = true;
                }

                if( level == 0 )
                {
                    _octree[ level ].at( index ).addEvent( event );
                }
                else
                    _octree[ level ].at( index ).addMaxValue( event.value );

                if(( level > 0 ) && ( child != nullptr ))
                    _octree[level].at( index ).setChild( child );

                if( newNode )
                    child = &(_octree[ level ].at( index ));
                else
                    child = nullptr;
            }
        }
        else
        {
            for( uint32_t level = 0; level < _octreeDepth; ++level )
            {
                const uint64_t divisor = std::pow( 2, level );
                const uint64_t nBlock = _octreeSize / divisor;
                const uint64_t index = std::floor( xpos / divisor ) +
                                       nBlock * std::floor( ypos / divisor ) +
                                       nBlock * nBlock * std::floor( zpos / divisor );
                if( level == 0 )
                {
                    _octree[ level ].at( index ).addEvent( event );
                }
                else
                    _octree[level].at(index).addMaxValue( event.value );
            }
        }
    }
    for( uint32_t i = 0; i < _octree.size(); ++i )
        std::cout<< "Number of leaves ["<<i<<"]: " << _octree[i].size() << std::endl;

    _volumeDim = glm::uvec3( std::ceil((max.x - min.x) / leafSize ),
                             std::ceil((max.y - min.y) / leafSize ),
                             std::ceil((max.z - min.z) / leafSize ));
    _volumeSize = (uint64_t)_volumeDim.x * (uint64_t)_volumeDim.y *
                  (uint64_t)_volumeDim.z;
}

Octree::~Octree()
{}

const OctreeNode* Octree::getRoot() const
{
    return &(_octree[ _octreeDepth - 1 ].at( 0 ));
}

uint32_t Octree::getOctreeSize() const
{
    return _octreeSize;
}
