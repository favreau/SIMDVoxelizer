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
#include <map>
#include <string.h>

#include "Octree.h"

typedef std::map< uint64_t, OctreeNode > OctreeLevelMap;

Octree::Octree( const std::vector<float>& events, float voxelSize )
    : _volumeDim( glm::uvec3( 0u, 0u, 0u ))
    , _volumeSize( 0u )
    , _offsetPerLevel( nullptr )
    , _flatIndexes( nullptr )
    , _flatData( nullptr )
{
    // Bounding box
    float minx = std::numeric_limits<float>::max();
    float miny = std::numeric_limits<float>::max();
    float minz = std::numeric_limits<float>::max();

    float maxx = -std::numeric_limits<float>::max();
    float maxy = -std::numeric_limits<float>::max();
    float maxz = -std::numeric_limits<float>::max();

    std::cout << "Nb of events      :" <<  events.size() / 5 << std::endl;
    for( size_t i = 0; i < events.size(); i += 5 )
    {
        float value = events[ i + 4 ];
        if( value != 0.f )
        {
            minx = std::min( minx, events[i] );
            miny = std::min( miny, events[i+1] );
            minz = std::min( minz, events[i+2] );

            maxx = std::max( maxx, events[i] );
            maxy = std::max( maxy, events[i+1] );
            maxz = std::max( maxz, events[i+2] );
        }
    }
    std::cout << "Bounding box      ["
              << minx << "," << miny << "," << minz << "] ["
              << maxx << "," << maxy << "," << maxz << "]" << std::endl;

    // **************** Octree creations *******************
    // *****************************************************
    glm::uvec3 octreeSize( _pow2roundup( std::ceil((maxx - minx) / voxelSize )),
                           _pow2roundup( std::ceil((maxy - miny) / voxelSize )),
                           _pow2roundup( std::ceil((maxz - minz) / voxelSize )));

    // This octree is always cubic
    _octreeSize = std::max( std::max( octreeSize.x,
                                      octreeSize.y),
                                      octreeSize.z);

    uint32_t octreeDepth = std::log2( _octreeSize ) + 1u;
    std::vector< OctreeLevelMap > octree( octreeDepth );

    std::cout<<"DEPTH: " << octreeDepth << " " << octree.size() << std::endl;

    for( uint32_t i = 0; i < events.size(); i += 5 )
    {
        const uint64_t xpos = std::floor(( events[i] - minx ) / voxelSize );
        const uint64_t ypos = std::floor(( events[i + 1] - miny ) / voxelSize );
        const uint64_t zpos = std::floor(( events[i + 2] - minz ) / voxelSize );
        const float value = events[ i + 4 ];

        const uint64_t indexX = xpos;
        const uint64_t indexY = ypos * (uint64_t)_octreeSize;
        const uint64_t indexZ = zpos * (uint64_t)_octreeSize * (uint64_t)_octreeSize;

        auto it = octree[ 0 ].find( indexX + indexY + indexZ );
        if( it == octree[ 0 ].end( ))
        {
            OctreeNode* child = nullptr;
            for( uint32_t level = 0; level < octreeDepth; ++level )
            {
                bool newNode = false;
                const uint64_t divisor = std::pow( 2, level );
                const glm::vec3 center( divisor * ( xpos / divisor + 0.5f ),
                                        divisor * ( ypos / divisor + 0.5f ),
                                        divisor * ( zpos / divisor + 0.5f ));

                const uint64_t nBlock = _octreeSize / divisor;
                const uint64_t index = std::floor( xpos / divisor ) +
                                       nBlock * std::floor( ypos / divisor ) +
                                       nBlock * nBlock * std::floor( zpos / divisor );

                const float size = voxelSize * ( level + 1u );

                if( octree[ level ].find( index ) == octree[ level ].end( ))
                {
                    octree[ level ].insert(
                                OctreeLevelMap::value_type(
                                    index, OctreeNode( center, size )));
                    newNode = true;
                }

                octree[ level ].at( index ).addValue( value );

                if(( level > 0 ) && ( child != nullptr ))
                    octree[level].at( index ).setChild( child );

                if( newNode )
                    child = &(octree[ level ].at( index ));
                else
                    child = nullptr;
            }
        }
        else
        {
            for( uint32_t level = 0; level < octreeDepth; ++level )
            {
                const uint64_t divisor = std::pow( 2, level );
                const uint64_t nBlock = _octreeSize / divisor;
                const uint64_t index = std::floor( xpos / divisor ) +
                                       nBlock * std::floor( ypos / divisor ) +
                                       nBlock * nBlock * std::floor( zpos / divisor );
                octree[level].at(index).addValue( value );
            }
        }
    }
    for( uint32_t i = 0; i < octree.size(); ++i )
        std::cout<< "Number of leaves ["<<i<<"]: " << octree[i].size() << std::endl;

    // **************** Octree flattening *******************
    // ******************************************************

    _offsetPerLevel = new uint32_t[ octreeDepth ];
    _offsetPerLevel[ octreeDepth - 1u ] = 0;
    uint32_t previousOffset = 0u;
    for( uint32_t i = octreeDepth - 1u; i > 0u; --i )
    {
        _offsetPerLevel[ i - 1u ] = previousOffset + octree[ i ].size();
        previousOffset = _offsetPerLevel[ i - 1u ];
    }

    uint32_t totalNodeNumber = 0;

    for( uint32_t i = 0; i < octree.size(); ++i )
        totalNodeNumber += octree[i].size();

    // need to be initialized with zeros
    _flatIndexes = new uint32_t[ totalNodeNumber * 2u ];
    memset( _flatIndexes, 0, totalNodeNumber * 8u );

    _flatData = new float[ totalNodeNumber * 4 ];

    //The root node
    _flattenChildren( &(octree[ octreeDepth - 1u ].at( 0 )), octreeDepth - 1u );

    // **************** Octree flattening end *****************
    // ********************************************************

    _volumeDim = glm::uvec3( std::ceil((maxx - minx) / voxelSize ),
                             std::ceil((maxy - miny) / voxelSize ),
                             std::ceil((maxz - minz) / voxelSize ));
    _volumeSize = (uint64_t)_volumeDim.x * (uint64_t)_volumeDim.y *
                  (uint64_t)_volumeDim.z;
}

Octree::~Octree()
{
    delete[] _flatData;
    delete[] _flatIndexes;
    delete[] _offsetPerLevel;
}

void Octree::_flattenChildren( const OctreeNode* node, uint32_t level )
{
    const std::vector< OctreeNode* > children = node->getChildren();

    if(( children.empty( )) || ( level == 0 ))
    {
        _flatData[ _offsetPerLevel[ level ] * 4u ] = node->getCenter().x;
        _flatData[ _offsetPerLevel[ level ] * 4u + 1 ] = node->getCenter().y;
        _flatData[ _offsetPerLevel[ level ] * 4u + 2 ] = node->getCenter().z;
        _flatData[ _offsetPerLevel[ level ] * 4u + 3 ] = node->getValue();

        _offsetPerLevel[ level ] += 1u;
        return;
    }
    _flatData[ _offsetPerLevel[ level ] * 4u ] = node->getCenter().x;
    _flatData[ _offsetPerLevel[ level ] * 4u + 1 ] = node->getCenter().y;
    _flatData[ _offsetPerLevel[ level ] * 4u + 2 ] = node->getCenter().z;
    _flatData[ _offsetPerLevel[ level ] * 4u + 3 ] = node->getValue();

    _flatIndexes[ _offsetPerLevel[ level ] * 2u ] = _offsetPerLevel[ level - 1 ];
    _flatIndexes[ _offsetPerLevel[ level ] * 2u + 1 ] = _offsetPerLevel[ level - 1 ] + children.size() - 1u;
    _offsetPerLevel[ level ] += 1u;

    for( const OctreeNode* child : children )
        _flattenChildren( child, level - 1u );
}

uint32_t Octree::getOctreeSize() const
{
    return _octreeSize;
}

uint32_t* Octree::getFlatIndexes() const
{
    return _flatIndexes;
}

float* Octree::getFlatData() const
{
    return _flatData;
}

const glm::uvec3& Octree::getVolumeDim() const
{
    return _volumeDim;
}

uint64_t Octree::getVolumeSize() const
{
    return _volumeSize;
}
