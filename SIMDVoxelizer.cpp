/* Copyright (c) 2015-2017, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Cyrille Favreau <cyrille.favreau@epfl.ch>
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

#include "SIMDSparseVoxelizer_ispc.h"
#include "Octree.h"
#include <string.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <map>
#include <limits>
#include <iostream>
#include <sstream>
#include <ctime>
#include <cstdlib>

using namespace ispc;

int span = 32;
float voxelSize = 2.f;
std::string inputFile;
std::string outputFile;

typedef std::map< uint64_t, OctreeNode > OctreeLevelMap;

inline uint32_t pow2roundup( uint32_t x )
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

void flattenChildren( const OctreeNode* node, uint32_t* offsetPerLevel,
                      uint32_t* flatOctreeIndex, float* flatOctreeData,
                      uint32_t level )
{
    /*std::cout<< "In level: " << level << " - "
             << offsetPerLevel[ 0 ] << " "
             << offsetPerLevel[ 1 ] << " "
             << offsetPerLevel[ 2 ] << " "
             << offsetPerLevel[ 3 ] << " "
             << offsetPerLevel[ 4 ] << " "
             << offsetPerLevel[ 5 ] << " "
             << offsetPerLevel[ 6 ] << " "
             << offsetPerLevel[ 7 ] << " "
             << offsetPerLevel[ 8 ] << " "
             << offsetPerLevel[ 9 ] << " "
             << offsetPerLevel[ 10 ]
             << std::endl;*/
    const std::vector< OctreeNode* > children = node->getChildren();

    if(( children.empty( )) || ( level == 0 ))
    {
        flatOctreeData[ offsetPerLevel[ level ] * 4u ] = node->getCenter().x;
        flatOctreeData[ offsetPerLevel[ level ] * 4u + 1 ] = node->getCenter().y;
        flatOctreeData[ offsetPerLevel[ level ] * 4u + 2 ] = node->getCenter().z;
        flatOctreeData[ offsetPerLevel[ level ] * 4u + 3 ] = node->getValue();

        /*if( node->getCenter().x > 256 )
        {
            std::cout<< "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
        }
        if( node->getCenter().y > 256 )
        {
            std::cout<< "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYY" << std::endl;
        }*/
        if( node->getCenter().z > 256 )
        {
            std::cout<< "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ" << std::endl;
        }

        offsetPerLevel[ level ] += 1u;
        return;
    }
    flatOctreeData[ offsetPerLevel[ level ] * 4u ] = node->getCenter().x;
    flatOctreeData[ offsetPerLevel[ level ] * 4u + 1 ] = node->getCenter().y;
    flatOctreeData[ offsetPerLevel[ level ] * 4u + 2 ] = node->getCenter().z;
    flatOctreeData[ offsetPerLevel[ level ] * 4u + 3 ] = node->getValue();

    flatOctreeIndex[ offsetPerLevel[ level ] * 2u ] = offsetPerLevel[ level - 1 ];
    flatOctreeIndex[ offsetPerLevel[ level ] * 2u + 1 ] = offsetPerLevel[ level - 1 ] + children.size() - 1u;
    offsetPerLevel[ level ] += 1u;

    for( const OctreeNode* child : children )
        flattenChildren( child, offsetPerLevel, flatOctreeIndex, flatOctreeData,
                         level - 1u );

    //std::cout<<"Ok" << std::endl;
}

int main( int argc, char* argv[] )
{
    if( argc != 5 )
    {
        std::cerr << "usage: SIMDVoxelizer <voxel_size> <cutoff_distance> "
                  << "<input_file> <output_file>" << std::endl;
        exit(1);
    }

    voxelSize = atof(argv[1]);
    span = atoi(argv[2]);
    inputFile = argv[3];
    outputFile = argv[4];
    std::vector<float> events;
    std::ifstream file(inputFile.c_str(), std::ios::in | std::ios::binary);
    if( file.is_open() )
    {
        while( !file.eof() )
        {
            float v;
            file.read( (char *)&v, sizeof(float));
            events.push_back(v);
        }
        file.close();
    }

    // Bounding box
    float minx = std::numeric_limits<float>::max();
    float miny = std::numeric_limits<float>::max();
    float minz = std::numeric_limits<float>::max();
    float minv = std::numeric_limits<float>::max();

    float maxx = -std::numeric_limits<float>::max();
    float maxy = -std::numeric_limits<float>::max();
    float maxz = -std::numeric_limits<float>::max();
    float maxv = -std::numeric_limits<float>::max();

    std::cout << "Nb of events      :" <<  events.size() / 5 << std::endl;
    for( size_t i = 0; i < events.size(); i += 5 )
    {
        float value = events[ i + 4 ];
        if( value != 0.f )
        {
            minx = std::min( minx, events[i] );
            miny = std::min( miny, events[i+1] );
            minz = std::min( minz, events[i+2] );
            minv = std::min( minv, events[i+3] );

            maxx = std::max( maxx, events[i] );
            maxy = std::max( maxy, events[i+1] );
            maxz = std::max( maxz, events[i+2] );
            maxv = std::max( maxv, events[i+3] );
        }
    }
    std::cout << "Bounding box      ["
              << minx << "," << miny << "," << minz << "] ["
              << maxx << "," << maxy << "," << maxz << "]" << std::endl;

    // !! DANGER !! THIS IS WRONG !!
    //size_t x = (maxx - minx) / voxelSize;
    //size_t y = (maxy - miny) / voxelSize;
    //size_t z = (maxz - minz) / voxelSize;

    // **************** Octree creations *******************
    // *****************************************************
    glm::uvec3 octreeSize( pow2roundup( std::ceil((maxx - minx) / voxelSize )),
                           pow2roundup( std::ceil((maxy - miny) / voxelSize )),
                           pow2roundup( std::ceil((maxz - minz) / voxelSize )));

    uint32_t maxOctreeSize = std::max( std::max( octreeSize.x,
                                                 octreeSize.y),
                                                 octreeSize.z);

    glm::uvec3 centerEvents = glm::uvec3( std::ceil((maxx - minx) / voxelSize ) / 2,
                                          std::ceil((maxy - miny) / voxelSize ) / 2,
                                          std::ceil((maxz - minz) / voxelSize ) / 2 );

    glm::uvec3 centerVoxels = glm::uvec3( maxOctreeSize / 2.0 );
    glm::uvec3 offsetToCenter = centerVoxels - centerEvents;

    std::cout<<" Offset to center: " << offsetToCenter.x << " " << offsetToCenter.y << " " << offsetToCenter.z << std::endl;

    size_t volume_size = maxOctreeSize * maxOctreeSize * maxOctreeSize;

    uint32_t octreeDepth = std::log2( maxOctreeSize ) + 1u;
    std::vector< OctreeLevelMap > octree( octreeDepth );

    std::cout<<"DEPTH: " << octreeDepth << " " << octree.size() << std::endl;

    for( uint32_t i = 0; i < events.size(); i += 5 )
    {
        const uint64_t xpos = std::floor(( events[i] - minx ) / voxelSize ) /*+ offsetToCenter.x*/;
        const uint64_t ypos = std::floor(( events[i + 1] - miny ) / voxelSize ) /*+ offsetToCenter.y*/;
        const uint64_t zpos = std::floor(( events[i + 2] - minz ) / voxelSize ) /*+ offsetToCenter.z*/;
        const float value = events[ i + 4 ];

        const uint64_t indexX = xpos;
        const uint64_t indexY = ypos * (uint64_t)maxOctreeSize;
        const uint64_t indexZ = zpos * (uint64_t)maxOctreeSize * (uint64_t)maxOctreeSize;

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

                /*if( level == 10 )
                {
                    std::cout<< "Center: " << center.x << " " << center.y << " " << center.z
                             << " pos: " << xpos << " " << ypos << " " << zpos << std::endl;
                }*/

                const uint64_t nBlock = maxOctreeSize / divisor;
                const uint64_t index = std::floor( xpos / divisor ) +
                                       nBlock * std::floor( ypos / divisor ) +
                                       nBlock * nBlock * std::floor( zpos / divisor );

                const float size = voxelSize * ( level + 1u );

                auto it = octree[ level ].find( index );
                if( it == octree[ level ].end( ))
                {
                    /*if( level == 2 )
                        std::cout << "ROOT: " << index << " " << divisor
                                  << " " << xpos << " " <<ypos << " " << zpos
                                  << " " << indexX << " " <<indexY << " " << indexZ
                                  << " center: " << center.x << " " << center.y << " " << center.z
                                  << std::endl;*/
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
                const uint64_t nBlock = maxOctreeSize / divisor;
                const uint64_t index = std::floor( xpos / divisor ) +
                                       nBlock * std::floor( ypos / divisor ) +
                                       nBlock * nBlock * std::floor( zpos / divisor );
                octree[level].at(index).addValue( value );
            }
        }
    }
    for( uint32_t i = 0; i < octree.size(); ++i )
        std::cout<< "Number of leaves ["<<i<<"]: " << octree[i].size() << std::endl;

    std::cout<<"MAX OCTREE SIZE: " << maxOctreeSize << std::endl;

    std::cout<<"x: " << ( maxx - minx ) / voxelSize
             <<" y: " << ( maxy - miny ) / voxelSize
             <<" z: " << ( maxz - minz ) / voxelSize << std::endl;
    std::cout<<"x: " << octreeSize.x
             <<" y: " << octreeSize.y
             <<" z: " << octreeSize.z << std::endl;

    // **************** Octree flattening *******************
    // ******************************************************

    uint32_t offsetPerLevel[ octreeDepth ];
    offsetPerLevel[ octreeDepth - 1u ] = 0;
    uint32_t previousOffset = 0u;
    for( uint32_t i = octreeDepth - 1u; i > 0u; --i )
    {
        offsetPerLevel[ i - 1u ] = previousOffset + octree[ i ].size();
        previousOffset = offsetPerLevel[ i - 1u ];
        std::cout<< "LEVEL " << i - 1 << " OFFSET: " << offsetPerLevel[i - 1] << std::endl;
    }

    uint32_t totalNodeNumber = 0;

    for( uint32_t i = 0; i < octree.size(); ++i )
        totalNodeNumber += octree[i].size();

    uint32_t flatOctreeIndex[ totalNodeNumber * 2u ];
    // need to be initialized with zeros
    memset( &flatOctreeIndex[0], 0, totalNodeNumber * 8u );

    float flatOctreeData[ totalNodeNumber * 4 ];

    //The root node
    flattenChildren( &(octree[ octreeDepth - 1u ].at( 0 )), offsetPerLevel, flatOctreeIndex, flatOctreeData, octreeDepth - 1u );

    //for( uint32_t i = 0; i < 100; ++i )
    //    std::cout<< i << " index: " << flatOctreeIndex[ i ] << std::endl;

    //for( uint32_t i = 0; i < 10000; ++i )
    //    std::cout<< i << " data: " << flatOctreeData[ i ] << std::endl;

    // **************** Octree flattening end *******************
    // **********************************************************

    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "SIMDVoxelizer" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Voxel size        : " << voxelSize << std::endl;
    std::cout << "Span              : " << span << std::endl;
    std::cout << "Input file        : " << inputFile << std::endl;
    std::cout << "Output file       : " << outputFile << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    /*std::cout << "Volume dimensions ["
              << x << " " << y << " " << z << "] "
              << x * y * z / 1024 / 1024 << " Mb"
              << std::endl;*/

    /*float* volume = new float[ volume_size ];

    memset(volume, 0, sizeof(float) * volume_size);

    SIMDVoxelizer_ispc(
        voxelSize, x, y, z, minx, miny, minz, cutoffDistance, events.data(), events.size(), volume );*/

    float* volume = new float[ volume_size ];
    //memset(volume, 0, sizeof(float) * volume_size);
    std::cout<<"TEST1" << std::endl;

    SIMDSparseVoxelizer_ispc( span, voxelSize, maxOctreeSize, flatOctreeIndex, flatOctreeData, volume );

    float minValue = std::numeric_limits<float>::max();
    float maxValue = -std::numeric_limits<float>::max();
    for( int i = 0; i < maxOctreeSize * maxOctreeSize * maxOctreeSize; ++i )
    {
        if( volume[i] > 1000000 )
            std::cout<< i <<": " << volume[i] << std::endl;

        if( volume[i] != 0.f )
        {
            minValue = std::min( minValue, volume[i] );
            maxValue = std::max( maxValue, volume[i] );
        }
    }

    float a = 255.f / ( maxValue - minValue );
    std::cout << "Normalization [" << minValue << " - " << maxValue << "] " << a << std::endl;

    char* volumeAsChar = new char[ volume_size ];
    for( int i = 0; i < volume_size; ++i )
    {
        float normalizedValue = (volume[i] - minValue) * a;
        volumeAsChar[i] = (uint8_t)normalizedValue;
    }
    std::cout<<"TEST4" << std::endl;
    std::cout << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    std::ofstream volumeFile(
        outputFile.c_str(), std::ios::out | std::ios::binary);
    volumeFile.write( (char*)&volumeAsChar[0], sizeof(char) * volume_size );
    volumeFile.close();
    std::cout<<"TEST5" << std::endl;
    delete [] volume;
    delete [] volumeAsChar;

    return 0;
}
