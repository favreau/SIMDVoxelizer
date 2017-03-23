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

#include <simdvoxelizer/SIMDSparseVoxelizer_ispc.h>
#include <simdvoxelizer/SIMDVoxelizer_ispc.h>
#include <simdvoxelizer/Octree.h>
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

int span = 32;
float voxelSize = 2.f;
std::string inputFile;
std::string outputFile;

typedef std::map< uint64_t, OctreeNode > OctreeLevelMap;

inline uint32_t pow2roundup( uint32_t x )
{
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
    const std::vector< OctreeNode* > children = node->getChildren();

    if(( children.empty( )) || ( level == 0 ))
    {
        flatOctreeData[ offsetPerLevel[ level ] * 4u ] = node->getCenter().x;
        flatOctreeData[ offsetPerLevel[ level ] * 4u + 1 ] = node->getCenter().y;
        flatOctreeData[ offsetPerLevel[ level ] * 4u + 2 ] = node->getCenter().z;
        flatOctreeData[ offsetPerLevel[ level ] * 4u + 3 ] = node->getMaxValue();

        offsetPerLevel[ level ] += 1u;
        return;
    }
    flatOctreeData[ offsetPerLevel[ level ] * 4u ] = node->getCenter().x;
    flatOctreeData[ offsetPerLevel[ level ] * 4u + 1 ] = node->getCenter().y;
    flatOctreeData[ offsetPerLevel[ level ] * 4u + 2 ] = node->getCenter().z;
    flatOctreeData[ offsetPerLevel[ level ] * 4u + 3 ] = node->getMaxValue();

    flatOctreeIndex[ offsetPerLevel[ level ] * 2u ] =
        offsetPerLevel[ level - 1 ];
    flatOctreeIndex[ offsetPerLevel[ level ] * 2u + 1 ] =
        offsetPerLevel[ level - 1 ] + children.size() - 1u;
    offsetPerLevel[ level ] += 1u;

    for( const OctreeNode* child : children )
        flattenChildren( child, offsetPerLevel, flatOctreeIndex, flatOctreeData,
                         level - 1u );
}

int main( /*int argc, char* argv[]*/ )
{
    /*if( argc != 5 )
    {
        std::cerr << "usage: SIMDVoxelizer <voxel_size> <span> "
                  << "<input_file> <output_file>" << std::endl;
        exit(1);
    }

    voxelSize = atof(argv[1]);
    span = atoi(argv[2]);
    inputFile = argv[3];
    outputFile = argv[4];
    std::vector<float> events;
    std::ifstream file(inputFile.c_str(), std::ios::in | std::ios::binary);
    if( file.is_open( ))
    {
        while( !file.eof( ))
        {
            float v;
            file.read( (char *)&v, sizeof(float));
            events.push_back(v);
        }
        file.close();
    }

    // Bounding box
    glm::vec3 min;
    min.x = std::numeric_limits<float>::max();
    min.y = std::numeric_limits<float>::max();
    min.z = std::numeric_limits<float>::max();

    glm::vec3 max;
    max.x = -std::numeric_limits<float>::max();
    max.y = -std::numeric_limits<float>::max();
    max.z = -std::numeric_limits<float>::max();

    std::cout << "Nb of events      :" <<  events.size() / 5 << std::endl;
    for( size_t i = 0; i < events.size(); i += 5 )
    {
        min.x = std::min( min.x, events[i] );
        min.y = std::min( min.y, events[i+1] );
        min.z = std::min( min.z, events[i+2] );

        max.x = std::max( max.x, events[i] );
        max.y = std::max( max.y, events[i+1] );
        max.z = std::max( max.z, events[i+2] );
    }

    //erase the soma position ( padded with zeros )-> posX posY posZ 0.0f 0.0f
    events.erase( events.begin(), events.begin() + 5 );

    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "SIMDVoxelizer" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Voxel size        : " << voxelSize << std::endl;
    std::cout << "Span              : " << span << std::endl;
    std::cout << "Input file        : " << inputFile << std::endl;
    std::cout << "Output file       : " << outputFile << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    Octree morphoOctree( events, voxelSize, min, max );

    uint64_t volumeSize = morphoOctree.getVolumeSize();
    glm::uvec3 volumeDim = morphoOctree.getVolumeDim();

    std::vector< float > volume( volumeSize, 0 );

    std::cout << "Volume dims: [" << volumeDim.x << ", " << volumeDim.y
              << ", " << volumeDim.z << "] "
              << volumeSize << " bytes" << std::endl;

    uint32_t zLenght = 32;
    for( uint32_t zOffset = 0; zOffset < volumeDim.z; zOffset += zLenght )
    {
        const size_t progress = float( zOffset ) / float( volumeDim.z ) * 100.f;
        std::cout << progress << "%\r";
        std::cout.flush();
        ispc::SIMDSparseVoxelizer_ispc(
            zOffset, zOffset + zLenght, span, voxelSize,
            morphoOctree.getOctreeSize(),
            volumeDim.x, volumeDim.y, volumeDim.z,
            morphoOctree.getFlatIndexes(),
            morphoOctree.getFlatData(), volume.data( ));
    }
    std::cout << std::endl;

    // Determine value range in volume
    float minValue = std::numeric_limits<float>::max();
    float maxValue = -std::numeric_limits<float>::max();
    for( size_t i = 0; i < volumeSize; ++i )
        if( volume[i] != 0.f )
        {
            minValue = std::min( minValue, volume[i] );
            maxValue = std::max( maxValue, volume[i] );
        }
    // Normalize volume
    float a = 255.f / ( maxValue - minValue );
    std::cout << "Normalization [" << minValue << " - "
              << maxValue << "] " << a << std::endl;

    std::vector< char > volumeAsChar( volumeSize );
    for( size_t i = 0; i < volumeSize; ++i )
        volumeAsChar[i] = static_cast< char >((volume[i] - minValue) * a);

    std::cout << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    std::ofstream volumeFile(
        outputFile.c_str(), std::ios::out | std::ios::binary);
    volumeFile.write( (char*)volumeAsChar.data(), sizeof(char) * volumeSize );
    volumeFile.close();*/

    return 0;
}
