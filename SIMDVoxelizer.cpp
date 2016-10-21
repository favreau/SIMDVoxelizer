/* Copyright (c) 2015-2016, EPFL/Blue Brain Project
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

#include "SIMDVoxelizer_ispc.h"
#include <string.h>
#include <fstream>
#include <vector>
#include <limits>
#include <iostream>
#include <sstream>
#include <ctime>
#include <cstdlib>

using namespace ispc;

float cutoffDistance = 50.f;
float voxelSize = 2.f;
std::string inputFile;
std::string outputFile;

int main( int argc, char* argv[] )
{
    if( argc != 5 )
    {
        std::cerr << "usage: SIMDVoxelizer <voxel_size> <cutoff_distance> "
                  << "<input_file> <output_file>" << std::endl;
        exit(1);
    }

    voxelSize = atof(argv[1]);
    cutoffDistance = atof(argv[2]);
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

    size_t x = (maxx - minx) / voxelSize;
    size_t y = (maxy - miny) / voxelSize;
    size_t z = (maxz - minz) / voxelSize;
    size_t volume_size = x * y * z;

    std::stringstream filename;
    filename << outputFile;
    filename << "." << x << "_" << y << "_" << z;

    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "SIMDVoxelizer" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Voxel size        : " << voxelSize << std::endl;
    std::cout << "Cutoff distance   : " << cutoffDistance << std::endl;
    std::cout << "Input file        : " << inputFile << std::endl;
    std::cout << "Output file       : " << filename.str() << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    std::cout << "Volume dimensions ["
              << x << " " << y << " " << z << "] "
              << x * y * z / 1024 / 1024 << " Mb"
              << std::endl;

    float* volume = new float[ volume_size ];
    memset(volume, 0, sizeof(float) * volume_size);

    SIMDVoxelizer_ispc(
        voxelSize, x, y, z, minx, miny, minz, cutoffDistance, events.data(), events.size(), volume );

    float minValue = std::numeric_limits<float>::max();
    float maxValue = -std::numeric_limits<float>::max();
    for( int i = 0; i < volume_size; ++i )
    {
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
    std::cout << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    std::ofstream volumeFile(filename.str().c_str(), std::ios::out | std::ios::binary);
    volumeFile.write( (char*)&volumeAsChar[0], sizeof(char) * volume_size );
    volumeFile.close();
    delete [] volume;
    delete [] volumeAsChar;

    return 0;
}
