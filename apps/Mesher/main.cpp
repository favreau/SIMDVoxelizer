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

#if 0 // WORK IN PROGRESS
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <fstream>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;
typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_3;

using namespace ispc;

int span = 32;
float voxelSize = 2.f;
float threshold = 15.f;
std::string inputFile;
std::string outputFile;
float triangleSize = 5;
int generateVolume = 1;

int main( int argc, char* argv[] )
{
    if( argc != 8 )
    {
        std::cerr << "usage: SIMDVoxelizer <voxel_size> <span> <input_file> <output_file> <threshold> <triangle_size> <generate_volume>" << std::endl;
        exit(1);
    }

    voxelSize = atof(argv[1]);
    span = atoi(argv[2]);
    inputFile = argv[3];
    outputFile = argv[4];
    threshold = atof(argv[5]);
    triangleSize = atof(argv[6]);
    generateVolume = atoi(argv[7]);

    glm::uvec3 volumeDim( 534, 829, 253);
    if( generateVolume == 1 )
    {
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

        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "SIMDVoxelizer" << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "Voxel size        : " << voxelSize << std::endl;
        std::cout << "Span              : " << span << std::endl;
        std::cout << "Input file        : " << inputFile << std::endl;
        std::cout << "Output file       : " << outputFile << std::endl;
        std::cout << "--------------------------------------------" << std::endl;

        Octree morphoOctree( events, voxelSize );
        uint64_t volumeSize = morphoOctree.getVolumeSize();
        volumeDim = morphoOctree.getVolumeDim();

        float* volume = new float[ volumeSize ];

        std::cout<< "Volume dims: " << volumeDim.x << " " << volumeDim.y << " " << volumeDim.z << " " << volumeSize << std::endl;

        uint32_t zLenght = 32;
        uint32_t nPasses = std::ceil( volumeDim.z / zLenght );

        for( uint32_t zOffset = 0; zOffset < volumeDim.z; zOffset += zLenght )
        {
            std::cout << "z: " << zOffset << std::endl;
            SIMDSparseVoxelizer_ispc( zOffset, zOffset + zLenght, span, voxelSize,
                                      morphoOctree.getOctreeSize(), volumeDim.x,
                                      volumeDim.y, volumeDim.z,
                                      morphoOctree.getFlatIndexes(),
                                      morphoOctree.getFlatData(), volume );
        }

        float minValue = std::numeric_limits<float>::max();
        float maxValue = -std::numeric_limits<float>::max();
        for( int i = 0; i < volumeSize; ++i )
        {
            if( volume[i] != 0.f )
            {
                minValue = std::min( minValue, volume[i] );
                maxValue = std::max( maxValue, volume[i] );
            }
        }

        //float a = 255.f / ( maxValue - minValue );
        float a = 9.f / ( maxValue - minValue );
        std::cout << "Normalization [" << minValue << " - " << maxValue << "] " << a << std::endl;

        char* volumeAsChar = new char[ volumeSize ];
        for( int i = 0; i < volumeSize; ++i )
        {
            //float normalizedValue = (volume[i] - minValue) * a;
            float normalizedValue = 255.0f * std::log((volume[i] - minValue) * a + 1.0);
            volumeAsChar[i] = (uint8_t)normalizedValue;
        }
        std::cout << std::endl;
        std::cout << "--------------------------------------------" << std::endl;

        std::ofstream volumeFile(
                    outputFile.c_str(), std::ios::out | std::ios::binary);
        std::stringstream header;
        header << "#INRIMAGE-4#{\n";
        header << "XDIM=" << volumeDim.x << "\n";
        header << "YDIM=" << volumeDim.y << "\n";
        header << "ZDIM=" << volumeDim.z << "\n";
        header << "VDIM=1\n";
        header << "TYPE=unsigned fixed\n";
        header << "SCALE=1**0\n";
        header << "PIXSIZE=8 bits\n";
        header << "CPU=pc\n";
        header << "VX=1\n";
        header << "VY=1\n";
        header << "VZ=1\n";
        const size_t len = 256 - 4 - header.str().length();
        for( size_t i =0 ; i < len; ++i )
            header << "\n";
        header << "##}\n";

        volumeFile.write( (char*)header.str().c_str(), sizeof(char) * header.str().length());
        volumeFile.write( (char*)&volumeAsChar[0], sizeof(char) * volumeSize );
        volumeFile.close();
        delete [] volume;
        delete [] volumeAsChar;
    }

    std::cout << "Meshing..." << std::endl;
    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
    // the 'function' is a 3D gray level image
    Gray_level_image image( outputFile.c_str(), threshold );

    // Carefully choosen bounding sphere: the center must be inside the
    // surface defined by 'image' and the radius must be high enough so that
    // the sphere actually bounds the whole image.
    GT::Point_3 bounding_sphere_center(
                0.498567f * volumeDim.x,
                0.199998f * volumeDim.y,
                0.609048f * volumeDim.z );
    const float radius =
        std::max( volumeDim.x, std::max( volumeDim.y, volumeDim.z ));
    std::cout << "Radius: " << radius << std::endl;
    GT::FT bounding_sphere_squared_radius = radius * radius;
    GT::Sphere_3 bounding_sphere(bounding_sphere_center,
                                 bounding_sphere_squared_radius);

    // definition of the surface, with 10^-5 as relative precision
    Surface_3 surface( image, bounding_sphere, 1e-7 );

    // defining meshing criteria
    std::cout << "Triangle Size: " << triangleSize << std::endl;
    CGAL::Surface_mesh_default_criteria_3<Tr> criteria( 30, triangleSize, triangleSize );

    // meshing surface, with the "manifold without boundary" algorithm
    //CGAL::make_surface_mesh( c2t3, surface, criteria, CGAL::Manifold_tag( ));
    CGAL::make_surface_mesh( c2t3, surface, criteria, CGAL::Non_manifold_tag());
    std::ofstream out( "out.off" );
    CGAL::output_surface_facets_to_off( out, c2t3 );
    std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
    return 0;
}
#else
int main( int , char* argv[] )
{
    std::cout << argv << std::endl;
    return 1;
}

#endif
