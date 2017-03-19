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

#include <boost/program_options.hpp>

#include <simdvoxelizer/Mesher.h>

int main( int argc, char* argv[] )
{
    float alpha = 3.0f;
    float angularBound = -1.0f;
    float radialBound = -1.0f;
    float distanceBound = -1.0f;
    float leafSize = -1.0f;
    float minRadius = 0.1f;
    float decimationRatio = 1.0f;
    float epsilon = 0.0000033f;
    uint32_t targetFaceNumber = 0u;

    std::string outputFile;
    std::string inputFile;

    namespace po = boost::program_options;
    po::options_description desc("");

    desc.add_options()
      ("help,h", "Print this help messages.")
      ("input,i", po::value< std::string >( &inputFile )->required(),
       "Input file containing the morphology.")
      ("output,o", po::value< std::string >( &outputFile )->required(),
       "Output file that will contain the mesh.\n")
      ("simplification-ratio,s", po::value< float >( &decimationRatio )->default_value( 1.0f ),
       "The decimation ratio by which the mesh must be decimated.\n")
      ("target-faces,t", po::value< uint32_t >( &targetFaceNumber ),
       "The target number of faces for the decimated mesh. This parameter overwrite the"
       " --simplification-ratio parameter.\n")
      ("min-radius", po::value< float >( &minRadius )->default_value( 0.1f ),
       "All radii in the morphology smaller than the specified value will be "
       "increased to that value.")
      ("angular-bound,a", po::value< float >( &angularBound ),
       "A lower bound in degrees for the angles of mesh facets. If not specified "
       "this parameter is automatically computed.")
      ("radial-bound,r", po::value< float >( &radialBound ),
       "An upper bound on the radii of surface Delaunay balls. A surface "
       "Delaunay ball is a ball circumscribing a mesh facet and centered on the "
       "surface. If not specified this parameter is automatically computed.")
      ("distance-bound,d", po::value< float >( &distanceBound ),
       "An upper bound for the distance between the circumcenter of a mesh facet "
       "and the center of a surface Delaunay ball of this facet. If not specified "
       "this parameter is automatically computed.\n")
      ("epsilon,e", po::value< float >( &epsilon )->default_value( 0.0000033f ),
       "Epsilon multiplyed by the radius of the bounding sphere is the precision at"
       "which the dichotomy will stop.\n")
      ("leaf-size,l", po::value< float >( &leafSize ),
       "The size of the octree's leaves in morphology space. If not specified "
       "this parameter is automatically computed.")
      ("alpha", po::value< float >( &alpha )->default_value( 3.0f ),
       "Octree's nodes further than closest radius * alpha are discarded");

    po::variables_map vm;

    try
    {
        po::store( po::parse_command_line( argc, argv, desc ), vm );

        if( vm.count( "help" ))
        {
            std::cout << desc << std::endl;
            return 0;
        }
        po::notify( vm );
    }
    catch( boost::program_options::required_option& e )
    {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }

    if( minRadius <= 0 )
    {
        std::cerr << " ERROR: min-radius is: " << minRadius
                  <<". The algorithm will not converge." << std::endl;
        return 1;
    }

    Mesher mesher( inputFile, alpha, angularBound, radialBound, distanceBound,
                   leafSize, minRadius, decimationRatio, epsilon, targetFaceNumber );
    mesher.mesh( outputFile );

    return 0;
}
