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
#include <fstream>
#include <algorithm>
#include <functional>

#include "Mesher.h"


namespace SMS = CGAL::Surface_mesh_simplification;

bool compareByValue( const Event& e1, const Event& e2 )
{
    return e1.value > e2.value;
}

Mesher::Mesher( const std::string& inputFile, const float alpha,
                const float angularBound, const float radialBound,
                const float distanceBound, const float leafSize,
                const float minRadiusForced, float decimationRatio )
    : _nCgalQuery( 0 )
    , _maxEventValue( -std::numeric_limits<float>::max( ))
    , _decimationRatio( decimationRatio )
    , _alpha( alpha )
    , _octree( nullptr )
{
    std::vector<float> events;
    std::ifstream file(inputFile.c_str(), std::ios::in | std::ios::binary);
    if( file.is_open( ))
    {
        while( true )
        {
            float v;
            file.read( (char *)&v, sizeof(float));
            if( file.eof( ))
                break;
            events.push_back( v );
        }
        file.close();
    }
    std::vector< Event > sortedEvents( events.size() / 5u );
    for( uint32_t i = 0; i < events.size(); i += 5 )
    {
        sortedEvents[ i / 5u ].position = glm::vec3( events[ i ],
                                                     events[ i + 1u ],
                                                     events[ i + 2u ] );
        if( minRadiusForced == -1.0f )
            sortedEvents[ i / 5u ].value = std::pow( events[ i + 3u ], 4 );
        else
        {
            if( events[ i + 3u ] < minRadiusForced )
                sortedEvents[ i / 5u ].value = std::pow( minRadiusForced, 4 );
            else
                sortedEvents[ i / 5u ].value = std::pow( events[ i + 3u ], 4 );
        }
    }
    std::sort( sortedEvents.begin(), sortedEvents.end(), compareByValue );
    _outLiners = _extractOutLiners( sortedEvents, 300 );

    // Bounding values
    glm::vec3 min;
    min.x = std::numeric_limits<float>::max();
    min.y = std::numeric_limits<float>::max();
    min.z = std::numeric_limits<float>::max();
    float minValue = std::numeric_limits<float>::max();

    glm::vec3 max;
    max.x = -std::numeric_limits<float>::max();
    max.y = -std::numeric_limits<float>::max();
    max.z = -std::numeric_limits<float>::max();

    for( const Event& event: sortedEvents )
    {
        min.x = std::min( min.x, event.position.x );
        min.y = std::min( min.y, event.position.y );
        min.z = std::min( min.z, event.position.z );
        minValue = std::min( minValue, event.value );

        max.x = std::max( max.x, event.position.x );
        max.y = std::max( max.y, event.position.y );
        max.z = std::max( max.z, event.position.z );
        _maxEventValue = std::max( _maxEventValue, event.value );
    }
    _boundingBoxDiameter = std::max( std::max( max.x - min.x,
                                               max.y - min.y ),
                                               max.z - min.z );
    std::cout << "Bounding diameter: " << _boundingBoxDiameter << std::endl;
    std::cout << "Bounding box      ["
              << min.x << "," << min.y << "," << min.z << "] ["
              << max.x << "," << max.y << "," << max.z << "]" << std::endl;

    std::cout << "Bounding values      [" << minValue << " ] ["
              << _maxEventValue << "]" << std::endl;

    _somaPosition = _outLiners[0].position;

    if( angularBound == -1.0f )
        _angularBound = 10.0f;
    else
        _angularBound = angularBound;

    if( radialBound == -1.0f )
        _radialBound = 1.0f;
    else
        _radialBound = radialBound;

    if( distanceBound == -1.0f )
        _distanceBound = 1.0f;
    else
        _distanceBound = distanceBound;

    if( leafSize == -1.0f )
        _octreeLeafSize = 4.0;
    else
        _octreeLeafSize = leafSize;

    //!!!!!!!! Need an optimal leaf size to be computed!!
    _octree.reset( new Octree( sortedEvents, _octreeLeafSize, min, max ));
}

void Mesher::mesh( const std::string& outputFile )
{
    float epsilon = 0.0000033;
    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

    GT::Point_3 sphereCenter( _somaPosition.x, _somaPosition.y,
                              _somaPosition.z );

    Function getFieldValue = [&]( const Point_3& p ){
       return _getFieldValue( p );
    };
    std::cout<< "Center of the sphere's value: "
             << _getFieldValue( sphereCenter ) << std::endl;
    Surface_3 surface( getFieldValue,
                       Sphere_3( sphereCenter,
                                 _boundingBoxDiameter * _boundingBoxDiameter ),
                                 epsilon);

    std::cout<< " Dicotomy stop value: " << _boundingBoxDiameter * epsilon << std::endl;

    CGAL::Surface_mesh_default_criteria_3<Tr> criteria( _angularBound,
                                                        _radialBound,
                                                        _distanceBound );
    // meshing surface
    std::cout<<"meshing..." << std::endl;
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
    std::cout<<"\nmeshing finished" << std::endl;

    std::ofstream output( outputFile + ".off" );
    CGAL::output_surface_facets_to_off( output, c2t3 );
    std::cout << "Final number of vertices: " << tr.number_of_vertices() << "\n";

    // ************* MESH SIMPLIFICATION **************

    if( _decimationRatio >= 1.0f )
        return;

    std::cout<<"decimating..." << std::endl;
    SurfaceMesh mesh;
    CGAL::output_surface_facets_to_polyhedron( c2t3, mesh );

    // This is a stop predicate (defines when the algorithm terminates).
    // In this example, the simplification stops when the number of undirected edges
    // left in the surface mesh drops below the specified number (1000)
    SMS::Count_ratio_stop_predicate< SurfaceMesh > stop( _decimationRatio );
    /*int r = SMS::edge_collapse
              ( mesh
               ,stop
               ,CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, mesh ))
                                 .halfedge_index_map  (get(CGAL::halfedge_external_index  , mesh ))
                                 .get_cost (SMS::Edge_length_cost <SurfaceMesh>())
                                 .get_placement(SMS::Midpoint_placement<SurfaceMesh>())
              );*/
    /*int r = SMS::edge_collapse
              ( mesh
               ,stop
               ,CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, mesh ))
                                 .halfedge_index_map  (get(CGAL::halfedge_external_index  , mesh ))
                                 .get_cost (SMS::LindstromTurk_cost<SurfaceMesh>())
                                 .get_placement(SMS::Midpoint_placement<SurfaceMesh>())
              );*/
    int r = SMS::edge_collapse
                  ( mesh
                   ,stop
                   ,CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, mesh ))
                                     .halfedge_index_map  (get(CGAL::halfedge_external_index  , mesh ))
                                     .get_cost (SMS::LindstromTurk_cost<SurfaceMesh>())
                                     .get_placement(SMS::LindstromTurk_placement<SurfaceMesh>())
                  );
    std::cout << "\nFinished...\n" << r << " edges removed.\n"
              << ( mesh.size_of_halfedges() / 2 ) << " final edges.\n" ;
    std::ofstream outputDecimated( outputFile + "_decimated.off" );
    outputDecimated << mesh;
}

FT Mesher::_getFieldValue( const Point_3 p )
{
    float threshold = 1.0f;

    glm::vec3 location( p.x(), p.y(), p.z( ));
    float value = 0;
    for( const Event& outLiner : _outLiners )
    {
        glm::vec3 delta = glm::vec3( outLiner.position.x, outLiner.position.y,
                                     outLiner.position.z ) - location;
        float distSquared = delta.x * delta.x +
                            delta.y * delta.y +
                            delta.z * delta.z;
        // Avoid division by 0
        if( distSquared < 0.01f )
            distSquared = 0.01f;

        value += outLiner.value / ( distSquared * distSquared );
    }
    value += _parseOctree( _octree->getRoot(), location );

    ++_nCgalQuery;
    if( _nCgalQuery % 50000 == 0 )
        std::cout <<"\rNumber of CGAL's oracle queries: "
                  << _nCgalQuery << std::flush;

    return threshold - value;
}

float Mesher::_parseOctree( const OctreeNode* node, glm::vec3 location )
{
    std::vector< OctreeNode* > children = node->getChildren();
    if( children.empty( ))
    {
        std::vector< Event > events = node->getData();
        float value = 0;
        for( const Event& event : events )
        {
            glm::vec3 delta = glm::vec3( event.position.x, event.position.y, event.position.z ) - location;
            float distSquared = delta.x * delta.x +
                                delta.y * delta.y +
                                delta.z * delta.z;
            // Avoid division by 0
            if( distSquared < 0.01f )
                distSquared = 0.01f;

            value += event.value / ( distSquared * distSquared );
        }
        return value;
    }

    float value = 0;
    for( OctreeNode* childNode : children )
    {
        glm::vec3 delta = childNode->getCenter() - location;
        float childDist = sqrt( delta.x * delta.x + delta.y * delta.y +
                                delta.z * delta.z ) - childNode->getHalfDiagonal();
        if( childDist < 0 )
            childDist = 0;

        if( childDist <= _alpha * _maxEventValue )
            value += _parseOctree( childNode, location );
    }
    return value;
}

std::vector< Event > Mesher::_extractOutLiners( std::vector< Event >& sortedEvents,
                                                uint32_t nOutLiners )
{
    std::vector< Event > outLiners;
    for( uint32_t i = 0; i < nOutLiners; ++i )
        outLiners.push_back( sortedEvents[i] );

    sortedEvents.erase( sortedEvents.begin(),
                        sortedEvents.begin() + nOutLiners );
    return outLiners;
}
