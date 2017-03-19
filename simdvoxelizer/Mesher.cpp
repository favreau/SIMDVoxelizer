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
namespace PMP = CGAL::Polygon_mesh_processing;

bool compareByValue( const Event& e1, const Event& e2 )
{
    return e1.value > e2.value;
}

Mesher::Mesher( const std::string& inputFile, const float alpha,
                const float angularBound, const float radialBound,
                const float distanceBound, const float leafSize,
                const float minRadiusForced, float decimationRatio,
                const float epsilon, const uint32_t targetFaceNumber )
    : _nCgalQuery( 0 )
    , _targetFaceNumber( targetFaceNumber )
    , _maxEventValue( -std::numeric_limits<float>::max( ))
    , _decimationRatio( decimationRatio )
    , _alpha( alpha )
    , _epsilon( epsilon )
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
    std::cout << "\nBounding box      ["
              << min.x << "," << min.y << "," << min.z << "] ["
              << max.x << "," << max.y << "," << max.z << "]" << std::endl;

    std::cout << "Bounding values      [" << minValue << " ] ["
              << _maxEventValue << "]" << std::endl;

    _somaPosition = _outLiners[0].position;

    if( angularBound == -1.0f )
        _angularBound = 22.0f;
    else
        _angularBound = angularBound;

    if( distanceBound == -1.0f )
        _distanceBound = minRadiusForced * 0.15;
    else
        _distanceBound = distanceBound;

    if( radialBound == -1.0f )
        _radialBound = 20 * minRadiusForced;
    else
        _radialBound = radialBound;

    if( leafSize == -1.0f )
        _octreeLeafSize = 4.0;
    else
        _octreeLeafSize = leafSize;

    std::cout << "CGAL's surface bounding diameter: " << _boundingBoxDiameter << std::endl;
    std::cout << "CGAL's facet angular bound: " << _angularBound << std::endl;
    std::cout << "CGAL's surface Delaunay balls upper bound radius: " << _radialBound << std::endl;
    std::cout << "CGAL's distance bound between Delaunay ball center"
                 " and facet circumcenter: " << _distanceBound << std::endl;
    std::cout << "CGAL's dicotomy precision: " << _boundingBoxDiameter * _epsilon << std::endl;

    std::cout << "Octree leaf size: " << _octreeLeafSize << std::endl;
    //!!!!!!!! Need an optimal leaf size to be computed!!
    _octree.reset( new Octree( sortedEvents, _octreeLeafSize, min, max ));
}

void Mesher::mesh( const std::string& outputFile )
{
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
                                 _epsilon);

    CGAL::Surface_mesh_default_criteria_3<Tr> criteria( _angularBound,
                                                        _radialBound,
                                                        _distanceBound );
    // meshing surface
    std::cout<<"meshing..." << std::endl;
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
    std::cout<<"\nmeshing finished" << std::endl;

    std::string fileOriginal = outputFile + ".off";
    std::ofstream output( fileOriginal );
    CGAL::output_surface_facets_to_off( output, c2t3 );

    // ************* MESH SIMPLIFICATION **************

    if( _targetFaceNumber > 0 )
    {
        _decimationRatio = (float)_targetFaceNumber / (float)c2t3.number_of_facets();
    }

    if( _decimationRatio >= 1.0f )
        return;

    std::cout<<"decimating with ratio "<< _decimationRatio <<"..." << std::endl;
    SurfaceMesh mesh;
    CGAL::output_surface_facets_to_polyhedron( c2t3, mesh );

    bool intersecting = PMP::does_self_intersect(mesh,
        PMP::parameters::vertex_point_map( get( CGAL::vertex_point, mesh )));

    std::cout << (intersecting ? "There are self-intersections " : "There is no self-intersection ")
              << "before the decimation process."<<std::endl;

    SMS::Count_ratio_stop_predicate< SurfaceMesh > stop( _decimationRatio );

    int r = SMS::edge_collapse
                  ( mesh
                   ,stop
                   ,CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, mesh ))
                                     .halfedge_index_map  (get(CGAL::halfedge_external_index  , mesh ))
                                     .get_cost (SMS::LindstromTurk_cost<SurfaceMesh>())
                                     .get_placement(SMS::LindstromTurk_placement<SurfaceMesh>())
                  );
    std::cout<<"decimation finished" << std::endl;
    std::string fileDecimated = outputFile + "_decimated.off";
    std::ofstream outputDecimated( fileDecimated );
    outputDecimated << mesh;

    intersecting = PMP::does_self_intersect(mesh,
        PMP::parameters::vertex_point_map( get( CGAL::vertex_point, mesh )));

    std::cout << (intersecting ? "There are self-intersections" : "There is no self-intersection")
              << " after the decimation process."<< std::endl;

    std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
    PMP::self_intersections( mesh, std::back_inserter( intersected_tris ),
                             PMP::parameters::vertex_point_map( get( CGAL::vertex_point, mesh )));
    std::cout << intersected_tris.size() << " pairs of triangles intersect." << std::endl;


    std::cout << "Done. 2 files written: " << std::endl;
    std::cout << "original  [V:" << tr.number_of_vertices() << " F:" << c2t3.number_of_facets()
              << "] -> " << fileOriginal << std::endl;
    std::cout << "decimated [V:" << mesh.size_of_vertices() << " F:" << mesh.size_of_facets()
              << "] -> " << fileDecimated << std::endl;
}

FT Mesher::_getFieldValue( const Point_3 p )
{
    const glm::vec3 location( p.x(), p.y(), p.z( ));
    float value = 0;
    for( const Event& outLiner : _outLiners )
    {
        const glm::vec3 delta = glm::vec3( outLiner.position.x, outLiner.position.y,
                                           outLiner.position.z ) - location;
        float distSquared = delta.x * delta.x +
                            delta.y * delta.y +
                            delta.z * delta.z;
        // Avoid division by 0
        if( distSquared < 0.01f )
            distSquared = 0.01f;

        value += outLiner.value / ( distSquared * distSquared );
    }
    value += _computeFieldFromOctree( _octree->getRoot(), location );

    ++_nCgalQuery;
    if( _nCgalQuery % 50000 == 0 )
        std::cout <<"\rNumber of CGAL's oracle queries: "
                  << _nCgalQuery << std::flush;

    const float threshold = 1.0f;
    return threshold - value;
}

float Mesher::_computeFieldFromOctree( const OctreeNode* node, glm::vec3 location )
{
    std::vector< OctreeNode* > children = node->getChildren();
    if( children.empty( ))
    {
        std::vector< Event > events = node->getData();
        float value = 0;
        for( const Event& event : events )
        {
            const glm::vec3 delta = glm::vec3( event.position.x, event.position.y,
                                               event.position.z ) - location;
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
        const glm::vec3 delta = childNode->getCenter() - location;
        float childDist = sqrt( delta.x * delta.x + delta.y * delta.y +
                                delta.z * delta.z ) - childNode->getHalfDiagonal();
        if( childDist < 0 )
            childDist = 0;

        if( childDist <= _alpha * _maxEventValue )
            value += _computeFieldFromOctree( childNode, location );
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
