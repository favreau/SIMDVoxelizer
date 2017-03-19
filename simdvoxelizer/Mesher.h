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

#ifndef _Mesher_h_
#define _Mesher_h_

#include <memory>

// 3d surface mesh generation
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

// Surface mesh simplificatoin
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

// Sefl-inrtersection test
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include "Octree.h"

// Surface mesh simplificatoin
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> SurfaceMesh;

// 3d surface mesh generation typedefs
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef std::function<FT ( const Point_3& )> Function;
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

// Self-inrtersection test
typedef boost::graph_traits<SurfaceMesh>::face_descriptor face_descriptor;

class Mesher
{
public:
    Mesher( const std::string& inputFile, float alpha, const float angularBound,
            const float radialBound, const float distanceBound,
            const float leafSize, const float minRadius, const float decimationRatio,
            const float epsilon, const uint32_t targetFaceNumber );

    void mesh( const std::string& outputFile );

private:

    std::vector< Event > _extractOutLiners( std::vector< Event >& sortedEvents,
                                            uint32_t nOutLiners );
    float _computeFieldFromOctree( const OctreeNode* node, glm::vec3 location );
    FT _getFieldValue( const Point_3 p );

    std::vector< Event > _outLiners;

    uint32_t _nCgalQuery;
    uint32_t _targetFaceNumber;

    float _maxEventValue;
    float _decimationRatio;
    float _boundingBoxDiameter;
    float _octreeLeafSize;
    float _alpha; // multiply the radius

    float _angularBound;
    float _radialBound;
    float _distanceBound;
    float _epsilon;

    glm::vec3 _somaPosition;

    std::string _outputFile;

    std::unique_ptr< Octree > _octree;

};

#endif //_Mesher_h_
