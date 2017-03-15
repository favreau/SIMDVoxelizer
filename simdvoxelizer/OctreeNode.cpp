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

#include "OctreeNode.h"
#include <stdint.h>
#include <iostream>

OctreeNode::OctreeNode( const glm::vec3 center, const float size )
    : _maxValue( 0 )
    , _numberOfEvent( 0 )
    , _center( center )
    , _size( size )
{
    _halfDiagonal = sqrt( _size.x * _size.x + _size.y * _size.y +
                          _size.z * _size.z ) / 2.0f;
}

void OctreeNode::setChild( OctreeNode* child )
{
    _children.push_back( child );
}

void OctreeNode::addMaxValue( const float value )
{
    if( value > _maxValue )
        _maxValue = value;

    ++_numberOfEvent;
}

void OctreeNode::addEvent(const Event& event )
{
    addMaxValue( event.value );
    ++_numberOfEvent;
    _data.push_back( event );
}

const glm::vec3& OctreeNode::getCenter() const
{
    return _center;
}

const std::vector< Event > OctreeNode::getData() const
{
    if( _data.empty( ))
    {
        Event event = { glm::vec3(_center.x, _center.y, _center.z), _maxValue };
        return std::vector< Event >( 1, event );
    }

    return _data;
}

float OctreeNode::getHalfDiagonal() const
{
    return _halfDiagonal;
}

float OctreeNode::getNumberOfEvents() const
{
    return _numberOfEvent;
}

float OctreeNode::getMaxValue() const
{
    return _maxValue;
}

const std::vector< OctreeNode* >& OctreeNode::getChildren() const
{
    return _children;
}
