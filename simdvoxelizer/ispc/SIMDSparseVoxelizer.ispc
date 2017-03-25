/* Copyright (c) 2015-2017, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Grigori Chevtchenko <grigori.chevtchenko@epfl.ch>
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

float treeWalker( const varying float x, const varying float y,
                  const varying float z, const varying float dist,
                  const varying float cutoff, const varying unsigned int64 index,
                  uniform unsigned int32 indexInput[],
                  uniform float dataInput[] )
{
    unsigned int32 begin = indexInput[ index * 2 ];
    unsigned int32 end = indexInput[ index * 2 + 1 ];

    if( begin == 0 && end == 0 )
    {
        // here we are garanteed to be in a leaf
        float d = dist;
        if( d < 0.5 )
            d = 0.5;

        return dataInput[ index * 4 + 3 ] / ( d * d );
    }

    float voxelValue = 0;
    for( unsigned int32 childIndex = begin; childIndex <= end; ++childIndex )
    {
        varying float centerX = dataInput[ index * 4 ];
        varying float centerY = dataInput[ index * 4 + 1];
        varying float centerZ = dataInput[ index * 4 + 2];
        varying float deltaX = x - centerX;
        varying float deltaY = y - centerY;
        varying float deltaZ = z - centerZ;

        varying const float childDist =
            sqrt( deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ );

        if( childDist >= cutoff )
            voxelValue += dataInput[ index * 4 + 3 ] / ( childDist * childDist );
            //voxelValue += dataInput[ index * 4 + 3 ] / childDist;
        else
            voxelValue += treeWalker( x, y, z, childDist, cutoff / 2.0,
                                      childIndex, indexInput, dataInput );
    }
    return voxelValue;
}

task void sparseVoxelize(
    uniform unsigned int32 zEnd,
    uniform float voxelSize,
    uniform unsigned int32 xspan, uniform unsigned int32 yspan, varying unsigned int32 z,
    uniform unsigned int32 octreeSize,
    uniform unsigned int32 volumeDimX,
    uniform unsigned int32 volumeDimY,
    uniform unsigned int32 volumeDimZ,
    uniform unsigned int32 indexInput[],
    uniform float dataInput[],
    uniform float output[])
{
    const uniform unsigned int32 xstart = taskIndex0 * xspan;
    const uniform unsigned int32 xend = min( xstart  + xspan, volumeDimX );
    const uniform unsigned int32 ystart = taskIndex1 * yspan;
    const uniform unsigned int32 yend = min( ystart  + yspan, volumeDimY );

    for( int y = ystart; y < yend; ++y )
    {
        for( int x = xstart; x < xend; ++x )
        {
            float cutoff = octreeSize * 5.0f;

            unsigned int64 index = x + y * volumeDimX + z * volumeDimX * volumeDimY;

            varying float dist = 0.f;
            varying unsigned int32 idx = 0;
            if(( index < volumeDimX * volumeDimY * zEnd ) && ( index < volumeDimX * volumeDimY * volumeDimZ ))
                output[ index ] = treeWalker( x, y, z, dist, cutoff, idx, indexInput, dataInput );
        }
    }
}

export void SIMDSparseVoxelizer_ispc(
    uniform unsigned int32 zBegin,
    uniform unsigned int32 zEnd,
    uniform unsigned int32 span,
    uniform float voxelSize,
    uniform unsigned int32 octreeSize,
    uniform unsigned int32 volumeDimX,
    uniform unsigned int32 volumeDimY,
    uniform unsigned int32 volumeDimZ,
    uniform unsigned int32 indexInput[],
    uniform float dataInput[],
    uniform float output[])
{
    foreach( z = zBegin ... zEnd + span - 1 )
    {
        const uniform unsigned int32 xspan = max( num_cores(), programCount );
        const uniform unsigned int32 yspan = span;

        launch[ ( volumeDimX + span - 1 ) / xspan, ( volumeDimY + span - 1) / yspan ]
            sparseVoxelize(
                zEnd,
                voxelSize,
                xspan, yspan, z,
                octreeSize,
                volumeDimX,
                volumeDimY,
                volumeDimZ,
                indexInput,
                dataInput,
                output );
    }
}
