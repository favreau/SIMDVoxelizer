task void voxelize(
    uniform float step,
    uniform int xspan, uniform int yspan, varying int z,
    uniform int volumex, uniform int volumey, uniform int volumez,
    uniform float minx, uniform float miny, uniform float minz,
    uniform float squaredCutoff,
    uniform float input[], uniform int inputSize,
    uniform float output[])
{
    const uniform int xstart = taskIndex0 * xspan;
    const uniform int xend   = min(xstart  + xspan, volumex);
    const uniform int ystart = taskIndex1 * yspan;
    const uniform int yend   = min(ystart  + yspan, volumey);

    for( int y = ystart; y < yend; ++y )
    {
        for( int x = xstart; x < xend; ++x )
        {
            float acc = 0.f;
            for( uniform int i = 0; i < inputSize; i += 5 )
            {
                const float dx = input[i  ] - ( x * step + minx );
                const float dy = input[i+1] - ( y * step + miny );
                const float dz = input[i+2] - ( z * step + minz );
                const float dist = rcp( dx*dx + dy*dy + dz*dz );

                if( dist >= squaredCutoff )
                {
                    const float radius = input[i+3];
                    const float value = input[i+4];
                    acc += ( dist > radius*radius ) ? value * radius : value * dist;
                }
            }
            const int index = x + y * volumex + z * volumex * volumey;
            output[ index ] = acc;
        }
    }
}

export void SIMDVoxelizer_ispc(
    uniform float step,
    uniform int volumex, uniform int volumey, uniform int volumez,
    uniform float minx, uniform float miny, uniform float minz,
    uniform float cutOffDistance,
    uniform float input[], uniform int inputSize,
    uniform float output[])
{
    uniform float squaredCutoff = rcp( cutOffDistance * cutOffDistance );

    foreach( z = 0 ... volumez )
    {
        const uniform int xspan = max( num_cores(), programCount );
        const uniform int yspan = 32;

        launch [volumex/xspan, volumey/yspan]
            voxelize(
                step,
                xspan, yspan, z,
                volumex, volumey, volumez,
                minx, miny, minz,
                squaredCutoff,
                input, inputSize,
                output );
    }
}
