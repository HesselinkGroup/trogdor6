/*
 *  EncodeOutput.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/23/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef T6_OUTPUTOPERATIONS
#define T6_OUTPUTOPERATIONS

#include "Pointer.h"
#include "operations/OutputEHDB.h"

class VoxelizedPartition;
class GridDescription;
class OutputDescription;

class EncodeOutput
{
public:
    EncodeOutput(const GridDescription & gridDesc,
        const VoxelizedPartition & vp);
    
    Pointer<OutputEHDB> encode(const OutputDescription & outputDesc) const;
    
    const GridDescription & gridDescription() const { return mGridDesc; }
    const VoxelizedPartition & voxelGrid() const { return mVoxelGrid; }
    
private:
    
    void encodeStandardOutput(OutputEHDB & output,
        const OutputDescription & outputDesc,
        const OutputField & field) const;
    
    // for J and M basically.
    void encodeSparseOutput(OutputEHDB & output,
        const OutputDescription & outputDesc,
        const OutputField & field) const;
    
    void encodeInterpolatedOutput(OutputEHDB & output,
        const OutputDescription & outputDesc,
        const OutputField & field) const;
    
    void encodeTensorialDiagonalOutput(OutputEHDB & output,
        const OutputDescription & outputDesc,
        const OutputField & field) const;
    
    void encodeTensorialOffDiagonalOutput(OutputEHDB & output,
        const OutputDescription & outputDesc,
        const OutputField & field) const;
    
    const GridDescription & mGridDesc;
    const VoxelizedPartition & mVoxelGrid;
};

Pointer<OutputEHDB> encodeOutput(const GridDescription & gridDesc,
    const VoxelizedPartition & voxelGrid,
    const OutputDescription & outputDesc);


#endif
