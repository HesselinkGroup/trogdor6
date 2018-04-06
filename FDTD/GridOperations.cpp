/*
 *  GridOperations.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/21/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "GridOperations.h"

using namespace std;


GridOperations::
GridOperations(GridDescPtr gridDesc,
    const std::map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids) :
    mGridDescription(gridDesc)
{
}


void GridOperations::
setPointers(map<int, Pointer<GridFields> > & fields)
{
    assert(fields.count(mGridDescription->id()));
    mFields = fields[mGridDescription->id()];
    
    handleSetPointers(fields);
}


