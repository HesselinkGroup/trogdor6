/*
 *  GridOperations.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/21/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _GRIDOPERATIONS_
#define _GRIDOPERATIONS_

#include <vector>
#include <map>
#include "SimulationDescription.h"
#include "GridFields.h"
#include "VoxelizedPartition.h"
#include "Pointer.h"

class GridOperations
{
public:
    GridOperations(GridDescPtr gridDesc,
        const std::map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids);
    
    virtual ~GridOperations() {}
    
    void setPointers(
        std::map<int, Pointer<GridFields> > & fields);
    
    virtual void allocate() = 0;
    
    virtual void firstHalfStep(long timestep, Precision::Float dt) = 0;
    virtual void secondHalfStep(long timestep, Precision::Float dt) = 0;
    virtual void output(long timestep) = 0;
    
    virtual long bytes() const = 0;
    virtual void printRunlines() const = 0;
    
    virtual void initPerformance() = 0;
    virtual void printPerformance(double numT) const = 0;
        
    GridFields & fields() { return *mFields; }
    GridDescPtr gridDescription() const { return mGridDescription; }
protected:
    virtual void handleSetPointers(
        std::map<int, Pointer<GridFields> > & fields) = 0;
    
private:
    Pointer<GridFields> mFields;
    GridDescPtr mGridDescription;
};
typedef Pointer<GridOperations> GridOperationsPtr;















#endif
