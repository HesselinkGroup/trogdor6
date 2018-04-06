/*
 *  AdjointGridOperations.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 7/7/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _ADJOINTGRIDOPERATIONS_
#define _ADJOINTGRIDOPERATIONS_

#include "GridOperations.h"
#include "Operation.h"

#include "operations/JMPML.h"

#include "operations/AdjointEH.h"
#include "operations/AdjointDB.h"

#include "operations/OutputEHDB.h"
#include "operations/Source.h"
#include "operations/PeriodicBC.h"

class AdjointGridOperations : public GridOperations
{
public:
    AdjointGridOperations(GridDescPtr gridDesc,
        const std::map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids);
    
    virtual void allocate();
    
    virtual void firstHalfStep(long timestep, Precision::Float dt);
    virtual void secondHalfStep(long timestep, Precision::Float dt);
    virtual void output(long timestep);
    
    virtual long bytes() const;
    virtual void printRunlines() const;
    
    virtual void initPerformance();
    virtual void printPerformance(double numT) const;
    
protected:
    virtual void handleSetPointers(
        std::map<int, Pointer<GridFields> > & fields);
private:
    void clearJ();
    void clearM();
    
//    void updateD(long timestep, Precision::Float dt);
//    void updateE(long timestep, Precision::Float dt);
    void updateMSource(long timestep, Precision::Float dt);
//    void updateMPML(long timestep, Precision::Float dt);
//    void updateMTFSF(long timestep, Precision::Float dt);
//    void updateGhostE(long timestep, Precision::Float dt);
//    void updateMPIE(long timestep, Precision::Float dt);
//    
//    void updateB(long timestep, Precision::Float dt);
//    void updateH(long timestep, Precision::Float dt);
    void updateJSource(long timestep, Precision::Float dt);
//    void updateJPML(long timestep, Precision::Float dt);
//    void updateJTFSF(long timestep, Precision::Float dt);
//    void updateGhostH(long timestep, Precision::Float dt);
//    void updateMPIH(long timestep, Precision::Float dt);
    
    void addBoundaryOutputs(GridDescription & gridDesc,
        const VoxelizedPartition & voxelGrid);
    void encodeOutputs(const GridDescription & gridDesc,
        const VoxelizedPartition & voxelGrid);
    void encodeCurrentSources(const GridDescription & gridDesc,
        const VoxelizedPartition & voxelGrid);
    void encodeHuygens(const GridDescription & gridDesc,
        const VoxelizedPartition & voxelGrid,
        const std::map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids);
    void encodePeriodicBCs(const VoxelizedPartition & voxelGrid);
    
    void encodeD(const VoxelizedPartition & voxelGrid);
    void encodeE(const VoxelizedPartition & voxelGrid);
    void encodeB(const VoxelizedPartition & voxelGrid);
    void encodeH(const VoxelizedPartition & voxelGrid);
    void encodePMLJ(const VoxelizedPartition & voxelGrid);
    void encodePMLM(const VoxelizedPartition & voxelGrid);
    
    template<class Op>
    void pushUpdateE(Op* op)
    {
        mUpdatesE.push_back(Pointer<Operation>((Operation*)op));
    }
    template<class Op>
    void pushUpdateH(Op* op)
    {
        mUpdatesH.push_back(Pointer<Operation>((Operation*)op));
    }
    std::vector<Pointer<Operation> > mUpdatesE;
    std::vector<Pointer<Operation> > mUpdatesH;
    
    std::vector<Pointer<SourceJM> > mSourceJM;
    
    std::vector<Pointer<PeriodicBC> > mPeriodicBCD;
    std::vector<Pointer<PeriodicBC> > mPeriodicBCB;
    
    std::vector<Pointer<OutputEHDB> > mOutputs;
    
    // Benchmarking information
    // All times are in microseconds.
    double mTimeClearJM;
    std::vector<double> mTimeE;
    std::vector<double> mTimeH;
    std::vector<double> mTimeSourceJM;
    std::vector<double> mTimeBCD;
    std::vector<double> mTimeBCB;
    std::vector<double> mTimeOutputs;
};




#endif
