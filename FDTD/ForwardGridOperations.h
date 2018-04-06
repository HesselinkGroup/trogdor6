/*
 *  ForwardGridOperations.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 7/7/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _FORWARDGRIDOPERATIONS_
#define _FORWARDGRIDOPERATIONS_

#include "GridOperations.h"
#include "Operation.h"

#include "operations/JMPML.h"
#include "operations/OutputEHDB.h"
#include "operations/Source.h"
#include "operations/HuygensSurface.h"
#include "operations/PeriodicBC.h"

class ForwardGridOperations : public GridOperations
{
public:
    ForwardGridOperations(GridDescPtr gridDesc,
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
    virtual void handleSetPointers(std::map<int, Pointer<GridFields> > & fields);
    
private:
    void clearJ();
    void clearM();
    
    void updateSourceE(long timestep, Precision::Float dt);
    void updateMSource(long timestep, Precision::Float dt);
    
    void updateSourceH(long timestep, Precision::Float dt);
    void updateJSource(long timestep, Precision::Float dt);
        
    void addBoundaryOutputs(GridDescription & gridDesc,
        const VoxelizedPartition & voxelGrid);
    void encodeOutputs(const GridDescription & gridDesc,
        const VoxelizedPartition & voxelGrid);
    void encodeSources(const GridDescription & gridDesc,
        const VoxelizedPartition & voxelGrid);
    void encodeCurrentSources(const GridDescription & gridDesc,
        const VoxelizedPartition & voxelGrid);
    void encodeHuygens(const GridDescription & gridDesc,
        const VoxelizedPartition & voxelGrid,
        const std::map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids);
    void encodePeriodicBCs(const VoxelizedPartition & voxelGrid);
    
    void encodeD(const VoxelizedPartition & voxelGrid);
    void encodeNondispersiveE(const VoxelizedPartition & voxelGrid);
    void encodeEi(const VoxelizedPartition & voxelGrid);
    void encodeEii(const VoxelizedPartition & voxelGrid);
    void encodeEij_0(const VoxelizedPartition & voxelGrid); // method "0"
    void encodeEij_i(const VoxelizedPartition & voxelGrid); // method "i"
    void encodeEij_j(const VoxelizedPartition & voxelGrid); // method "j"
    void encodeB(const VoxelizedPartition & voxelGrid);
    void encodeNondispersiveH(const VoxelizedPartition & voxelGrid);
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
    std::vector<Pointer<SourceEH> > mSourceEH;
    
    std::vector<Pointer<PeriodicBC> > mPeriodicBCE;
    std::vector<Pointer<PeriodicBC> > mPeriodicBCH;
    
    std::vector<Pointer<OutputEHDB> > mOutputs;
    
    // Benchmarking information
    // All times are in microseconds.
    double mTimeClearJM;
    std::vector<double> mTimeE;
    std::vector<double> mTimeH;
    std::vector<double> mTimeSourceJM;
    std::vector<double> mTimeSourceEH;
    std::vector<double> mTimeBCE;
    std::vector<double> mTimeBCH;
    std::vector<double> mTimeOutputs;
};




#endif
