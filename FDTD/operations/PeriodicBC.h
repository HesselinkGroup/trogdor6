/*
 *  PeriodicBC.h
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 2/28/12.
 *  Copyright 2012 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef PERIODICBC_H
#define PERIODICBC_H

#include "../Operation.h"
#include "../FieldEnum.h"

struct RunlinePeriodicBC
{
    RunlinePeriodicBC() {}
    RunlinePeriodicBC(long read, long write, long length) :
        read(read),
        write(write),
        mLength(length)
    {
        assert(length > 0);
    }
    
    long length() const { return mLength; }
    
    long read, write, mLength;
};

class PeriodicBC : public Operation
{
public:
    PeriodicBC() {}
    PeriodicBC(Field f, const std::vector<RunlinePeriodicBC> & runlines) :
        mField(f),
        mHead(0L),
        mRunlines(runlines)
    {
    }
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields)
    {
        mHead = currentGridFields.field(mField);
    }
    
    unsigned long bytes() const
    {
        long totalBytes = 0;
        totalBytes += mRunlines.size() * sizeof(RunlinePeriodicBC);
        return totalBytes;
    }
    
    void apply(long timestep, Precision::Float dt)
    {
        for (int rr = 0; rr < mRunlines.size(); rr++)
        {
            Precision::Float* write = mHead + mRunlines[rr].write;
            const Precision::Float* read = mHead + mRunlines[rr].read;
            
            for (int ll = 0; ll < mRunlines[rr].length(); ll++)
                *write++ = *read++;
        }
    }
    
    long numCells() const
    {
        long cells = 0;
        for (int nn = 0; nn < mRunlines.size(); nn++)
            cells += mRunlines[nn].length();
        return cells;
    }
    
    void allocate() {};
    
private:
    Field mField;
    Precision::Float* mHead;
    std::vector<RunlinePeriodicBC> mRunlines;
};







#endif
