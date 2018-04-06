/*
 *  ForwardDB_2d.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/30/18.
 *  Copyright 2018 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _FORWARDDB_2D_
#define _FORWARDDB_2D_

#include <iostream>
#include <map>
#include "Pointer.h"
#include "../FieldEnum.h"
#include "../GridFields.h"
#include "../Operation.h"

struct RunlineDB_2d
{
    RunlineDB_2d() {}
    RunlineDB_2d(long f, long gj1, long gj2, long length) :
        fi(f),
        gjLow(gj1),
        gjHigh(gj2),
        mLength(length)
    {
        assert(length > 0);
    }
    long length() const { return mLength; }
    
    long gjLow, gjHigh;
    long fi;
    long mLength;
};
inline std::ostream & operator<<(std::ostream & str, const RunlineDB & rl)
{
    str << "[fi " << rl.fi << " gj " << rl.gjLow << " " << rl.gjHigh
        << " length " << rl.length() << "]";
    return str;
}

struct RunlineDB_JM_2d : public RunlineDB_2d
{
    RunlineDB_JM_2d() {}
    RunlineDB_JM_2d(long f, long gj1, long gj2, long curr, long length) :
        RunlineDB(f, gj1, gj2, length),
        current(curr)
    {
        assert(length > 0);
    }
    long length() const { return mLength; }
    
    long current;
};
inline std::ostream & operator<<(std::ostream & str, const RunlineDB_JM & rl)
{
    str << "[fi " << rl.fi << " gj " << rl.gjLow << " " << rl.gjHigh
        << " current " << rl.current
        << " length " << rl.length()<< "]";
    return str;
}


class ForwardDB_2d : public Operation
{
public:
    ForwardDB_2d()
    {
    }
    
    ForwardDB_2d(Field db, Field hej, Precision::Float curlCoeff1, const std::vector<RunlineDB> & runlines) :
        mFieldDB(db), mFieldEHj(hej),
        mRunlines(runlines),
        m_curlCoeff(curlCoeff1)
    {
        name("DB 2d");
    }
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields)
    {
        mHeadDBUpdate = currentGridFields.field(mFieldDB);
        mHeadEHNeighbors_j = currentGridFields.field(mFieldEHj);
        
        assert(mHeadDBUpdate);
        assert(mHeadEHNeighbors_j);
    }
    
    unsigned long bytes() const
    {
        return mRunlines.size()*sizeof(RunlineDB_2d);
    }
    
    void apply(long timestep, Precision::Float dt)
    {
        for (int rr = 0; rr < mRunlines.size(); rr++)
        {
            Precision::Float* f = mHeadDBUpdate + mRunlines[rr].fi;
            const Precision::Float* gjLow = mHeadEHNeighbors_j + mRunlines[rr].gjLow;
            const Precision::Float* gjHigh = mHeadEHNeighbors_j + mRunlines[rr].gjHigh;
            
            Precision::Float c = m_curlCoeff;
            for (int ll = 0; ll < mRunlines[rr].length(); ll++)
            {
//                *f = *f + c1*(*gkHigh - *gkLow) + c2*(*gjHigh - *gjLow);
                *f = *f + c*(*gjHigh - *gjLow);
                f++;
                gjLow++;
                gjHigh++;
            }
        }
    }
    
    void printRunlines(std::ostream & str) const
    {
        for (int nn = 0; nn < mRunlines.size(); nn++)
            str << mRunlines[nn] << "\n";
    }
    
    long numCells() const
    {
        long cells = 0;
        for (int nn = 0; nn < mRunlines.size(); nn++)
            cells += mRunlines[nn].length();
        return cells;
    }
    
    void allocate()
    {
    }

private:
    Field mFieldDB, mFieldEHj, mFieldEHk;
    std::vector<RunlineDB_2d> mRunlines;
    Precision::Float m_curlCoeff;
    Precision::Float* mHeadEHNeighbors_j; // for updating Dx, this is Hy or Hk
    Precision::Float* mHeadDBUpdate; // for updating D, this is D
};

class ForwardDB_JM_2d : public Operation
{
public:
    ForwardDB_JM_2d()
    {
    }
    
    ForwardDB_JM_2d(Field db, Field hej, Field jm,
        Precision::Float curlCoeff, Precision::Float currentCoeff,
        const std::vector<RunlineDB_JM_2d> & runlines) :
        mFieldDB(db), mFieldEHj(hej), mFieldCurrent(jm),
        mRunlines(runlines),
        m_curlCoeff(curlCoeff),
        m_currentCoeff(currentCoeff),
        mHeadEHNeighbors_j(0L),
        mHeadDBUpdate(0L),
        mHeadCurrent(0L)
    {
        name("DB+JM 2d");
    }
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields)
    {
        mHeadDBUpdate = currentGridFields.field(mFieldDB);
        mHeadEHNeighbors_j = currentGridFields.field(mFieldEHj);
        mHeadCurrent = currentGridFields.field(mFieldCurrent);
    }
    
    unsigned long bytes() const
    {
        return mRunlines.size()*sizeof(RunlineDB_JM_2d);
    }
    
    void apply(long timestep, Precision::Float dt)
    {
        for (int rr = 0; rr < mRunlines.size(); rr++)
        {
            Precision::Float* f = mHeadDBUpdate + mRunlines[rr].fi;
            const Precision::Float* gjLow = mHeadEHNeighbors_j + mRunlines[rr].gjLow;
            const Precision::Float* gjHigh = mHeadEHNeighbors_j + mRunlines[rr].gjHigh;
            const Precision::Float* current = mHeadCurrent + mRunlines[rr].current;
            
            Precision::Float c = m_curlCoeff, cCurrent = m_currentCoeff;
            for (int ll = 0; ll < mRunlines[rr].length(); ll++)
            {
                *f = *f + c*(*gjHigh - *gjLow) + cCurrent*(*current);
                
                f++;
                gjLow++;
                gjHigh++;
                current++;
            }
        }
    }
    
    void printRunlines(std::ostream & str) const
    {
        for (int nn = 0; nn < mRunlines.size(); nn++)
            str << mRunlines[nn] << "\n";
    }
    
    long numCells() const
    {
        long cells = 0;
        for (int nn = 0; nn < mRunlines.size(); nn++)
            cells += mRunlines[nn].length();
        return cells;
    }
    
    void allocate()
    {
    }

private:
    Field mFieldDB, mFieldEHj, mFieldEHk, mFieldCurrent;
    std::vector<RunlineDB_JM_2d> mRunlines;
    Precision::Float m_curlCoeff;
    Precision::Float m_currentCoeff;
    Precision::Float* mHeadEHNeighbors_j; // for updating Dx, this is Hy
    Precision::Float* mHeadEHNeighbors_k;
    Precision::Float* mHeadDBUpdate; // for updating H, this is the neighbor E
    Precision::Float* mHeadCurrent;
};



#endif
