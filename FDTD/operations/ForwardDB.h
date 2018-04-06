/*
 *  ForwardDB.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/22/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _UPDATEDB_
#define _UPDATEDB_

#include <iostream>
#include <map>
#include "Pointer.h"
#include "../FieldEnum.h"
#include "../GridFields.h"
#include "../Operation.h"

struct RunlineDB
{
    RunlineDB() {}
    RunlineDB(long f, long gj1, long gj2, long gk1, long gk2, long length) :
        fi(f),
        gjLow(gj1),
        gjHigh(gj2),
        gkLow(gk1),
        gkHigh(gk2),
        mLength(length)
    {
        assert(length > 0);
    }
    long length() const { return mLength; }
    
    long gjLow, gjHigh, gkLow, gkHigh;
    long fi;
    long mLength;
};
inline std::ostream & operator<<(std::ostream & str, const RunlineDB & rl)
{
    str << "[fi " << rl.fi << " gj " << rl.gjLow << " " << rl.gjHigh
        << " gk " << rl.gkLow << " " << rl.gkHigh
        << " length " << rl.length() << "]";
    return str;
}

struct RunlineDB_JM : public RunlineDB
{
    RunlineDB_JM() {}
    RunlineDB_JM(long f, long gj1, long gj2, long gk1, long gk2, long curr,
        long length) :
        RunlineDB(f, gj1, gj2, gk1, gk2, length),
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
        << " gk " << rl.gkLow << " " << rl.gkHigh << " current " << rl.current
        << " length " << rl.length()<< "]";
    return str;
}


class ForwardDB : public Operation
{
public:
    ForwardDB()
    {
    }
    
    ForwardDB(Field db, Field hej, Field hek,
        Precision::Float curlCoeff1, Precision::Float curlCoeff2,
        const std::vector<RunlineDB> & runlines) :
        mFieldDB(db), mFieldEHj(hej), mFieldEHk(hek),
        mRunlines(runlines),
        m_c1(curlCoeff1),
        m_c2(curlCoeff2)
    {
        name("DB");
    }
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields)
    {
        mHeadDBUpdate = currentGridFields.field(mFieldDB);
        mHeadEHNeighbors_j = currentGridFields.field(mFieldEHj);
        mHeadEHNeighbors_k = currentGridFields.field(mFieldEHk);
        
        assert(mHeadDBUpdate);
        assert(mHeadEHNeighbors_j);
        assert(mHeadEHNeighbors_k);
    }
    
    unsigned long bytes() const
    {
        return mRunlines.size()*sizeof(RunlineDB);
    }
    
    void apply(long timestep, Precision::Float dt)
    {
        for (int rr = 0; rr < mRunlines.size(); rr++)
        {
            Precision::Float* f = mHeadDBUpdate + mRunlines[rr].fi;
            const Precision::Float* gjLow = mHeadEHNeighbors_j + mRunlines[rr].gjLow;
            const Precision::Float* gjHigh = mHeadEHNeighbors_j + mRunlines[rr].gjHigh;
            const Precision::Float* gkLow = mHeadEHNeighbors_k+ mRunlines[rr].gkLow;
            const Precision::Float* gkHigh = mHeadEHNeighbors_k + mRunlines[rr].gkHigh;
            
            Precision::Float c1 = m_c1, c2 = m_c2;
            for (int ll = 0; ll < mRunlines[rr].length(); ll++)
            {
                *f = *f + c1*(*gkHigh - *gkLow) + c2*(*gjHigh - *gjLow);
                f++;
                gjLow++;
                gjHigh++;
                gkLow++;
                gkHigh++;
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
    std::vector<RunlineDB> mRunlines;
    Precision::Float m_c1;
    Precision::Float m_c2;
    Precision::Float* mHeadEHNeighbors_j; // for updating Dx, this is Hy
    Precision::Float* mHeadEHNeighbors_k;
    Precision::Float* mHeadDBUpdate; // for updating D, this is D
};

class ForwardDB_JM : public Operation
{
public:
    ForwardDB_JM()
    {
    }
    
    ForwardDB_JM(Field db, Field hej, Field hek, Field jm, 
        Precision::Float curlCoeff1, Precision::Float curlCoeff2, Precision::Float currentCoeff,
        const std::vector<RunlineDB_JM> & runlines) :
        mFieldDB(db), mFieldEHj(hej), mFieldEHk(hek), mFieldCurrent(jm),
        mRunlines(runlines),
        m_c1(curlCoeff1),
        m_c2(curlCoeff2),
        m_c3(currentCoeff),
        mHeadEHNeighbors_j(0L),
        mHeadEHNeighbors_k(0L),
        mHeadDBUpdate(0L),
        mHeadCurrent(0L)
    {
        name("DB+JM");
    }
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields)
    {
        mHeadDBUpdate = currentGridFields.field(mFieldDB);
        mHeadEHNeighbors_j = currentGridFields.field(mFieldEHj);
        mHeadEHNeighbors_k = currentGridFields.field(mFieldEHk);
        mHeadCurrent = currentGridFields.field(mFieldCurrent);
    }
    
    unsigned long bytes() const
    {
        return mRunlines.size()*sizeof(RunlineDB_JM);
    }
    
    void apply(long timestep, Precision::Float dt)
    {
        for (int rr = 0; rr < mRunlines.size(); rr++)
        {
            Precision::Float* f = mHeadDBUpdate + mRunlines[rr].fi;
            const Precision::Float* gjLow = mHeadEHNeighbors_j + mRunlines[rr].gjLow;
            const Precision::Float* gjHigh = mHeadEHNeighbors_j + mRunlines[rr].gjHigh;
            const Precision::Float* gkLow = mHeadEHNeighbors_k + mRunlines[rr].gkLow;
            const Precision::Float* gkHigh = mHeadEHNeighbors_k + mRunlines[rr].gkHigh;
            const Precision::Float* current = mHeadCurrent + mRunlines[rr].current;
            
            Precision::Float c1 = m_c1, c2 = m_c2, c3 = m_c3;
            for (int ll = 0; ll < mRunlines[rr].length(); ll++)
            {
                *f = *f + c1*(*gkHigh - *gkLow) + c2*(*gjHigh - *gjLow)
                    + c3*(*current);
                
                f++;
                gjLow++;
                gjHigh++;
                gkLow++;
                gkHigh++;
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
    std::vector<RunlineDB_JM> mRunlines;
    Precision::Float m_c1;
    Precision::Float m_c2;
    Precision::Float m_c3;
    Precision::Float* mHeadEHNeighbors_j; // for updating Dx, this is Hy
    Precision::Float* mHeadEHNeighbors_k;
    Precision::Float* mHeadDBUpdate; // for updating H, this is the neighbor E
    Precision::Float* mHeadCurrent;
};



#endif
