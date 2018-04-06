/*
 *  AdjointDB.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 7/8/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _UPDATEADJOINTDB_
#define _UPDATEADJOINTDB_

#include <iostream>
#include <vector>
#include <map>
#include "Pointer.h"
#include "../PrecisionRationalFunction.h"
#include "../FieldEnum.h"
#include "../GridFields.h"
#include "../Operation.h"

struct RunlineAdjointDB
{
    RunlineAdjointDB() {}
    RunlineAdjointDB(long inDB, long inEH, long inAux, long inCoeff,
        long inCoeffStride,
        long length) :
        db(inDB),
        eh(inEH),
        aux(inAux),
        coeff(inCoeff),
        coeffStride(inCoeffStride),
        mLength(length)
    {
        assert(length > 0);
    }
    
    long length() const { return mLength; }
    
    long db, eh, aux, coeff, coeffStride, mLength;
};
inline std::ostream & operator<<(std::ostream & str,
    const RunlineAdjointDB & rl)
{
    str << "[db " << rl.db << " eh " << rl.eh << " aux " << rl.aux
        << " coeff " << rl.coeff << "x" << rl.coeffStride
        << " length " << rl.length() << "]";
    return str;
}

struct RunlineAdjointDB_JM
{
    RunlineAdjointDB_JM() {}
    RunlineAdjointDB_JM(long inDB, long inEH, long inAux, long inCoeff,
        long inCoeffStride, long inCurrent, long length) :
        db(inDB),
        eh(inEH),
        aux(inAux),
        coeff(inCoeff),
        coeffStride(inCoeffStride),
        current(inCurrent),
        mLength(length)
    {
        assert(length > 0);
    }
    
    long length() const { return mLength; }
    
    long db, eh, aux, coeff, coeffStride, current, mLength;
};
inline std::ostream & operator<<(std::ostream & str,
    const RunlineAdjointDB_JM & rl)
{
    str << "[db " << rl.db << " eh " << rl.eh << " aux " << rl.aux
        << " coeff " << rl.coeff << "x" << rl.coeffStride
        << " current " << rl.current << " length " << rl.length() << "]";
    return str;
}


// AdjointDB:
//
//  D(n) - D(n+1) = b0*E(n) + b1*E(n+1) + ...
//
// bn are the *numerator* of the inverse permittivity

class AdjointDB : public Operation
{
public:
    AdjointDB();
    
    AdjointDB(Field db, Field eh,
        int ehOrder,
        const std::vector<RunlineAdjointDB> & runlines,
        const std::vector<Precision::RationalFunction> & coefficients,
        Precision::Float eps0_or_mu0);
    
    void allocate();
    
//    void setPointers(Precision::Float* headEH, Precision::Float* headDB);
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields);
    
    unsigned long bytes() const;
    
    void apply(long timestep, Precision::Float dt);
    void printRunlines(std::ostream & str) const;
    
    long numCells() const
    {
        long cells = 0;
        for (int nn = 0; nn < mRunlines.size(); nn++)
            cells += mRunlines[nn].length();
        return cells;
    }

private:
    /**
     * Return the buffer currently representing the given lag.  Because the
     * buffers are circular, the correct buffer will depend on the current
     * timestep.
     */
    std::vector<Precision::Float> & bufferEH(unsigned long lag)
    {
        return mBufferEH[(mCircularBufferLag - lag + mBufferEH.size()) % mBufferEH.size()];
    }
    unsigned long mCircularBufferLag;
    
    Field mFieldDB, mFieldEH;
    std::vector<RunlineAdjointDB> mRunlines;
    std::vector<std::vector<Precision::Float> > mCoefficientsEH; // [order][run]
    std::vector<std::vector<Precision::Float> > mBufferEH;
    Precision::Float* mHeadDB;
    Precision::Float* mHeadEH;
//    Precision::Float m_eps0_or_mu0;
    int mOrderEH;
};


// AdjointDB with J_D:
//
//  D(n) - D(n+1) = b0*E(n) + b1*E(n+1) + J_D
//
// bn are the *numerator* of the inverse permittivity

class AdjointDB_JM : public Operation
{
public:
    AdjointDB_JM();
    
    AdjointDB_JM(Field db, Field eh, Field jm,
        int ehOrder,
        Precision::Float currentCoefficient,
        const std::vector<RunlineAdjointDB_JM> & runlines,
        const std::vector<Precision::RationalFunction> & coefficients,
        Precision::Float eps0_or_mu0);
    
    void allocate();
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields);
    
    unsigned long bytes() const;
    
    void apply(long timestep, Precision::Float dt);
    
    void printRunlines(std::ostream & str) const;
    
    long numCells() const
    {
        long cells = 0;
        for (int nn = 0; nn < mRunlines.size(); nn++)
            cells += mRunlines[nn].length();
        return cells;
    }

private:
    /**
     * Return the buffer currently representing the given lag.  Because the
     * buffers are circular, the correct buffer will depend on the current
     * timestep.
     */
    std::vector<Precision::Float> & bufferEH(unsigned long lag)
    {
        return mBufferEH[(mCircularBufferLag - lag + mBufferEH.size()) % mBufferEH.size()];
    }
    unsigned long mCircularBufferLag;
    
    Field mFieldDB, mFieldEH, mFieldCurrent;
    std::vector<RunlineAdjointDB_JM> mRunlines;
    std::vector<std::vector<Precision::Float> > mCoefficientsEH; // [order][run]
    std::vector<std::vector<Precision::Float> > mBufferEH;
    Precision::Float* mHeadDB;
    Precision::Float* mHeadEH;
    Precision::Float* mHeadJM;
    Precision::Float mCurrentCoefficient;
//    Precision::Float m_eps0_or_mu0;
    int mOrderEH;
};



#endif
