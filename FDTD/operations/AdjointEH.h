/*
 *  AdjointEH.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 7/8/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _UPDATEADJOINTEH_
#define _UPDATEADJOINTEH_

#include <iostream>
#include <vector>
#include <map>
#include "Pointer.h"
#include "../PrecisionRationalFunction.h"
#include "../FieldEnum.h"
#include "../GridFields.h"
#include "../Operation.h"

struct RunlineAdjointEH
{
    RunlineAdjointEH() {}
    RunlineAdjointEH(long inEH, long gj1, long gj2, long gk1, long gk2,
        long inAux, long inCoeff, long inCoeffStride, long length) :
        eh(inEH),
        gjLow(gj1), 
        gjHigh(gj2),
        gkLow(gk1),
        gkHigh(gk2),
        aux(inAux),
        coeff(inCoeff),
        coeffStride(inCoeffStride),
        mLength(length)
    {
        assert(length > 0);
    }
    
    long length() const { return mLength; }
    
    long eh, gjLow, gjHigh, gkLow, gkHigh, aux, coeff, coeffStride, mLength;
};
inline std::ostream & operator<<(std::ostream & str,
    const RunlineAdjointEH & rl)
{
    str << "[eh " << rl.eh << " gj " << rl.gjLow << " " << rl.gjHigh
        << " gk " << rl.gkLow << " " << rl.gkHigh << " aux " << rl.aux
        << " coeff " << rl.coeff << "x" << rl.coeffStride
        << " length " << rl.length() << "]";
    return str;
}

struct RunlineAdjointEH_JM : public RunlineAdjointEH
{
    RunlineAdjointEH_JM() {}
    RunlineAdjointEH_JM(long inEH, long gj1, long gj2, long gk1, long gk2,
        long inAux, long inCoeff, long inCoeffStride, long curr, long length) :
        RunlineAdjointEH(inEH, gj1, gj2, gk1, gk2, inAux, inCoeff,
            inCoeffStride, length),
        current(curr)
    {
    }
    
    long current;
};
inline std::ostream & operator<<(std::ostream & str,
    const RunlineAdjointEH_JM & rl)
{
    str << "[eh " << rl.eh << " gj " << rl.gjLow << " " << rl.gjHigh
        << " gk " << rl.gkLow << " " << rl.gkHigh << " aux " << rl.aux
        << " coeff " << rl.coeff << "x" << rl.coeffStride << " current "
        << rl.current << " length " << rl.length() << "]";
    return str;
}

class AdjointEH : public Operation
{
public:
    AdjointEH();
    
    AdjointEH(Field eh, Field dbj, Field dbk,
        int ehOrder,
        Precision::Float curlCoeff1, Precision::Float curlCoeff2,
        const std::vector<RunlineAdjointEH> & runlines,
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
    
    Field mField, mFieldNeighborDBj, mFieldNeighborDBk;
    std::vector<RunlineAdjointEH> mRunlines;
    std::vector<std::vector<Precision::Float> > mCoefficientsEH;
    std::vector<std::vector<Precision::Float> > mBufferEH;
    Precision::Float m_c1;
    Precision::Float m_c2;
    Precision::Float* mHeadEHUpdate;
    Precision::Float* mHeadDBNeighbors_j;
    Precision::Float* mHeadDBNeighbors_k;
//    Precision::Float m_eps0_or_mu0;
    int mOrderEH;
};


class AdjointEH_JM : public Operation
{
public:
    AdjointEH_JM();
    
    AdjointEH_JM(Field eh, Field dbj, Field dbk, Field jm,
        int ehOrder,
        Precision::Float curlCoeff1, Precision::Float curlCoeff2,
        Precision::Float currentCoefficient,
        const std::vector<RunlineAdjointEH_JM> & runlines,
        const std::vector<Precision::RationalFunction> & coefficients,
        Precision::Float eps0_or_mu0);
    
    void allocate();
    
//    void setPointers(Precision::Float* headEH, Precision::Float* headNeighbors_j,
//        Precision::Float* headNeighbors_k, Precision::Float* headCurrent);
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
        unsigned long bufferIndex =
            (mCircularBufferLag - lag + mBufferEH.size()) % mBufferEH.size();
//        std::cerr << "lag " << lag << " index " << bufferIndex << "\n";
        return mBufferEH[bufferIndex];
    }
    unsigned long mCircularBufferLag;
    
    Field mField, mFieldNeighborDBj, mFieldNeighborDBk, mFieldCurrent;
    std::vector<RunlineAdjointEH_JM> mRunlines;
    std::vector<std::vector<Precision::Float> > mCoefficientsEH;
    std::vector<std::vector<Precision::Float> > mBufferEH;
    Precision::Float m_c1;
    Precision::Float m_c2;
    Precision::Float m_c3;
    Precision::Float* mHeadEHUpdate;
    Precision::Float* mHeadDBNeighbors_j;
    Precision::Float* mHeadDBNeighbors_k;
    Precision::Float* mHeadCurrent;
//    Precision::Float m_eps0_or_mu0;
    int mOrderEH;
};




#endif
