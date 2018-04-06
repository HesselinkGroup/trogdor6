/*
 *  ForwardEH.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/22/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _UPDATEEH_
#define _UPDATEEH_

#include "../PrecisionRationalFunction.h"
#include "../FieldEnum.h"
#include "../GridFields.h"
#include "../Operation.h"
#include <vector>
#include <map>
#include "Pointer.h"

struct RunlineEH
{
    RunlineEH() {}
    RunlineEH(long inEH, long inDB, long inAux, long inCoeff, long inCoeffStride, long length) :
        eh(inEH),
        db(inDB),
        aux(inAux),
        coeff(inCoeff),
        coeffStride(inCoeffStride),
        mLength(length)
    {
        assert(length > 0);
    }
    
    long length() const { return mLength; }
    
    long eh;
    long db;
    long aux;
    long coeff;
    long coeffStride;
    long mLength;
};
std::ostream & operator<<(std::ostream & str, const RunlineEH & rl);

class ForwardEH : public Operation
{
public:
    ForwardEH();
    
    // I take the coefficients to be permittivity or permeability.
    // The numerator is associated with E or H.
    // The denominator is associated with D or B.
    ForwardEH(Field eh, Field db, 
        int ehOrder, int dbOrder,
        const std::vector<RunlineEH> & runlines,
        const std::vector<Precision::RationalFunction> & coefficients,
        Precision::Float eps0_or_mu0);
    
    void allocate();
//    void setPointers(Precision::Float* headDB, Precision::Float* headEH);
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
    std::vector<Precision::Float> & bufferDB(unsigned long lag)
    {
        assert(lag > 0); // lag 0 is always available on the grid!
        return mBufferDB[(mCircularBufferLag - lag + mBufferDB.size()) % mBufferDB.size()];
    }
    std::vector<Precision::Float> & bufferEH(unsigned long lag)
    {
        return mBufferEH[(mCircularBufferLag - lag + mBufferEH.size()) % mBufferEH.size()];
    }
    unsigned long mCircularBufferLag;
    
    Field mField;
    Field mAuxField;
    
    std::vector<RunlineEH> mRunlines;
    std::vector<std::vector<Precision::Float> > mCoefficientsDB; // [order][coeff run]
    std::vector<std::vector<Precision::Float> > mCoefficientsEH; // [order][coeff run]
    std::vector<std::vector<Precision::Float> > mBufferDB; // [lag][cell]
    std::vector<std::vector<Precision::Float> > mBufferEH; // [lag][cell]
    Precision::Float* mHeadDB;
    Precision::Float* mHeadEH;
    int mOrderEH;
    int mOrderDB;
//    Precision::Float m_eps0_or_mu0; // either eps0 or mu0
};





#endif
