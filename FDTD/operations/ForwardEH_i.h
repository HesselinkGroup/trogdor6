/*
 *  ForwardEH.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/22/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _UPDATEEHI_
#define _UPDATEEHI_

#include "../PrecisionRationalFunction.h"
#include "../FieldEnum.h"
#include "../GridFields.h"
#include "../Operation.h"
#include <vector>
#include <map>
#include "Pointer.h"

struct RunlineSumEH_i
{
    RunlineSumEH_i(long inEHi, long inEHii, long inEHij, long inEHik, long len) :
        ehi(inEHi), ehii(inEHii), ehij(inEHij), ehik(inEHik), mLength(len)
    {}
    
    long length() const { return mLength; }
    long ehi, ehii, ehij, ehik, mLength;
};
std::ostream & operator<<(std::ostream & str, const RunlineSumEH_i & rl);

// This operation handles the summation Ex = Exx + Exy + Exz.
class SumEH_i : public Operation
{
public:
    SumEH_i();
    SumEH_i(Field ehi, Field ehii, Field ehij, Field ehik,
        const std::vector<RunlineSumEH_i> & runlines);
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields);
    unsigned long bytes() const;
    void apply(long timestep, Precision::Float dt);
    void printRunlines(std::ostream & str) const;
    void allocate() {}
    
    long numCells() const
    {
        long cells = 0;
        for (int nn = 0; nn < mRunlines.size(); nn++)
            cells += mRunlines[nn].length();
        return cells;
    }
    
private:
    Field mField_i, mField_ii, mField_ij, mField_ik;
    std::vector<RunlineSumEH_i> mRunlines;
    Precision::Float *mHeadEHi, *mHeadEHii, *mHeadEHij, *mHeadEHik;
};

struct RunlineAnisotropicEH_i
{
    RunlineAnisotropicEH_i() {}
    RunlineAnisotropicEH_i(long inEH,
        long inDB0, long inDB1, long inDB2, long inDB3,
        long inAux, long inCoeff, long inCoeffStride, long length) :
        eh(inEH),
        db0(inDB0),
        db1(inDB1),
        db2(inDB2),
        db3(inDB3),
        aux(inAux),
        coeff(inCoeff),
        coeffStride(inCoeffStride),
        mLength(length)
    {
        assert(length > 0);
    }
    
    long length() const { return mLength; }
    
    long eh;
    long db0, db1, db2, db3;
    long aux;
    long coeff;
    long coeffStride;
    long mLength;
};
std::ostream & operator<<(std::ostream & str, const RunlineAnisotropicEH_i & rl);

// TODO: This class is almost entirely identical to ForwardEH, and so should
// TODO: probably be somehow consolidated with it.
class UpdateAnisotropicEH_i : public Operation
{
public:
    UpdateAnisotropicEH_i();
    UpdateAnisotropicEH_i(Field eh, Field db,
        int ehOrder, int dbOrder,
        const std::vector<RunlineAnisotropicEH_i> & runlines,
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
    
    std::vector<RunlineAnisotropicEH_i> mRunlines;
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
