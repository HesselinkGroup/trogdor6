/*
 *  ForwardEH.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/22/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _UPDATEEH0_
#define _UPDATEEH0_

#include "../PrecisionRationalFunction.h"
#include "../FieldEnum.h"
#include "../GridFields.h"
#include "../Operation.h"
#include <vector>
#include <map>
#include "Pointer.h"

struct RunlineSumEH_0
{
    RunlineSumEH_0(long inEHi, long inEHii, long inEHij0, long inEHij1,
        long inEHik0, long inEHik1, long len) :
        ehi(inEHi), ehii(inEHii), ehij0(inEHij0), ehij1(inEHij1), 
        ehik0(inEHik0), ehik1(inEHik1), mLength(len)
    {}
    
    long length() const { return mLength; }
    long ehi, ehii, ehij0, ehij1, ehik0, ehik1, mLength;
};
std::ostream & operator<<(std::ostream & str, const RunlineSumEH_0 & rl);

// This operation handles the summation Ex = Exx + Exy + Exz.
class SumEH_0 : public Operation
{
public:
    SumEH_0();
    SumEH_0(Field ehi, Field ehii, Field ehij, Field ehik,
        const std::vector<RunlineSumEH_0> & runlines);
    
//    void setPointers(Precision::Float* headEHi, Precision::Float* headEHii,
//        Precision::Float* headEHij, Precision::Float* headEHik);
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
    std::vector<RunlineSumEH_0> mRunlines;
    Precision::Float *mHeadEHi, *mHeadEHii, *mHeadEHij, *mHeadEHik;
};



struct RunlineAnisotropicEH_0
{
    RunlineAnisotropicEH_0() {}
    RunlineAnisotropicEH_0(long inEH, long inDB0, long inDB1, long inAux,
        long inCoeff, long inCoeffStride, long length) :
        eh(inEH),
        db0(inDB0),
        db1(inDB1),
        aux(inAux),
        coeff(inCoeff),
        coeffStride(inCoeffStride),
        mLength(length)
    {
        assert(length > 0);
    }
    
    long length() const { return mLength; }
    
    long eh;
    long db0, db1;
    long aux;
    long coeff;
    long coeffStride;
    long mLength;
};
std::ostream & operator<<(std::ostream & str, const RunlineAnisotropicEH_0 & rl);

// TODO: This class is almost entirely identical to ForwardEH, and so should
// TODO: probably be somehow consolidated with it.
class UpdateAnisotropicEH_0 : public Operation
{
public:
    UpdateAnisotropicEH_0();
    UpdateAnisotropicEH_0(Field eh, Field db,
        int ehOrder, int dbOrder,
        const std::vector<RunlineAnisotropicEH_0> & runlines,
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
    
    std::vector<RunlineAnisotropicEH_0> mRunlines;
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
