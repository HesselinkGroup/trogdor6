/*
 *  AdjointDB.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 7/15/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "AdjointDB.h"

#include "Log.h"
#include <cmath>

using namespace std;

AdjointDB::
AdjointDB() :
    mCircularBufferLag(0)
{
}

AdjointDB::
AdjointDB(Field db, Field eh,
    int ehOrder,
    const std::vector<RunlineAdjointDB> & runlines,
    const std::vector<Precision::RationalFunction> & coefficients,
    Precision::Float dt) :
    mCircularBufferLag(0),
    mFieldDB(db),
    mFieldEH(eh),
    mRunlines(runlines),
    mOrderEH(ehOrder)
{
    mCoefficientsEH.resize(mOrderEH+1);
    for (int eh = 0; eh <= mOrderEH; eh++)
    {
        // Numerator of 1/eps is [b0/a0, b1/a0, ...].
        // The coefficients should be divided by eps0 or mu0.
        mCoefficientsEH[eh] = std::vector<Precision::Float>(coefficients.size(), 0.0);
        for (int mm = 0; mm < coefficients.size(); mm++)
        {
            mCoefficientsEH[eh][mm] = coefficients[mm].numerator()[eh] * dt;
        }
    }
    
    // Check that a0 == 1.
    for (int mm = 0; mm < coefficients.size(); mm++)
        assert(fabs(coefficients[mm].denominator()[0] - 1.0) < 1e-4);
}

void AdjointDB::
allocate()
{
    long totalLength = 0;
    for (int nn = 0; nn < mRunlines.size(); nn++)
        totalLength += mRunlines[nn].length();
    
    mBufferEH.resize(mOrderEH);
    for (int eh = 0; eh < mBufferEH.size(); eh++)
        mBufferEH[eh] = std::vector<Precision::Float>(totalLength, 0.0);
    
//    LOG << "allocated.\n";
}

void AdjointDB::
setPointers(GridFields & currentGridFields,
    std::map<int, Pointer<GridFields> > & allGridFields)
{
    mHeadDB = currentGridFields.field(mFieldDB);
    mHeadEH = currentGridFields.field(mFieldEH);
    
    assert(mHeadDB != 0L);
    assert(mHeadEH != 0L);
}

unsigned long AdjointDB::
bytes() const
{
    long totalBytes = 0;
    totalBytes += mRunlines.size()*sizeof(RunlineAdjointDB);
    
    totalBytes += sizeof(std::vector<Precision::Float>)*(
        mCoefficientsEH.size() + mBufferEH.size());
    
    for (int nn = 0; nn < mCoefficientsEH.size(); nn++)
        totalBytes += sizeof(Precision::Float)*mCoefficientsEH[nn].size();
    for (int nn = 0; nn < mBufferEH.size(); nn++)
        totalBytes += sizeof(Precision::Float)*mBufferEH[nn].size();
    
    return totalBytes;
}

void AdjointDB::
apply(long timestep, Precision::Float dt)
{
    // only handling zeroth-order material now
    for (int rr = 0; rr < mRunlines.size(); rr++)
    {
        Precision::Float* db = mHeadDB + mRunlines[rr].db;
        Precision::Float* eh = mHeadEH + mRunlines[rr].eh;
        const Precision::Float* constEH[10];
        Precision::Float* oldEH[9];
        
        for (int nn = 0; nn <= mOrderEH; nn++)
            constEH[nn] = &(mCoefficientsEH[nn].at(mRunlines[rr].coeff));
        for (int nn = 0; nn < mOrderEH; nn++)
        {
            oldEH[nn] = &(bufferEH(nn).at(mRunlines[rr].aux));
            assert(mRunlines[rr].aux + mRunlines[rr].length() - 1 <
                bufferEH(nn).size());
        }
        
        int stride = mRunlines[rr].coeffStride;
        
        for (int ll = 0; ll < mRunlines[rr].length(); ll++)
        {
            // "D old" = "D new" + ... (counterintuitive!)
            
            *db = *db + *constEH[0] * (*eh);
            for (int nn = 1; nn <= mOrderEH; nn++)
                *db += *constEH[nn] * (*oldEH[nn-1]);
            
            if (mOrderEH > 0)
                *(oldEH[mOrderEH-1]) = *eh;
            for (int nn = 0; nn <= mOrderEH; nn++)
                constEH[nn] += stride;
            for (int nn = 0; nn < mOrderEH; nn++)
                oldEH[nn]++;
            
            db++;
            eh++;
        }
    }
    mCircularBufferLag++;
}

void AdjointDB::
printRunlines(std::ostream & str) const
{
    for (int nn = 0; nn < mRunlines.size(); nn++)
        str << mRunlines[nn] << "\n";
}



AdjointDB_JM::
AdjointDB_JM() :
    mCircularBufferLag(0)
{
}

AdjointDB_JM::
AdjointDB_JM(Field db, Field eh, Field jm,
    int ehOrder,
    Precision::Float currentCoefficient,
    const std::vector<RunlineAdjointDB_JM> & runlines,
    const std::vector<Precision::RationalFunction> & coefficients,
    Precision::Float dt) :
    mCircularBufferLag(0),
    mFieldDB(db),
    mFieldEH(eh),
    mFieldCurrent(jm),
    mRunlines(runlines),
    mCurrentCoefficient(currentCoefficient),
    mOrderEH(ehOrder)
{
    mCoefficientsEH.resize(mOrderEH+1);
    for (int eh = 0; eh <= mOrderEH; eh++)
    {
        // Numerator of 1/eps is [b0/a0, b1/a0, ...].
        // The coefficients should be divided by eps0 or mu0.
        mCoefficientsEH[eh] = std::vector<Precision::Float>(coefficients.size(), 0.0);
        for (int mm = 0; mm < coefficients.size(); mm++)
        {
            mCoefficientsEH[eh][mm] = coefficients[mm].numerator()[eh] * dt;
        }
    }
    
    // Check that a0 == 1.
    for (int mm = 0; mm < coefficients.size(); mm++)
        assert(fabs(coefficients[mm].denominator()[0] - 1.0) < 1e-4);
}

void AdjointDB_JM::
allocate()
{
    long totalLength = 0;
    for (int nn = 0; nn < mRunlines.size(); nn++)
        totalLength += mRunlines[nn].length();
    
    mBufferEH.resize(mOrderEH);
    for (int eh = 0; eh < mOrderEH; eh++)
        mBufferEH[eh] = std::vector<Precision::Float>(totalLength, 0.0);
    
//    LOG << "allocated.\n";
}

void AdjointDB_JM::
setPointers(GridFields & currentGridFields,
    std::map<int, Pointer<GridFields> > & allGridFields)
{
    mHeadEH = currentGridFields.field(mFieldEH);
    mHeadDB = currentGridFields.field(mFieldDB);
    mHeadJM = currentGridFields.field(mFieldCurrent);
}

unsigned long AdjointDB_JM::
bytes() const
{
    long totalBytes = 0;
    totalBytes += mRunlines.size()*sizeof(RunlineAdjointDB);
    
    totalBytes += sizeof(std::vector<Precision::Float>)*(
        mCoefficientsEH.size() + mBufferEH.size());
    
    for (int nn = 0; nn < mCoefficientsEH.size(); nn++)
        totalBytes += sizeof(Precision::Float)*mCoefficientsEH[nn].size();
    for (int nn = 0; nn < mBufferEH.size(); nn++)
        totalBytes += sizeof(Precision::Float)*mBufferEH[nn].size();
    
    return totalBytes;
}

void AdjointDB_JM::
apply(long timestep, Precision::Float dt)
{
    // only handling zeroth-order material now
    for (int rr = 0; rr < mRunlines.size(); rr++)
    {
        Precision::Float* db = mHeadDB + mRunlines[rr].db;
        const Precision::Float* eh = mHeadEH + mRunlines[rr].eh;
        const Precision::Float* current = mHeadJM + mRunlines[rr].current;
        
        const Precision::Float* constEH[10];
        Precision::Float* oldEH[9];
        
        for (int nn = 0; nn <= mOrderEH; nn++)
            constEH[nn] = &(mCoefficientsEH[nn].at(mRunlines[rr].coeff));
        for (int nn = 0; nn < mOrderEH; nn++)
        {
            oldEH[nn] = &(bufferEH(nn).at(mRunlines[rr].aux));
            assert(mRunlines[rr].aux + mRunlines[rr].length() - 1 <
                bufferEH(nn).size());
        }
        
        int stride = mRunlines[rr].coeffStride;
        
        Precision::Float currentCoeff = mCurrentCoefficient;
        for (int ll = 0; ll < mRunlines[rr].length(); ll++)
        {
//            Precision::Float prevDB_unused = *db;
            
            // currentCoeff should be +1 in an all cases.
            // constEH is b/a0/eps0 or d/c0/mu0.
            *db = *db + *constEH[0] * (*eh) + currentCoeff*(*current);
            
            for (int nn = 1; nn <= mOrderEH; nn++)
                *db += *constEH[nn] * (*oldEH[nn-1]);
            
            if (mOrderEH > 0)
                *(oldEH[mOrderEH-1]) = *eh;
            for (int nn = 0; nn <= mOrderEH; nn++)
                constEH[nn] += stride;
            for (int nn = 0; nn < mOrderEH; nn++)
                oldEH[nn]++;
            
            db++;
            eh++;
            current++;
        }
    }
    mCircularBufferLag++;
}

void AdjointDB_JM::
printRunlines(std::ostream & str) const
{
    for (int nn = 0; nn < mRunlines.size(); nn++)
        str << mRunlines[nn] << "\n";
}

