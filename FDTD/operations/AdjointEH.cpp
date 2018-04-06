/*
 *  AdjointEH.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 7/15/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "AdjointEH.h"

#include "Log.h"
#include <cmath>

using namespace std;

AdjointEH::
AdjointEH() :
    mCircularBufferLag(0)
{
}

AdjointEH::
AdjointEH(Field eh, Field dbj, Field dbk,
    int ehOrder,
    Precision::Float curlCoeff1, Precision::Float curlCoeff2,
    const std::vector<RunlineAdjointEH> & runlines,
    const std::vector<Precision::RationalFunction> & coefficients,
    Precision::Float eps0_or_mu0) :
    mCircularBufferLag(0),
    mField(eh),
    mFieldNeighborDBj(dbj),
    mFieldNeighborDBk(dbk),
    mRunlines(runlines),
    m_c1(curlCoeff1/eps0_or_mu0),
    m_c2(curlCoeff2/eps0_or_mu0),
    mOrderEH(ehOrder)
//    m_eps0_or_mu0(eps0_or_mu0),
{
    mCoefficientsEH.resize(mOrderEH+1);
    for (int eh = 0; eh <= mOrderEH; eh++)
    {
        mCoefficientsEH[eh] = std::vector<Precision::Float>(coefficients.size(), 0.0);
        
        for (int mm = 0; mm < coefficients.size(); mm++)
            mCoefficientsEH[eh][mm] = coefficients[mm].denominator()[eh];
    }
    
    // Check that a0 == 1.
    for (int mm = 0; mm < coefficients.size(); mm++)
        assert(fabs(coefficients[mm].denominator()[0] - 1.0) < 1e-4);
}

void AdjointEH::
allocate()
{
    long totalLength = 0;
    for (int nn = 0; nn < mRunlines.size(); nn++)
        totalLength += mRunlines[nn].length();
    
    if (mOrderEH > 1)
    {
        mBufferEH.resize(mOrderEH-1);
        for (int eh = 0; eh < mBufferEH.size(); eh++)
            mBufferEH[eh] = std::vector<Precision::Float>(totalLength, 0.0);
    }
    
//    LOG << "allocated.\n";
}

void AdjointEH::
setPointers(GridFields & currentGridFields,
    std::map<int, Pointer<GridFields> > & allGridFields)
{
    mHeadEHUpdate = currentGridFields.field(mField);
    mHeadDBNeighbors_j = currentGridFields.field(mFieldNeighborDBj);
    mHeadDBNeighbors_k = currentGridFields.field(mFieldNeighborDBk);
}

unsigned long AdjointEH::
bytes() const
{
    long totalBytes = 0;
    totalBytes += mRunlines.size()*sizeof(RunlineAdjointEH);
    
    totalBytes += sizeof(std::vector<Precision::Float>)*(
        mCoefficientsEH.size() + mBufferEH.size());
    
    for (int nn = 0; nn < mCoefficientsEH.size(); nn++)
        totalBytes += sizeof(Precision::Float)*mCoefficientsEH[nn].size();
    for (int nn = 0; nn < mBufferEH.size(); nn++)
        totalBytes += sizeof(Precision::Float)*mBufferEH[nn].size();
    
    return totalBytes;
}

void AdjointEH::
apply(long timestep, Precision::Float dt)
{
    for (int rr = 0; rr < mRunlines.size(); rr++)
    {
        Precision::Float* eh = mHeadEHUpdate + mRunlines[rr].eh;
        const Precision::Float* gjLow = mHeadDBNeighbors_j + mRunlines[rr].gjLow;
        const Precision::Float* gjHigh = mHeadDBNeighbors_j + mRunlines[rr].gjHigh;
        const Precision::Float* gkLow = mHeadDBNeighbors_k + mRunlines[rr].gkLow;
        const Precision::Float* gkHigh = mHeadDBNeighbors_k + mRunlines[rr].gkHigh;
        
        const Precision::Float* constEH[10];
        Precision::Float* oldEH[9];
        for (int nn = 0; nn <= mOrderEH; nn++)
            constEH[nn] = &(mCoefficientsEH[nn].at(mRunlines[rr].coeff));
        for (int nn = 0; nn < mOrderEH-1; nn++)
        {
            oldEH[nn] = &(bufferEH(nn).at(mRunlines[rr].aux));
            assert(mRunlines[rr].aux + mRunlines[rr].length() - 1 <
                bufferEH(nn).size());
        }
        
        int stride = mRunlines[rr].coeffStride;
        
        Precision::Float c1 = m_c1, c2 = m_c2;
        for (int ll = 0; ll < mRunlines[rr].length(); ll++)
        {
            Precision::Float prevEH = *eh;
            
            Precision::Float sumE = 0;
            if (mOrderEH > 0)
                sumE -= (*constEH[1])*prevEH;
            for (int nn = 2; nn <= mOrderEH; nn++)
                sumE -= (*constEH[nn])*(*oldEH[nn-2]);
            
            // c1 is 1/dj, and c2 is 1/dk.
            *eh = c1*(*gkHigh - *gkLow) + c2*(*gjHigh - *gjLow) + sumE;
            
            if (mOrderEH > 1)
                *(oldEH[mOrderEH-2]) = prevEH;
            for (int nn = 0; nn <= mOrderEH; nn++)
                constEH[nn] += stride;
            for (int nn = 0; nn < mOrderEH-1; nn++)
                oldEH[nn]++;
            
            eh++;
            gjLow++;
            gjHigh++;
            gkLow++;
            gkHigh++;
        }
    }
    mCircularBufferLag++;
}

void AdjointEH::
printRunlines(std::ostream & str) const
{
    for (int nn = 0; nn < mRunlines.size(); nn++)
        str << mRunlines[nn] << "\n";
}





AdjointEH_JM::
AdjointEH_JM() :
    mCircularBufferLag(0)
{
}

AdjointEH_JM::
AdjointEH_JM(Field eh, Field dbj, Field dbk, Field jm,
    int ehOrder,
    Precision::Float curlCoeff1, Precision::Float curlCoeff2, Precision::Float currentCoefficient,
    const std::vector<RunlineAdjointEH_JM> & runlines,
    const std::vector<Precision::RationalFunction> & coefficients,
    Precision::Float eps0_or_mu0) :
    mCircularBufferLag(0),
    mField(eh), mFieldNeighborDBj(dbj), mFieldNeighborDBk(dbk),
    mFieldCurrent(jm),
    mRunlines(runlines),
    mOrderEH(ehOrder),
    m_c1(curlCoeff1/eps0_or_mu0),
    m_c2(curlCoeff2/eps0_or_mu0),
    m_c3(currentCoefficient/eps0_or_mu0)
//    m_eps0_or_mu0(eps0_or_mu0)
{
    mCoefficientsEH.resize(mOrderEH+1);
    for (int eh = 0; eh <= mOrderEH; eh++)
    {
        mCoefficientsEH[eh] = std::vector<Precision::Float>(coefficients.size(), 0.0);
        
        for (int mm = 0; mm < coefficients.size(); mm++)
            mCoefficientsEH[eh][mm] = coefficients[mm].denominator()[eh];
    }
    
    // Check that a0 == 1.
    for (int mm = 0; mm < coefficients.size(); mm++)
        assert(fabs(coefficients[mm].denominator()[0] - 1.0) < 1e-4);
}

void AdjointEH_JM::
allocate()
{
    long totalLength = 0;
    for (int nn = 0; nn < mRunlines.size(); nn++)
        totalLength += mRunlines[nn].length();
    
    if (mOrderEH > 1)
    {
        mBufferEH.resize(mOrderEH-1);
        for (int eh = 0; eh < mBufferEH.size(); eh++)
            mBufferEH[eh] = std::vector<Precision::Float>(totalLength, 0.0);
    }
    
//    LOG << "allocated.\n";
}

void AdjointEH_JM::
setPointers(GridFields & currentGridFields,
    std::map<int, Pointer<GridFields> > & allGridFields)
{
    mHeadEHUpdate = currentGridFields.field(mField);
    mHeadDBNeighbors_j = currentGridFields.field(mFieldNeighborDBj);
    mHeadDBNeighbors_k = currentGridFields.field(mFieldNeighborDBk);
    mHeadCurrent = currentGridFields.field(mFieldCurrent);
}

unsigned long AdjointEH_JM::
bytes() const
{
    long totalBytes = 0;
    totalBytes += mRunlines.size()*sizeof(RunlineAdjointEH);
    
    totalBytes += sizeof(std::vector<Precision::Float>)*(
        mCoefficientsEH.size() + mBufferEH.size());
    
    for (int nn = 0; nn < mCoefficientsEH.size(); nn++)
        totalBytes += sizeof(Precision::Float)*mCoefficientsEH[nn].size();
    for (int nn = 0; nn < mBufferEH.size(); nn++)
        totalBytes += sizeof(Precision::Float)*mBufferEH[nn].size();
    
    return totalBytes;
}

void AdjointEH_JM::
apply(long timestep, Precision::Float dt)
{
//    std::cerr << "Lag " << mCircularBufferLag << "\n";
    for (int rr = 0; rr < mRunlines.size(); rr++)
    {
        Precision::Float* eh = mHeadEHUpdate + mRunlines[rr].eh;
        const Precision::Float* gjLow = mHeadDBNeighbors_j + mRunlines[rr].gjLow;
        const Precision::Float* gjHigh = mHeadDBNeighbors_j + mRunlines[rr].gjHigh;
        const Precision::Float* gkLow = mHeadDBNeighbors_k + mRunlines[rr].gkLow;
        const Precision::Float* gkHigh = mHeadDBNeighbors_k + mRunlines[rr].gkHigh;
        const Precision::Float* current = mHeadCurrent + mRunlines[rr].current;
        
        const Precision::Float* constEH[10];
        Precision::Float* oldEH[9];
        for (int lag = 0; lag <= mOrderEH; lag++)
            constEH[lag] = &(mCoefficientsEH[lag].at(mRunlines[rr].coeff));
        for (int lag = 2; lag <= mOrderEH; lag++)
        {
//            std::cerr << "setting up buffer for lag " << lag << "\n";
            oldEH[lag-2] = &(bufferEH(lag).at(mRunlines[rr].aux));
            assert(mRunlines[rr].aux + mRunlines[rr].length() - 1 <
                bufferEH(lag).size());
        }
        
//        for (int nn = 0; nn < mOrderEH-1; nn++)
//        {
//            oldEH[nn] = &(bufferEH(nn).at(mRunlines[rr].aux));
//            assert(mRunlines[rr].aux + mRunlines[rr].length() - 1 <
//                bufferEH(nn).size());
//        }
        
        int stride = mRunlines[rr].coeffStride;
        
        Precision::Float c1 = m_c1, c2 = m_c2, currCoeff = m_c3;
        for (int ll = 0; ll < mRunlines[rr].length(); ll++)
        {
            Precision::Float prevEH = *eh;

//            if (mOrderEH > 1 && *eh != 0)
//            {
//                std::cerr << "ll = " << ll << " eh " << *eh;
//                for (int lag = 2; lag <= mOrderEH; lag++)
//                {
//                    std::cerr << " " << *oldEH[lag-2];
//                }
//                std::cerr << "\n";
//            }
            
            Precision::Float sumE = 0;
            if (mOrderEH > 0)
                sumE -= (*constEH[1])*prevEH;
            for (int lag = 2; lag <= mOrderEH; lag++)
                sumE -= (*constEH[lag])*(*oldEH[lag-2]);
            
            *eh = c1*(*gkHigh - *gkLow) + c2*(*gjHigh - *gjLow) + sumE
                + currCoeff* (*current);
            
            if (mOrderEH > 1)
                *(oldEH[mOrderEH-2]) = prevEH;
            for (int lag = 0; lag <= mOrderEH; lag++)
                constEH[lag] += stride;
            for (int nn = 0; nn < mOrderEH-1; nn++)
                oldEH[nn]++;
            
            eh++;
            gjLow++;
            gjHigh++;
            gkLow++;
            gkHigh++;
            current++;
        }
    }
    mCircularBufferLag++;
}

void AdjointEH_JM::
printRunlines(std::ostream & str) const
{
    for (int nn = 0; nn < mRunlines.size(); nn++)
        str << mRunlines[nn] << "\n";
}

