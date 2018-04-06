/*
 *  ForwardEH_j.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 7/13/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "ForwardEH_j.h"

#include "Log.h"
#include <iostream>
#include <cmath>

using namespace std;


std::ostream & operator<<(std::ostream & str, const RunlineSumEH_j & rl)
{
    str << "[ehi " << rl.ehi << " ehii " << rl.ehii
        << " ehij " << rl.ehij0 << "," << rl.ehij1 << "," << rl.ehij2 << "," << rl.ehij3
        << " ehik " << rl.ehik0 << "," << rl.ehik1 << "," << rl.ehik2 << "," << rl.ehik3
        << " length " << rl.length() << "]";
    return str;
}

std::ostream & operator<<(std::ostream & str, const RunlineAnisotropicEH_j & rl)
{
    str << "[eh " << rl.eh << " db " << rl.db << " aux " << rl.aux
        << " coeff " << rl.coeff << "x" << rl.coeffStride
        << " length " << rl.length() << "]";
    return str;
}


#pragma mark *** SumEH ***


SumEH_j::
SumEH_j()
{
}

SumEH_j::
SumEH_j(Field ehi, Field ehii, Field ehij, Field ehik,
    const vector<RunlineSumEH_j> & runlines) :
    mField_i(ehi), mField_ii(ehii), mField_ij(ehij), mField_ik(ehik),
    mRunlines(runlines)
{
}

void SumEH_j::
setPointers(GridFields & currentGridFields,
    std::map<int, Pointer<GridFields> > & allGridFields)
{
    mHeadEHi = currentGridFields.field(mField_i);
    mHeadEHii = currentGridFields.field(mField_ii);
    mHeadEHij = currentGridFields.field(mField_ij);
    mHeadEHik = currentGridFields.field(mField_ik);
}


unsigned long SumEH_j::
bytes() const
{
    return mRunlines.size() * sizeof(RunlineSumEH_j);
}

void SumEH_j::
apply(long timestep, Precision::Float dt)
{
    for (int rr = 0; rr < mRunlines.size(); rr++)
    {
        Precision::Float* ehi = mHeadEHi + mRunlines[rr].ehi;
        Precision::Float* ehii = mHeadEHii + mRunlines[rr].ehii;
        Precision::Float* ehij0 = mHeadEHij + mRunlines[rr].ehij0;
        Precision::Float* ehik0 = mHeadEHik + mRunlines[rr].ehik0;
        Precision::Float* ehij1 = mHeadEHij + mRunlines[rr].ehij1;
        Precision::Float* ehik1 = mHeadEHik + mRunlines[rr].ehik1;
        Precision::Float* ehij2 = mHeadEHij + mRunlines[rr].ehij2;
        Precision::Float* ehik2 = mHeadEHik + mRunlines[rr].ehik2;
        Precision::Float* ehij3 = mHeadEHij + mRunlines[rr].ehij3;
        Precision::Float* ehik3 = mHeadEHik + mRunlines[rr].ehik3;
        
        for (int ll = 0; ll < mRunlines[rr].length(); ll++)
        {
            ehi[ll] = ehii[ll] + 0.25*(
                ehij0[ll] + ehij1[ll] + ehij2[ll] + ehij3[ll] +
                ehik0[ll] + ehik1[ll] + ehik2[ll] + ehik3[ll]);
        }
    }
}

void SumEH_j::
printRunlines(std::ostream & str) const
{
    for (int nn = 0; nn < mRunlines.size(); nn++)
        str << mRunlines[nn] << "\n";
}


#pragma mark *** UpdateAnisotropicEH ***


UpdateAnisotropicEH_j::
UpdateAnisotropicEH_j() :
    mCircularBufferLag(0)
{
}

// I take the coefficients to be inverse permittivity or permeability.
// The denominator is associated with E or H.
// The numerator is associated with D or B.
// Zeroth order in numerator and denominator means a dielectric.
UpdateAnisotropicEH_j::
UpdateAnisotropicEH_j(Field eh, Field db,
    int ehOrder, int dbOrder,
    const std::vector<RunlineAnisotropicEH_j> & runlines,
    const std::vector<Precision::RationalFunction> & coefficients,
    Precision::Float eps0_or_mu0) :
    mCircularBufferLag(0),
    mField(eh), mAuxField(db),
    mRunlines(runlines),
    mOrderEH(ehOrder),
    mOrderDB(dbOrder)
//    m_eps0_or_mu0(eps0_or_mu0)
{
    mCoefficientsDB.resize(mOrderDB+1);
    for (int db = 0; db <= mOrderDB; db++)
    {
        mCoefficientsDB.at(db) = std::vector<Precision::Float>(coefficients.size(), 0.0);
        
        // Numerator of 1/eps is [b0/a0, b1/a0, b2/a0, ...]
        // The update constants are [b0/a0/eps0, b1/a0/eps0, ...]
        for (int mm = 0; mm < coefficients.size(); mm++)
        {
            mCoefficientsDB.at(db).at(mm) = coefficients[mm].numerator()[db] /
                eps0_or_mu0;
        }
    }
    
    mCoefficientsEH.resize(mOrderEH+1);
    for (int eh = 0; eh <= mOrderEH; eh++)
    {
        mCoefficientsEH.at(eh) = std::vector<Precision::Float>(coefficients.size(), 0.0);
        
        // Denominator of 1/eps is [1, a1/a0, a2/a0, ...]
        // The update constants are the same.
        for (int mm = 0; mm < coefficients.size(); mm++)
            mCoefficientsEH.at(eh).at(mm) = coefficients[mm].denominator()[eh];
    }
    
    // Check that a0 == 1.
    for (int mm = 0; mm < coefficients.size(); mm++)
        assert(fabs(coefficients[mm].denominator()[0] - 1.0) < 1e-4);
}

void UpdateAnisotropicEH_j::
allocate()
{
    long totalLength = 0;
    for (int nn = 0; nn < mRunlines.size(); nn++)
        totalLength += mRunlines[nn].length();
    
    mBufferDB.resize(mOrderDB);
    for (int db = 0; db < mBufferDB.size(); db++)
        mBufferDB[db] = std::vector<Precision::Float>(totalLength, 0.0);
    
    if (mOrderEH > 1)
    {
        mBufferEH.resize(mOrderEH-1);
        for (int eh = 0; eh < mBufferEH.size(); eh++)
            mBufferEH[eh] = std::vector<Precision::Float>(totalLength, 0.0);
    }
    
//    LOG << "allocated.\n";
}

void UpdateAnisotropicEH_j::
setPointers(GridFields & currentGridFields,
    std::map<int, Pointer<GridFields> > & allGridFields)
{
    mHeadEH = currentGridFields.field(mField);
    mHeadDB = currentGridFields.field(mAuxField);
}

unsigned long UpdateAnisotropicEH_j::
bytes() const
{
    long totalBytes = 0;
    totalBytes += mRunlines.size()*sizeof(RunlineAnisotropicEH_j);
    
    totalBytes += sizeof(std::vector<Precision::Float>)*(
        mCoefficientsDB.size() + mCoefficientsEH.size() +
        mBufferDB.size() + mBufferEH.size());
    
    for (int nn = 0; nn < mCoefficientsDB.size(); nn++)
        totalBytes += sizeof(Precision::Float)*mCoefficientsDB[nn].size();
    for (int nn = 0; nn < mCoefficientsEH.size(); nn++)
        totalBytes += sizeof(Precision::Float)*mCoefficientsEH[nn].size();
    for (int nn = 0; nn < mBufferDB.size(); nn++)
        totalBytes += sizeof(Precision::Float)*mBufferDB[nn].size();
    for (int nn = 0; nn < mBufferEH.size(); nn++)
        totalBytes += sizeof(Precision::Float)*mBufferEH[nn].size();
    
    return totalBytes;
}

void UpdateAnisotropicEH_j::
apply(long timestep, Precision::Float dt)
{
    for (int rr = 0; rr < mRunlines.size(); rr++)
    {
        Precision::Float* db = mHeadDB + mRunlines[rr].db;
        Precision::Float* eh = mHeadEH + mRunlines[rr].eh;
        
        const Precision::Float* constDB[10];
        const Precision::Float* constEH[10];
        Precision::Float* oldEH[9]; // oldEH[0] is the array of EH from two timesteps ago
        Precision::Float* oldDB[9]; // oldDB[0] is the array of DB at the previous timestep
        
        // remember: lag 0 is always available on the grid.
        for (int lag = 0; lag <= mOrderDB; lag++)
            constDB[lag] = &(mCoefficientsDB[lag].at(mRunlines[rr].coeff));
        for (int lag = 1; lag <= mOrderDB; lag++)
        {
            oldDB[lag-1] = &(bufferDB(lag).at(mRunlines[rr].aux));
            assert(mRunlines[rr].aux + mRunlines[rr].length() - 1 <
                bufferDB(lag).size());
        }
        
        for (int lag = 0; lag <= mOrderEH; lag++)
            constEH[lag] = &(mCoefficientsEH[lag].at(mRunlines[rr].coeff));
        for (int lag = 2; lag <= mOrderEH; lag++)
        {
            oldEH[lag-2] = &(bufferEH(lag).at(mRunlines[rr].aux));
            assert(mRunlines[rr].aux + mRunlines[rr].length() - 1 <
                bufferEH(lag).size());
//            cerr << "oldEH[" << nn << "] = " << oldEH[nn] << ", use "
//                << mRunlines[rr].aux << "\n";
        }
        
        //std::cerr << "Current runline " << mRunlines[rr] << "\n";
//            Precision::Float* constDB = &(mCoefficientsDB[0].at(mRunlines[rr].coeff));
//            Precision::Float* constEH = &(mCoefficientsEH[0].at(mRunlines[rr].coeff));
        int stride = mRunlines[rr].coeffStride;
        
        for (int ll = 0; ll < mRunlines[rr].length(); ll++)
        {
            Precision::Float prevEH = *eh;
            
            // constDB are the "b" and "d" coefficients, the denominators
            // of the permittivity or permeability.
            
            Precision::Float sumD = (*constDB[0]) * (*db);
            for (int lag = 1; lag <= mOrderDB; lag++)
                sumD += (*constDB[lag])*(*oldDB[lag-1]);
            
            Precision::Float sumE = 0;
            if (mOrderEH > 0)
                sumE -= (*constEH[1])*prevEH;
            for (int lag = 2; lag <= mOrderEH; lag++)
                sumE -= (*constEH[lag])*(*oldEH[lag-2]);
            
            //*eh = ( sumD + sumE )/(*constEH[0]);
            //cerr << "*constEH[0] = " << *constEH[0] << "\n";
            *eh = sumD + sumE; // CAREFUL.  I assume *constEH[0] == 1.
            
            if (mOrderEH > 1)
                *(oldEH[mOrderEH-2]) = prevEH;
            if (mOrderDB > 0)
                *(oldDB[mOrderDB-1]) = *db;
            
            eh++;
            db++;
            
            for (int lag = 0; lag <= mOrderDB; lag++)
                constDB[lag] += stride;
            for (int nn = 0; nn < mOrderDB; nn++)
                oldDB[nn]++;
            for (int lag = 0; lag <= mOrderEH; lag++)
                constEH[lag] += stride;
            for (int nn = 0; nn < mOrderEH-1; nn++)
                oldEH[nn]++;
        }
    }
    mCircularBufferLag++;
}

void UpdateAnisotropicEH_j::
printRunlines(std::ostream & str) const
{
    for (int nn = 0; nn < mRunlines.size(); nn++)
        str << mRunlines[nn] << "\n";
}



