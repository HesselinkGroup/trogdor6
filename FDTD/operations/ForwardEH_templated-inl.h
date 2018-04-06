/*
 *  ForwardEH_templated-inl.h
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 2/11/12.
 *  Copyright 2012 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef FORWARDEH_TEMPLATED_H

#include <sstream>

template<int ORDER>
ForwardEHT<ORDER>::
ForwardEHT() :
    mCircularBufferLag(0)
{
}

// I take the coefficients to be inverse permittivity or permeability.
// The denominator is associated with E or H.
// The numerator is associated with D or B.
// Zeroth order in numerator and denominator means a dielectric.
template<int ORDER>
ForwardEHT<ORDER>::
ForwardEHT(Field eh, Field db,
    const std::vector<RunlineEH> & runlines,
    const std::vector<Precision::RationalFunction> & coefficients,
    Precision::Float eps0_or_mu0) :
    mCircularBufferLag(0),
    mField(eh), mAuxField(db),
    mRunlines(runlines),
    m_eps0_or_mu0(eps0_or_mu0)
{
    mCoefficientsDB.resize(ORDER+1);
    for (int db = 0; db <= ORDER; db++)
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
    
    mCoefficientsEH.resize(ORDER+1);
    for (int eh = 0; eh <= ORDER; eh++)
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
        
    std::ostringstream str;
    str << "EH<T> O(" << ORDER << ")";
    name(str.str());
}

template<int ORDER>
void ForwardEHT<ORDER>::
allocate()
{
    long totalLength = 0;
    for (int nn = 0; nn < mRunlines.size(); nn++)
        totalLength += mRunlines[nn].length();
    
    mBufferDB.resize(ORDER);
    for (int db = 0; db < mBufferDB.size(); db++)
        mBufferDB[db] = std::vector<Precision::Float>(totalLength, 0.0);
    
    if (ORDER > 1)
    {
        mBufferEH.resize(ORDER-1);
        for (int eh = 0; eh < mBufferEH.size(); eh++)
            mBufferEH[eh] = std::vector<Precision::Float>(totalLength, 0.0);
    }
    
//    LOG << "allocated.\n";
}

template<int ORDER>
void ForwardEHT<ORDER>::
setPointers(GridFields & currentGridFields,
    std::map<int, Pointer<GridFields> > & allGridFields)
{
    mHeadEH = currentGridFields.field(mField);
    mHeadDB = currentGridFields.field(mAuxField);
}

template<int ORDER>
unsigned long ForwardEHT<ORDER>::
bytes() const
{
    long totalBytes = 0;
    totalBytes += mRunlines.size()*sizeof(RunlineEH);
    
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

template<int ORDER>
void ForwardEHT<ORDER>::
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
        for (int lag = 0; lag <= ORDER; lag++)
            constDB[lag] = &(mCoefficientsDB[lag].at(mRunlines[rr].coeff));
        for (int lag = 1; lag <= ORDER; lag++)
        {
            oldDB[lag-1] = &(bufferDB(lag).at(mRunlines[rr].aux));
            assert(mRunlines[rr].aux + mRunlines[rr].length() - 1 <
                bufferDB(lag).size());
        }
        
        for (int lag = 0; lag <= ORDER; lag++)
            constEH[lag] = &(mCoefficientsEH[lag].at(mRunlines[rr].coeff));
        for (int lag = 2; lag <= ORDER; lag++)
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
            
            Precision::Float sumD = (*constDB[0])*(*db);
            for (int lag = 1; lag <= ORDER; lag++)
                sumD += (*constDB[lag])*(*oldDB[lag-1]);
            
            Precision::Float sumE = 0;
            if (ORDER > 0)
                sumE -= (*constEH[1])*prevEH;
            for (int lag = 2; lag <= ORDER; lag++)
                sumE -= (*constEH[lag])*(*oldEH[lag-2]);
            
            //*eh = ( sumD + sumE )/(*constEH[0]);
            //cerr << "*constEH[0] = " << *constEH[0] << "\n";
            *eh = sumD + sumE; // CAREFUL.  I assume *constEH[0] == 1.
            
            if (ORDER > 1)
                *(oldEH[ORDER-2]) = prevEH;
            if (ORDER > 0)
                *(oldDB[ORDER-1]) = *db;
            
            eh++;
            db++;
            
            for (int lag = 0; lag <= ORDER; lag++)
                constDB[lag] += stride;
            for (int nn = 0; nn < ORDER; nn++)
                oldDB[nn]++;
            for (int lag = 0; lag <= ORDER; lag++)
                constEH[lag] += stride;
            for (int nn = 0; nn < ORDER-1; nn++)
                oldEH[nn]++;
        }
    }
    mCircularBufferLag++;
}

template<int ORDER>
void ForwardEHT<ORDER>::
printRunlines(std::ostream & str) const
{
    for (int nn = 0; nn < mRunlines.size(); nn++)
        str << mRunlines[nn] << "\n";
}




#endif
