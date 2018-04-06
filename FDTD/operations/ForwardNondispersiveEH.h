/*
 *  ForwardNondispersiveEH.h
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 3/4/12.
 *  Copyright 2012 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef FORWARDNONDISPERSIVEEH_H
#define FORWARDNONDISPERSIVEEH_H

#include <iostream>
#include <map>
#include "Pointer.h"
#include "../FieldEnum.h"
#include "../GridFields.h"
#include "../Operation.h"

struct RunlineNondispersiveEH
{
    RunlineNondispersiveEH() {}
    RunlineNondispersiveEH(long f, long gj1, long gj2, long gk1, long gk2,
        long inCoeff, long inCoeffStride, long length) :
        fi(f),
        gjLow(gj1),
        gjHigh(gj2),
        gkLow(gk1),
        gkHigh(gk2),
        coeff(inCoeff),
        coeffStride(inCoeffStride),
        mLength(length)
    {
        assert(length > 0);
    }
    long length() const { return mLength; }
    long fi;
    long gjLow, gjHigh, gkLow, gkHigh;
    long coeff;
    long coeffStride;
    long mLength;
};
inline std::ostream & operator<<(std::ostream & str,
    const RunlineNondispersiveEH & rl)
{
    str << "[eh " << rl.fi << " hej " << rl.gjLow << " " << rl.gjHigh
        << " hek " << rl.gkLow << " " << rl.gkHigh
        << " coeff " << rl.coeff << "x" << rl.coeffStride
        << " length " << rl.length() << "]";
    return str;
}

struct RunlineNondispersiveEH_JM : public RunlineNondispersiveEH
{
    RunlineNondispersiveEH_JM() {}
    RunlineNondispersiveEH_JM(long f, long gj1, long gj2, long gk1, long gk2,
        long inCoeff, long inCoeffStride,
        long curr, long length) :
        RunlineNondispersiveEH(f, gj1, gj2, gk1, gk2, inCoeff, inCoeffStride,
            length),
        current(curr)
    {
    }
    
    long current;
};
inline std::ostream & operator<<(std::ostream & str,
    const RunlineNondispersiveEH_JM & rl)
{
    str << "[eh " << rl.fi << " hej " << rl.gjLow << " " << rl.gjHigh
        << " hek " << rl.gkLow << " " << rl.gkHigh
        << " coeff " << rl.coeff << "x" << rl.coeffStride
        << " current " << rl.current
        << " length " << rl.length() << "]";
    return str;
}

class ForwardNondispersiveEH : public Operation
{
public:
    ForwardNondispersiveEH()
    {
    }
    
    ForwardNondispersiveEH(Field eh, Field hej, Field hek,
        Precision::Float curlCoeff1, Precision::Float curlCoeff2,
        const std::vector<RunlineNondispersiveEH> & runlines,
        const std::vector<Precision::RationalFunction> & coefficients,
        Precision::Float eps0_or_mu0) :
        mFieldEH(eh), mFieldHEj(hej), mFieldHEk(hek),
        mRunlines(runlines),
        m_c1(curlCoeff1),
        m_c2(curlCoeff2)
    {
        mCoefficients.resize(coefficients.size());
        for (int mm = 0; mm < coefficients.size(); mm++)
        {
            assert(coefficients[mm].numerator().order() == 0);
            assert(coefficients[mm].denominator().order() == 0);
            mCoefficients[mm] = coefficients[mm].numerator()[0] /
                coefficients[mm].denominator()[0] /
                eps0_or_mu0;
        }
        
        name("EH nondispersive");
    }
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields)
    {
        mHeadEHUpdate = currentGridFields.field(mFieldEH);
        mHeadHENeighbors_j = currentGridFields.field(mFieldHEj);
        mHeadHENeighbors_k = currentGridFields.field(mFieldHEk);
    }
    
    unsigned long bytes() const
    {
        return mRunlines.size() * sizeof(RunlineNondispersiveEH) +
            mCoefficients.size() * sizeof(Precision::Float);
    }
    
    void apply(long timestep, Precision::Float dt)
    {
        for (int rr = 0; rr < mRunlines.size(); rr++)
        {
            Precision::Float* f = mHeadEHUpdate + mRunlines[rr].fi;
            const Precision::Float* gjLow = mHeadHENeighbors_j + mRunlines[rr].gjLow;
            const Precision::Float* gjHigh = mHeadHENeighbors_j + mRunlines[rr].gjHigh;
            const Precision::Float* gkLow = mHeadHENeighbors_k + mRunlines[rr].gkLow;
            const Precision::Float* gkHigh = mHeadHENeighbors_k + mRunlines[rr].gkHigh;
            const Precision::Float* coeff = &(mCoefficients.at(mRunlines[rr].coeff));
            
            Precision::Float c1 = m_c1, c2 = m_c2;
            int stride = mRunlines[rr].coeffStride;
            
            for (int ll = 0; ll < mRunlines[rr].length(); ll++)
            {
                *f = *f + *coeff * (c1*(*gkHigh - *gkLow) + c2*(*gjHigh - *gjLow));
                f++;
                gjLow++;
                gjHigh++;
                gkLow++;
                gkHigh++;
                coeff += stride;
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
    Field mFieldEH, mFieldHEj, mFieldHEk;
    std::vector<RunlineNondispersiveEH> mRunlines;
    std::vector<Precision::Float> mCoefficients;
    Precision::Float m_c1, m_c2;
    Precision::Float* mHeadEHUpdate;
    Precision::Float* mHeadHENeighbors_j;
    Precision::Float* mHeadHENeighbors_k;
};


class ForwardNondispersiveEH_JM : public Operation
{
public:
    ForwardNondispersiveEH_JM()
    {
    }
    
    ForwardNondispersiveEH_JM(Field eh, Field hej, Field hek, Field jm,
        Precision::Float curlCoeff1, Precision::Float curlCoeff2,
        Precision::Float currentCoeff,
        const std::vector<RunlineNondispersiveEH_JM> & runlines,
        const std::vector<Precision::RationalFunction> & coefficients,
        Precision::Float eps0_or_mu0) :
        mFieldEH(eh), mFieldHEj(hej), mFieldHEk(hek), mFieldCurrent(jm),
        mRunlines(runlines),
        m_c1(curlCoeff1),
        m_c2(curlCoeff2),
        m_c3(currentCoeff)
    {
        mCoefficients.resize(coefficients.size());
        for (int mm = 0; mm < coefficients.size(); mm++)
        {
            assert(coefficients[mm].numerator().order() == 0);
            assert(coefficients[mm].denominator().order() == 0);
            mCoefficients[mm] = coefficients[mm].numerator()[0] /
                coefficients[mm].denominator()[0] /
                eps0_or_mu0;
        }
        name("EH+JM nondispersive");
    }
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields)
    {
        mHeadEHUpdate = currentGridFields.field(mFieldEH);
        mHeadHENeighbors_j = currentGridFields.field(mFieldHEj);
        mHeadHENeighbors_k = currentGridFields.field(mFieldHEk);
        mHeadCurrent = currentGridFields.field(mFieldCurrent);
    }
    
    unsigned long bytes() const
    {
        return mRunlines.size() * sizeof(RunlineNondispersiveEH_JM) +
            mCoefficients.size() * sizeof(Precision::Float);
    }
    
    void apply(long timestep, Precision::Float dt)
    {
        for (int rr = 0; rr < mRunlines.size(); rr++)
        {
            Precision::Float* f = mHeadEHUpdate + mRunlines[rr].fi;
            const Precision::Float* gjLow = mHeadHENeighbors_j + mRunlines[rr].gjLow;
            const Precision::Float* gjHigh = mHeadHENeighbors_j + mRunlines[rr].gjHigh;
            const Precision::Float* gkLow = mHeadHENeighbors_k + mRunlines[rr].gkLow;
            const Precision::Float* gkHigh = mHeadHENeighbors_k + mRunlines[rr].gkHigh;
            const Precision::Float* coeff = &(mCoefficients.at(mRunlines[rr].coeff));
            const Precision::Float* current = mHeadCurrent + mRunlines[rr].current;
            
            Precision::Float c1 = m_c1, c2 = m_c2, c3 = m_c3;
            int stride = mRunlines[rr].coeffStride;
            
            for (int ll = 0; ll < mRunlines[rr].length(); ll++)
            {
                *f = *f + *coeff*(c1*(*gkHigh - *gkLow) + c2*(*gjHigh - *gjLow)
                    + c3*(*current));
                f++;
                gjLow++;
                gjHigh++;
                gkLow++;
                gkHigh++;
                coeff += stride;
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
    Field mFieldEH, mFieldHEj, mFieldHEk, mFieldCurrent;
    std::vector<RunlineNondispersiveEH_JM> mRunlines;
    std::vector<Precision::Float> mCoefficients;
    Precision::Float m_c1, m_c2, m_c3;
    Precision::Float* mHeadEHUpdate;
    Precision::Float* mHeadHENeighbors_j;
    Precision::Float* mHeadHENeighbors_k;
    Precision::Float* mHeadCurrent;
};





#endif
