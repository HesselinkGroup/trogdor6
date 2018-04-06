/*
 *  ForwardNondispersiveEH_2d.h
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 3/30/18.
 *  Copyright 2018 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef FORWARDNONDISPERSIVEEH_2D_H
#define FORWARDNONDISPERSIVEEH_2D_H

#include <iostream>
#include <map>
#include "Pointer.h"
#include "../FieldEnum.h"
#include "../GridFields.h"
#include "../Operation.h"

struct RunlineNondispersiveEH_2d
{
    RunlineNondispersiveEH_2d() {}
    RunlineNondispersiveEH_2d(long f, long gj1, long gj2,
        long inCoeff, long inCoeffStride, long length) :
        fi(f),
        gjLow(gj1),
        gjHigh(gj2),
        coeff(inCoeff),
        coeffStride(inCoeffStride),
        mLength(length)
    {
        assert(length > 0);
    }
    long length() const { return mLength; }
    long fi;
    long gjLow, gjHigh;
    long coeff;
    long coeffStride;
    long mLength;
};
inline std::ostream & operator<<(std::ostream & str,
    const RunlineNondispersiveEH & rl)
{
    str << "[eh " << rl.fi << " hej " << rl.gjLow << " " << rl.gjHigh
        << " coeff " << rl.coeff << "x" << rl.coeffStride
        << " length " << rl.length() << "]";
    return str;
}

struct RunlineNondispersiveEH_JM_2d : public RunlineNondispersiveEH_2d
{
    RunlineNondispersiveEH_JM_2d() {}
    RunlineNondispersiveEH_JM_2d(long f, long gj1, long gj2,
        long inCoeff, long inCoeffStride,
        long curr, long length) :
        RunlineNondispersiveEH(f, gj1, gj2, inCoeff, inCoeffStride, length),
        current(curr)
    {
    }
    
    long current;
};
inline std::ostream & operator<<(std::ostream & str,
    const RunlineNondispersiveEH_JM_2d & rl)
{
    str << "[eh " << rl.fi << " hej " << rl.gjLow << " " << rl.gjHigh
        << " coeff " << rl.coeff << "x" << rl.coeffStride
        << " current " << rl.current
        << " length " << rl.length() << "]";
    return str;
}

class ForwardNondispersiveEH_2d : public Operation
{
public:
    ForwardNondispersiveEH_2d()
    {
    }
    
    ForwardNondispersiveEH_2d(Field eh, Field hej,
        Precision::Float curlCoeff,
        const std::vector<RunlineNondispersiveEH> & runlines,
        const std::vector<Precision::RationalFunction> & coefficients,
        Precision::Float eps0_or_mu0) :
        mFieldEH(eh), mFieldHEj(hej),
        mRunlines(runlines),
        m_curlCoeff(curlCoeff),
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
        
        name("EH nondispersive 2d");
    }
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields)
    {
        mHeadEHUpdate = currentGridFields.field(mFieldEH);
        mHeadHENeighbors_j = currentGridFields.field(mFieldHEj);
    }
    
    unsigned long bytes() const
    {
        return mRunlines.size() * sizeof(RunlineNondispersiveEH_2d) +
            mCoefficients.size() * sizeof(Precision::Float);
    }
    
    void apply(long timestep, Precision::Float dt)
    {
        for (int rr = 0; rr < mRunlines.size(); rr++)
        {
            Precision::Float* f = mHeadEHUpdate + mRunlines[rr].fi;
            const Precision::Float* gjLow = mHeadHENeighbors_j + mRunlines[rr].gjLow;
            const Precision::Float* gjHigh = mHeadHENeighbors_j + mRunlines[rr].gjHigh;
            const Precision::Float* coeff = &(mCoefficients.at(mRunlines[rr].coeff));
            
            Precision::Float curlCoeff = m_curlCoeff;
            int stride = mRunlines[rr].coeffStride;
            
            for (int ll = 0; ll < mRunlines[rr].length(); ll++)
            {
                *f = *f + *coeff * curlCoeff * (*gjHigh - *gjLow);
                f++;
                gjLow++;
                gjHigh++;
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
    Precision::Float m_curlCoeff;
    Precision::Float* mHeadEHUpdate;
    Precision::Float* mHeadHENeighbors_j;
};


class ForwardNondispersiveEH_JM_2d : public Operation
{
public:
    ForwardNondispersiveEH_JM_2d()
    {
    }
    
    ForwardNondispersiveEH_JM_2d(Field eh, Field hej, Field jm,
        Precision::Float curlCoeff,
        Precision::Float currentCoeff,
        const std::vector<RunlineNondispersiveEH_JM_2d> & runlines,
        const std::vector<Precision::RationalFunction> & coefficients,
        Precision::Float eps0_or_mu0) :
        mFieldEH(eh), mFieldHEj(hej), mFieldCurrent(jm),
        mRunlines(runlines),
        m_curlCoeff(curlCoeff),
        m_currentCoeff(currentCoeff)
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
        name("EH+JM nondispersive 2d");
    }
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields)
    {
        mHeadEHUpdate = currentGridFields.field(mFieldEH);
        mHeadHENeighbors_j = currentGridFields.field(mFieldHEj);
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
            const Precision::Float* coeff = &(mCoefficients.at(mRunlines[rr].coeff));
            const Precision::Float* current = mHeadCurrent + mRunlines[rr].current;
            
            Precision::Float curlCoeff = m_curlCoeff, currentCoeff = m_currentCoeff;
            int stride = mRunlines[rr].coeffStride;
            
            for (int ll = 0; ll < mRunlines[rr].length(); ll++)
            {
                *f = *f + *coeff*(curlCoeff*(*gjHigh - *gjLow) + currentCoeff*(*current));
                f++;
                gjLow++;
                gjHigh++;
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
    Field mFieldEH, mFieldHEj, mFieldCurrent;
    std::vector<RunlineNondispersiveEH_JM_2d> mRunlines;
    std::vector<Precision::Float> mCoefficients;
    Precision::Float m_curlCoeff, m_currentCoeff
    Precision::Float* mHeadEHUpdate;
    Precision::Float* mHeadHENeighbors_j;
    Precision::Float* mHeadCurrent;
};





#endif
