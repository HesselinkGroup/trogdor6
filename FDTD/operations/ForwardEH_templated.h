/*
 *  ForwardEH_templated.h
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 2/11/12.
 *  Copyright 2012 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef FORWARDEH_TEMPLATED_H
#define FORWARDEH_TEMPLATED_H

#include "ForwardEH.h" // for RunlineEH

template<int ORDER>
class ForwardEHT : public Operation
{
public:
    ForwardEHT();
    
    ForwardEHT(Field eh, Field db, const std::vector<RunlineEH> & runlines,
        const std::vector<Precision::RationalFunction> & coefficients,
        Precision::Float eps0_or_mu0);
    
    void allocate();
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields);
    
    unsigned long bytes() const;
    virtual void apply(long timestep, Precision::Float dt);
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
    Precision::Float m_eps0_or_mu0; // either eps0 or mu0
};


class FwdEHFactory
{
public:
    static Operation* newForwardEH(Field eh, Field db,
        int order,
        const std::vector<RunlineEH> & runlines,
        const std::vector<Precision::RationalFunction> & coefficients,
        Precision::Float eps0_or_mu0)
    {
        if (order == 0)
            return (Operation*)(new ForwardEHT<0>(eh, db, runlines,
                coefficients, eps0_or_mu0));
        else if (order == 1)
            return (Operation*)(new ForwardEHT<1>(eh, db, runlines,
                coefficients, eps0_or_mu0));
        else if (order == 2)
            return (Operation*)(new ForwardEHT<2>(eh, db, runlines,
                coefficients, eps0_or_mu0));
        else if (order == 3)
            return (Operation*)(new ForwardEHT<3>(eh, db, runlines,
                coefficients, eps0_or_mu0));
        else if (order == 4)
            return (Operation*)(new ForwardEHT<4>(eh, db, runlines,
                coefficients, eps0_or_mu0));
        else if (order == 5)
            return (Operation*)(new ForwardEHT<5>(eh, db, runlines,
                coefficients, eps0_or_mu0));
        else if (order == 6)
            return (Operation*)(new ForwardEHT<6>(eh, db, runlines,
                coefficients, eps0_or_mu0));
        else if (order == 7)
            return (Operation*)(new ForwardEHT<7>(eh, db, runlines,
                coefficients, eps0_or_mu0));
        else if (order == 8)
            return (Operation*)(new ForwardEHT<8>(eh, db, runlines,
                coefficients, eps0_or_mu0));
        else if (order == 9)
            return (Operation*)(new ForwardEHT<9>(eh, db, runlines,
                coefficients, eps0_or_mu0));
        else
        {
            throw(std::logic_error("Order must be between 0 and 9 inclusive"));
        }
    }
};




#include "ForwardEH_templated-inl.h"

#endif
