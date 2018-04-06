/*
 *  UpdateJPML.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/28/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _UPDATEJPML_
#define _UPDATEJPML_

#include "../PhysicalConstants.h"
#include "../FieldEnum.h"
#include "../GridFields.h"
#include "../Operation.h"
#include <vector>


struct RunlineJMPML
{
    RunlineJMPML() {}
    RunlineJMPML(long current, long fieldHigh, long fieldLow, long inAux,
        long inConsts, long inStride, long length) :
        mCurrent(current),
        mFieldHigh(fieldHigh),
        mFieldLow(fieldLow),
        mAux(inAux),
        mConsts(inConsts),
        mStride(inStride),
        mLength(length)
    {
        assert(length > 0);
    }
    
    long length() const { return mLength; }
    
    long mCurrent;
    long mFieldHigh;
    long mFieldLow;
    long mAux;
    long mConsts;
    long mStride;
    long mLength;
};
inline std::ostream & operator<<(std::ostream & str, const RunlineJMPML & rl)
{
    str << "[curr " << rl.mCurrent << " fields " << rl.mFieldLow << ", "
        << rl.mFieldHigh << " aux " << rl.mAux << " consts " << rl.mConsts
            << "x" << rl.mStride << " length "
        << rl.mLength << "]";
    return str;
}


class UpdateJMPML : public Operation
{
public:
    UpdateJMPML()
    {
    }
    
    UpdateJMPML(Field current, Field neighborField,
        Precision::Float leadingSign, const std::vector<RunlineJMPML> & runlines,
        const std::vector<PMLParameters> & params,
        Precision::Float dt, Precision::Float dxyz) :
        mFieldCurrent(current),
        mFieldNeighbor(neighborField),
        mRunlines(runlines),
        m_c_hj(params.size()),
        m_c_phij(params.size()),
        m_c_hphi(params.size()),
        m_c_jphi(params.size()),
        mSign(leadingSign)
    {
        for (int nn = 0; nn < params.size(); nn++)
        {
            Precision::Float alpha = params[nn].alpha();
            Precision::Float kappa = params[nn].kappa();
            Precision::Float sigma = params[nn].sigma();
            
            m_c_hj[nn] = ((2*Constants::eps0 + dt*alpha)/
                (2*Constants::eps0*kappa + dt*(alpha*kappa + sigma)) - 1.0)
                / dxyz; // ALERT: the 1/dx went here!
            
            m_c_phij[nn] = 2*Constants::eps0*kappa / (
                2*Constants::eps0*kappa + dt*(alpha*kappa+sigma));
            
            m_c_hphi[nn] = dt*(alpha - alpha*kappa - sigma) /
                (Constants::eps0*kappa) / dxyz; // ALERT: the 1/dx went here!
            
            m_c_jphi[nn] = dt*(alpha*kappa + sigma) /
                (Constants::eps0*kappa);
            
//            std::cerr << "nn " << nn << " hj " << m_c_hj[nn]
//                << " phij " << m_c_phij[nn]
//                << " hphi " << m_c_hphi[nn]
//                << " jphi " << m_c_jphi[nn]
//                << "\n";
        }
        name("JMPML");
    }
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields)
    {
        mHeadCurrent = currentGridFields.field(mFieldCurrent);
        mHeadNeighborField = currentGridFields.field(mFieldNeighbor);
    }
    
    void allocate()
    {
        long totalLength = 0;
        for (int nn = 0; nn < mRunlines.size(); nn++)
            totalLength += mRunlines[nn].length();
        
        mPhi = std::vector<Precision::Float>(totalLength, 0.0);
//        LOG << "allocated.\n";
    }
    
    
    unsigned long bytes() const
    {
        long totalBytes = 0;
        totalBytes += mRunlines.size()*sizeof(RunlineJMPML);
        totalBytes += sizeof(Precision::Float)*m_c_hj.size();
        totalBytes += sizeof(Precision::Float)*m_c_phij.size();
        totalBytes += sizeof(Precision::Float)*m_c_hphi.size();
        totalBytes += sizeof(Precision::Float)*m_c_jphi.size();
        totalBytes += sizeof(Precision::Float)*mPhi.size();
        
        return totalBytes;
    }
    
    void apply(long timestep, Precision::Float dt)
    {
        for (int rr = 0; rr < mRunlines.size(); rr++)
        {
            Precision::Float* jm = mHeadCurrent + mRunlines[rr].mCurrent;
            Precision::Float* phi = &(mPhi.at(mRunlines[rr].mAux));
            const Precision::Float* ehLow = mHeadNeighborField + mRunlines[rr].mFieldLow;
            const Precision::Float* ehHigh = mHeadNeighborField + mRunlines[rr].mFieldHigh;
            const Precision::Float* c_hj = &(m_c_hj.at(mRunlines[rr].mConsts));
            const Precision::Float* c_phij = &(m_c_phij.at(mRunlines[rr].mConsts));
            const Precision::Float* c_hphi = &(m_c_hphi.at(mRunlines[rr].mConsts));
            const Precision::Float* c_jphi = &(m_c_jphi.at(mRunlines[rr].mConsts));
            unsigned long stride = mRunlines[rr].mStride;
            
            const Precision::Float sign = mSign;
            for (int ll = 0; ll < mRunlines[rr].length(); ll++)
            {
                Precision::Float diffEH = *ehHigh - *ehLow;
                Precision::Float curr = *c_hj * diffEH + 
                    *c_phij * (*phi);
//                if (diffEH != 0)
//                {
//                std::cerr << "chj "  << *c_hj << " diffEH " << diffEH
//                    << " c_phij " << *c_phij << " phi " << *phi << " curr "
//                    << sign*curr << "\n";
//                std::cerr << " adding " << sign*curr << "\n";
//                }
                *jm += sign*curr;
                *phi += *c_hphi * diffEH - *c_jphi * curr;
                
                jm++;
                phi++;
                ehLow++;
                ehHigh++;
                c_hj += stride;
                c_phij += stride;
                c_hphi += stride;
                c_jphi += stride;
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
    
private:
    Field mFieldCurrent, mFieldNeighbor;
    
    std::vector<RunlineJMPML> mRunlines;
    std::vector<Precision::Float> m_c_hj;
    std::vector<Precision::Float> m_c_phij;
    std::vector<Precision::Float> m_c_hphi;
    std::vector<Precision::Float> m_c_jphi;
    
    std::vector<Precision::Float> mPhi;
    Precision::Float mSign;
    
    Precision::Float* mHeadCurrent;
    Precision::Float* mHeadNeighborField;
};






#endif
