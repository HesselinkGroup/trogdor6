/*
 *  HuygensSurface.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 11/21/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "HuygensSurface.h"

#include <sstream>

using namespace std;

HuygensSurface::
HuygensSurface(GridDescPtr sourceGrid, double factor,
    Field writeField, Field readField,
    const vector<RunlineHuygens> & runlines) :
    mSourceGrid(sourceGrid),
    mWriteField(writeField),
    mReadField(readField),
    mRunlines(runlines),
    mFactor(factor)
{
    std::ostringstream str;
    
    if (writeField.isElectric())
    {
        str << "Huygens E" << char('x' + writeField.xyz());
    }
    else
    {
        str << "Huygens H" << char('x' + writeField.xyz());
    }
    name(str.str());
}

void HuygensSurface::
setPointers(GridFields & currentGridFields,
    std::map<int, Pointer<GridFields> > & allGridFields)
{
    GridFields & sourceGridFields = *(allGridFields.find(sourceGrid()->id())->second);
    mHeadJM = currentGridFields.field(mWriteField);
    mHeadHE = sourceGridFields.field(mReadField);
}

//void HuygensSurface::
//setPointers(Precision::Float* headJM, Precision::Float* headHE)
//{
//    mHeadJM = headJM;
//    mHeadHE = headHE;
//}

unsigned long HuygensSurface::
bytes() const
{
    long totalBytes = 0;
    totalBytes += mRunlines.size()*sizeof(RunlineHuygens);
    
    return totalBytes;
}

void HuygensSurface::
apply(long timestep, Precision::Float dt)
{
    for (int rr = 0; rr < mRunlines.size(); rr++)
    {
//        LOG << "Runline " << rr << "\n";
        Precision::Float* jm = mHeadJM + mRunlines[rr].jm;
        Precision::Float* he = mHeadHE + mRunlines[rr].he;
        
        const int srcStride = mRunlines[rr].srcStride;
        const int destStride = mRunlines[rr].destStride;
        
        const Precision::Float factor = mFactor;
        
        for (int ll = 0; ll < mRunlines[rr].length(); ll++)
        {
            *jm += *he * factor;
            he += srcStride;
            jm += destStride;
        }
    }
}

void HuygensSurface::
printRunlines(std::ostream & str) const
{
    for (int nn = 0; nn < mRunlines.size(); nn++)
        str << mRunlines.at(nn) << "\n";
}
