/*
 *  GridFields.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/18/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "GridFields.h"
#include "ElectromagneticFieldIndices.h"
#include "Log.h"

using namespace std;

GridFields::
GridFields(const ElectromagneticFieldIndices & fieldIndices)
{
    // Allocate the E and H tensor fields
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
        mEE[i][j] = vector<Precision::Float>(fieldIndices.ee(i,j).length(), 0.0);
        mHH[i][j] = vector<Precision::Float>(fieldIndices.hh(i,j).length(), 0.0);
    }
    
    // Allocate the xyz fields
    for (int xyz = 0; xyz < 3; xyz++)
    {
        mE[xyz] = vector<Precision::Float>(fieldIndices.e(xyz).length(), 0.0);
        mH[xyz] = vector<Precision::Float>(fieldIndices.h(xyz).length(), 0.0);
        mD[xyz] = vector<Precision::Float>(fieldIndices.d(xyz).length(), 0.0);
        mB[xyz] = vector<Precision::Float>(fieldIndices.b(xyz).length(), 0.0);
        mJ[xyz] = vector<Precision::Float>(fieldIndices.j(xyz).length(), 0.0);
        mM[xyz] = vector<Precision::Float>(fieldIndices.m(xyz).length(), 0.0);
        mJ_E[xyz] = vector<Precision::Float>(fieldIndices.je(xyz).length(), 0.0);
        mM_H[xyz] = vector<Precision::Float>(fieldIndices.mh(xyz).length(), 0.0);
    }
}

Precision::Float* GridFields::
field(FieldName fName, int xyz)
{
    switch (fName)
    {
    case kE: return e(xyz);
    case kD: return d(xyz);
    case kJ: return j(xyz);
    case kH: return h(xyz);
    case kB: return b(xyz);
    case kM: return m(xyz);
    default: throw(std::logic_error("Bad field"));
    };
}

Precision::Float* GridFields::
field(Field f)
{
    switch (f.name())
    {
    case kE: return e(f.xyz());
    case kEE: return ee(f.i(), f.j());
    case kD: return d(f.xyz());
    case kJ: return j(f.xyz());
    case kJE: return je(f.xyz());
    case kH: return h(f.xyz());
    case kHH: return hh(f.i(), f.j());
    case kB: return b(f.xyz());
    case kM: return m(f.xyz());
    case kMH: return mh(f.xyz());
    default: throw(std::logic_error("Bad field"));
    };
}

Precision::Float* GridFields::
field(FieldName fName, int i, int j)
{
    switch(fName)
    {
    case kEE: return ee(i,j);
    case kHH: return hh(i,j);
    default: throw(std::logic_error("Bad field"));
    };
}

long GridFields::
bytes() const
{
    long eh = 0,
        db = 0,
        ee = 0,
        hh = 0,
        current = 0;
    long b = sizeof(Precision::Float);
    for (int xyz = 0; xyz < 3; xyz++)
    {
        eh += b*mE[xyz].size();
        eh += b*mH[xyz].size();
        db += b*mD[xyz].size();
        db += b*mB[xyz].size();
        current += b*mJ[xyz].size();
        current += b*mM[xyz].size();
        current += b*mJ_E[xyz].size();
        current += b*mM_H[xyz].size();
        
        for (int ijk = 0; ijk < 3; ijk++)
        {
            ee += b*mEE[xyz][ijk].size();
            hh += b*mHH[xyz][ijk].size();
        }
    }
    
    long totalBytes = eh + db + ee + hh + current;
    
    LOG << "\teh " << eh/1024.0 << " kB\n"
        << "\tdb " << db/1024.0 << " kB\n"
        << "\tee " << ee/1024.0 << " kB\n"
        << "\thh " << hh/1024.0 << " kB\n"
        << "\tjm " << current/1024.0 << " kB\n";
    
    return totalBytes;
}





