/*
 *  ElectromagneticFieldIndices.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 5/15/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _ELECTROMAGNETICFIELDINDICES_
#define _ELECTROMAGNETICFIELDINDICES_

#include "rle/SupportRegion3.h"
#include "rle/IndexArray3.h"
#include "VectorMatrix2.h"
#include "FieldEnum.h"
#include <map>

class NodeExtents;
class ConstitutiveTensorField;
class PMLField;

class IndexTensorField
{
public:
    IndexTensorField() {}
    
    RLE::IndexArray3 & operator()(int i, int j)
    {
        assert(i >= 0 && j >= 0 && i < 3 && j < 3);
        return mIndices[i][j];
    }
    
    const RLE::IndexArray3 & operator()(int i, int j) const
    {
        assert(i >= 0 && j >= 0 && i < 3 && j < 3);
        return mIndices[i][j];
    }
    
    long bytes() const
    {
        long total = 0;
        for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
            total += mIndices[ii][jj].bytes();
        return total;
    }
    
private:
    RLE::IndexArray3 mIndices[3][3];
};

class ElectromagneticFieldIndices
{
public:
    ElectromagneticFieldIndices() {}
    
    void index(const NodeExtents & extents,
        const ConstitutiveTensorField & inversePermittivity,
        const ConstitutiveTensorField & inversePermeability);
    
    void ee(int i, int j, const RLE::SupportRegion3 & support);
    void hh(int i, int j, const RLE::SupportRegion3 & support);
    void e(int xyz, const RLE::SupportRegion3 & support);
    void h(int xyz, const RLE::SupportRegion3 & support);
    void d(int xyz, const RLE::SupportRegion3 & support);
    void b(int xyz, const RLE::SupportRegion3 & support);
    void j(int xyz, const RLE::SupportRegion3 & support);
    void m(int xyz, const RLE::SupportRegion3 & support);
    void je(int xyz, const RLE::SupportRegion3 & support);
    void mh(int xyz, const RLE::SupportRegion3 & support);
    
    const RLE::IndexArray3 & ee(int i, int j) const { return mIndicesEE(i,j); }
    const RLE::IndexArray3 & hh(int i, int j) const { return mIndicesHH(i,j); }
    const RLE::IndexArray3 & e(int xyz) const { return mIndicesE[xyz]; }
    const RLE::IndexArray3 & h(int xyz) const { return mIndicesH[xyz]; }
    const RLE::IndexArray3 & d(int xyz) const { return mIndicesD[xyz]; }
    const RLE::IndexArray3 & b(int xyz) const { return mIndicesB[xyz]; }
    const RLE::IndexArray3 & j(int xyz) const { return mIndicesSourceJ[xyz]; }
    const RLE::IndexArray3 & m(int xyz) const { return mIndicesSourceM[xyz]; }
    const RLE::IndexArray3 & je(int xyz) const { return mIndicesSourceJ_E[xyz];}
    const RLE::IndexArray3 & mh(int xyz) const { return mIndicesSourceM_H[xyz];}
    
    const RLE::IndexArray3 & field(FieldName fName, int xyz) const;
    const RLE::IndexArray3 & field(FieldName fName, int ii, int jj) const;
    
    long bytes() const;
private:
    IndexTensorField mIndicesEE;
    IndexTensorField mIndicesHH;
    RLE::IndexArray3 mIndicesE[3];
    RLE::IndexArray3 mIndicesH[3];
    RLE::IndexArray3 mIndicesD[3];
    RLE::IndexArray3 mIndicesB[3];
    RLE::IndexArray3 mIndicesSourceJ[3];
    RLE::IndexArray3 mIndicesSourceM[3];
    RLE::IndexArray3 mIndicesSourceJ_E[3]; // Adjoint source for E update
    RLE::IndexArray3 mIndicesSourceM_H[3]; // Adjoint source for H update
};

class MaterialIndices
{
public:
    MaterialIndices() {}
    MaterialIndices(const ConstitutiveTensorField & inversePermittivity,
        const ConstitutiveTensorField & inversePermeability);
    
    void index(const ConstitutiveTensorField & inversePermittivity,
        const ConstitutiveTensorField & inversePermeability);
    
    void permittivity(int i, int j, int numerOrder, int denomOrder,
        const RLE::SupportRegion3 & support);
    void permeability(int i, int j, int numerOrder, int denomOrder,
        const RLE::SupportRegion3 & support);
    
    // Warning: returned by value (thanks a lot, std::map)
    RLE::IndexArray3 permittivity(int i, int j, int numerOrder,
        int denomOrder) const;
    RLE::IndexArray3 permeability(int i, int j, int numerOrder,
        int denomOrder) const;
    
    long bytes() const;
private:
    std::map<Vector2i, RLE::IndexArray3> mPermittivityIndices[3][3];
    std::map<Vector2i, RLE::IndexArray3> mPermeabilityIndices[3][3];
};

class PMLIndices
{
public:
    PMLIndices() {}
    
    void index(const NodeExtents & extents,
        const ConstitutiveTensorField & inversePermittivity,
        const ConstitutiveTensorField & inversePermeability,
        const PMLField & pmlE,
        const PMLField & pmlH);
    
    void e(int fieldXYZ, int absorptionXYZ,
        const RLE::SupportRegion3 & support);
    void h(int fieldXYZ, int absorptionXYZ,
        const RLE::SupportRegion3 & support);
    
    const RLE::IndexArray3 e(int fieldXYZ, int absorptionXYZ) const;
    const RLE::IndexArray3 h(int fieldXYZ, int absorptionXYZ) const;
    
    long bytes() const;
private:
    RLE::IndexArray3 mPMLE[3][3]; // [field][absorption direction]
    RLE::IndexArray3 mPMLH[3][3]; // [field][absorption direction]
};

class GhostIndices
{
public:
    GhostIndices() {}
    
private:
};



#endif

