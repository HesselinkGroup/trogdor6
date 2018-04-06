/*
 *  ElectromagneticFieldIndices.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 5/15/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "ElectromagneticFieldIndices.h"

#include "UserPreferences.h"
#include "ConstitutiveTensorField.h"
#include "PMLField.h"
#include "Log.h"
#include "Extents.h"
#include "YeeUtilities.h"
#include "TimeWrapper.h"

using namespace RLE;
using namespace std;
using namespace YeeUtilities;

using TimeWrapper::now;
using TimeWrapper::TimePoint;
using TimeWrapper::elapsedMicroseconds;
using TimeWrapper::elapsedSeconds;

static SupportRegion3 sSupport(Rect3i r, Vector3b dims)
{
    SupportRegion3 supp(dims[0], dims[1], dims[2]);
    supp.mark(r.p1[0], r.p1[1], r.p1[2], r.p2[0], r.p2[1], r.p2[2]);
    return supp;
}

void ElectromagneticFieldIndices::
index(const NodeExtents & extents,
    const ConstitutiveTensorField & inversePermittivity,
    const ConstitutiveTensorField & inversePermeability)
{
//    Vector3b hasDim(extents.allocatedYeeCells().size());
//    Vector3b hasDim(extents.allocatedYeeCells().num(0) != 1,
//        extents.allocatedYeeCells().num(1) != 1, 
//        extents.allocatedYeeCells().num(2) != 1);
//    Vector3b hasDim(1,1,1); // because nnT and F begin as 3D, everything's 3D.
    Vector3b hasDim(extents.dimensions());
//    std::cout << "Dimensions: " << hasDim << "\n";
    
    // Forward case:    need E and H in all allocated cells
    //                  need D and B in all calculated cells (fewer)
    // Adjoint case:    need E and H in all calculated cells
    //                  need D and B in all allocated cells
    SupportRegion3 allYeeSupport = sSupport(extents.allocatedYeeCells(), hasDim);
    SupportRegion3 calcSupport;
    
    // D, B
    // The only cells that don't store D are the pure O(0,0) dielectrics in
    // regions that don't need D to be measured.
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (UserPreferences::defines("adjoint"))
        {
            d(xyz, allYeeSupport);
        }
        else
        {
            //calcSupport = sSupport(halfToYee(extents.calcHalfCells(), octantE(xyz)), hasDim);
            calcSupport = sSupport(halfToYee(extents.calcHalfCells(), octantE(xyz)), Vector3b(1,1,1));
            
            // Obtain the support of D that's necessary due to anisotropic epsilon.
            int jj = (xyz+1)%3, kk = (xyz+2)%3;
            SupportRegion3 offDiag = SupportRegion3(inversePermittivity(xyz,jj)) +
                SupportRegion3(inversePermittivity(xyz,kk));
            SupportRegion3 anisotropicSupport(offDiag), tx;
            translate(offDiag, tx, (-Vector3i::unit(xyz)).asArray());
            anisotropicSupport = anisotropicSupport + tx;
            
            SupportRegion3 support;
            if (anisotropicSupport.numRuns() > 0)
            {
                support = calcSupport - inversePermittivity.support(xyz, xyz, 0, 0) + anisotropicSupport;
            }
            else
            {
                support = calcSupport - inversePermittivity.support(xyz, xyz, 0, 0);
            }
            
            // AT THIS POINT support MUST HAVE THE RIGHT SYMMETRIC DIMENSIONS INDICATED
//            if (support.numRuns() > 0)
//            {
//                std::cout << "I HAVE DATA and my dimensions are " << support.numDimensions() << "\n";
//            }
            support.changeDimensions(hasDim.asArray());
            d(xyz, support);
        }
    }
    
    for (int xyz = 0; xyz < 3; xyz++)
    {
        calcSupport = sSupport(halfToYee(extents.calcHalfCells(), octantH(xyz)), Vector3b(1,1,1));
        if (UserPreferences::defines("adjoint"))
            b(xyz, allYeeSupport);
        else
        {
            calcSupport = sSupport(halfToYee(extents.calcHalfCells(), octantH(xyz)), Vector3b(1,1,1));
            SupportRegion3 support = calcSupport -
                inversePermeability.support(xyz, xyz, 0, 0);
            
            // AT THIS POINT support MUST HAVE THE RIGHT SYMMETRIC DIMENSIONS INDICATED
//            if (support.numRuns() > 0)
//            {
//                std::cout << "I HAVE DATA and my dimensions are " << support.numDimensions() << "\n";
//            }
            support.changeDimensions(hasDim.asArray());
            b(xyz, support);
        }
    }
    
    // E, H
    for (int xyz = 0; xyz < 3; xyz++)
    {
//        if (UserPreferences::defines("adjoint"))
//        {
//            calcSupport = sSupport(halfToYee(extents.calcHalfCells(),
//                octantE(xyz)), hasDim);
//            
//            SupportRegion3 support = calcSupport -
//                inversePermittivity.support(xyz, xyz, 0, 0);
//            e(xyz, support);
//        }
//        else
            e(xyz, allYeeSupport);
    }
    for (int xyz = 0; xyz < 3; xyz++)
    {
//        if (UserPreferences::defines("adjoint"))
//        {
//            calcSupport = sSupport(halfToYee(extents.calcHalfCells(),
//                octantH(xyz)), hasDim);
//            
//            SupportRegion3 support = calcSupport -
//                inversePermeability.support(xyz, xyz, 0, 0);
//            h(xyz, support);
//        }
//        else
            h(xyz, allYeeSupport);
    }
    
    // We need the Exx, Eyy etc. components whenever adjacent cells at offset
    // (0,0,0) have off-diagonal terms.  Same with H (also using (0,0,0)).
    
    // Obtain Exx, Eyy and Ezz.
    // Secondarily obtain support for Exy, Exz, etc.
    for (int xyz = 0; xyz < 3; xyz++)
    {
        int jj = (xyz+1)%3, kk = (xyz+2)%3;
        SupportRegion3 offDiag = SupportRegion3(inversePermittivity(xyz,jj)) +
            SupportRegion3(inversePermittivity(xyz,kk));
        SupportRegion3 wholeSupport(offDiag), tx;
        translate(offDiag, tx, (-Vector3i::unit(xyz)).asArray());
        wholeSupport = wholeSupport + tx;
        
        ee(xyz, xyz, wholeSupport);
        
        translate(offDiag, tx, Vector3i::unit(xyz).asArray());
        wholeSupport = wholeSupport + tx;
        
        // AT THIS POINT wholeSupport MUST HAVE THE RIGHT SYMMETRIC DIMENSIONS INDICATED
        
//        if (wholeSupport.numRuns() > 0)
//        {
//            std::cout << "I HAVE DATA and my dimensions are " << wholeSupport.numDimensions() << "\n";
//        }
        wholeSupport.changeDimensions(hasDim.asArray()); // such a stupid hack
        
        ee(xyz, jj, wholeSupport);
        ee(xyz, kk, wholeSupport);
    }
    
    // Obtain Hxx, Hyy and Hzz.
    // Secondarily obtain support for Hxy, Hxz, etc.
    for (int xyz = 0; xyz < 3; xyz++)
    {
        int jj = (xyz+1)%3, kk = (xyz+2)%3;
        SupportRegion3 offDiag = SupportRegion3(inversePermeability(xyz,jj)) +
            SupportRegion3(inversePermeability(xyz,kk));
        SupportRegion3 wholeSupport(offDiag), tx;
        translate(offDiag, tx, Vector3i::unit(xyz).asArray());
        wholeSupport = wholeSupport + tx;
        
        hh(xyz, xyz, wholeSupport);
        
        translate(offDiag, tx, (-Vector3i::unit(xyz)).asArray());
        wholeSupport = wholeSupport + tx;
        
        // AT THIS POINT wholeSupport MUST HAVE THE RIGHT SYMMETRIC DIMENSIONS INDICATED
        wholeSupport.changeDimensions(hasDim.asArray()); // STUPID HATEFUL HACK
        
        hh(xyz, jj, wholeSupport);
        hh(xyz, kk, wholeSupport);
    }
    
//    std::cout << "\n==== Indexing a grid:\n";
//    std::cout << extents.allocatedYeeCells() << "\n";
//    //std::cout << allYeeSupport << "\n";
//    std::cout << e(0) << "\n";
//    std::cout << e(0).numDimensions() << "dimensions \n";
    
}

void ElectromagneticFieldIndices::
ee(int i, int j, const RLE::SupportRegion3 & support)
{
    mIndicesEE(i,j) = IndexArray3(support);
}

void ElectromagneticFieldIndices::
hh(int i, int j, const RLE::SupportRegion3 & support)
{
    mIndicesHH(i,j) = IndexArray3(support);
}

void ElectromagneticFieldIndices::
e(int xyz, const RLE::SupportRegion3 & support)
{
    mIndicesE[xyz] = IndexArray3(support);
}

void ElectromagneticFieldIndices::
h(int xyz, const RLE::SupportRegion3 & support)
{
    mIndicesH[xyz] = IndexArray3(support);
}

void ElectromagneticFieldIndices::
d(int xyz, const RLE::SupportRegion3 & support)
{
    mIndicesD[xyz] = IndexArray3(support);
}

void ElectromagneticFieldIndices::
b(int xyz, const RLE::SupportRegion3 & support)
{
    mIndicesB[xyz] = IndexArray3(support);
}

void ElectromagneticFieldIndices::
j(int xyz, const RLE::SupportRegion3 & support)
{
    mIndicesSourceJ[xyz] = IndexArray3(support);
}

void ElectromagneticFieldIndices::
je(int xyz, const RLE::SupportRegion3 & support)
{
    mIndicesSourceJ_E[xyz] = IndexArray3(support);
}

void ElectromagneticFieldIndices::
m(int xyz, const RLE::SupportRegion3 & support)
{
    mIndicesSourceM[xyz] = IndexArray3(support);
}

void ElectromagneticFieldIndices::
mh(int xyz, const RLE::SupportRegion3 & support)
{
    mIndicesSourceM_H[xyz] = IndexArray3(support);
}

const IndexArray3 & ElectromagneticFieldIndices::
field(FieldName fName, int xyz) const
{
    switch (fName)
    {
        case kD: return d(xyz);
        case kE: return e(xyz);
        case kJ: return j(xyz);
        case kB: return b(xyz);
        case kH: return h(xyz);
        case kM: return m(xyz);
        default:
            throw(std::logic_error("Bad field"));
    };
}

const IndexArray3 & ElectromagneticFieldIndices::
field(FieldName fName, int ii, int jj) const
{
    switch (fName)
    {
        case kEE: return ee(ii,jj);
        case kHH: return hh(ii,jj);
        default:
            throw(std::logic_error("Bad field"));
    };
}

long ElectromagneticFieldIndices::
bytes() const
{
    long total = 0;
    total += mIndicesEE.bytes();
    total += mIndicesHH.bytes();
    
    for (int xyz = 0; xyz < 3; xyz++)
    {
        total += mIndicesD[xyz].bytes();
        total += mIndicesB[xyz].bytes();
        total += mIndicesE[xyz].bytes();
        total += mIndicesH[xyz].bytes();
        total += mIndicesSourceJ[xyz].bytes();
        total += mIndicesSourceM[xyz].bytes();
        total += mIndicesSourceJ_E[xyz].bytes();
        total += mIndicesSourceM_H[xyz].bytes();
    }
    return total;
}

MaterialIndices::
MaterialIndices(const ConstitutiveTensorField & inversePermittivity,
    const ConstitutiveTensorField & inversePermeability)
{
    index(inversePermittivity, inversePermeability);
}

void MaterialIndices::
index(const ConstitutiveTensorField & inversePermittivity,
    const ConstitutiveTensorField & inversePermeability)
{
    
    TimePoint t0 = now();
    LOG << "Sorting permittivities...\n";
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    {
        for (int numer = 0; numer < 20; numer++)
        for (int denom = 0; denom < 20; denom++)
        {
            SupportRegion3 supp = inversePermittivity.support(
                ii, jj, numer, denom);
            
            permittivity(ii, jj, numer, denom, supp);
        }
    }
    LOG << "Elapsed time " << elapsedSeconds(t0, now()) << " seconds.\n";
    
    t0 = now();
    LOG << "Sorting permeabilities...\n";
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    {
        for (int numer = 0; numer < 20; numer++)
        for (int denom = 0; denom < 20; denom++)
        {
            SupportRegion3 supp = inversePermeability.support(
                ii, jj, numer, denom);
            
            permeability(ii, jj, numer, denom, supp);
        }
    }
    LOG << "Elapsed time " << elapsedSeconds(t0, now()) << " seconds.\n";
}

void MaterialIndices::
permittivity(int i, int j, int numerOrder, int denomOrder,
    const RLE::SupportRegion3 & support)
{
    mPermittivityIndices[i][j][Vector2i(numerOrder, denomOrder)] =
        IndexArray3(support);
}

void MaterialIndices::
permeability(int i, int j, int numerOrder, int denomOrder,
    const RLE::SupportRegion3 & support)
{
    mPermeabilityIndices[i][j][Vector2i(numerOrder, denomOrder)] =
        IndexArray3(support);
}

RLE::IndexArray3 MaterialIndices::
permittivity(int i, int j, int numerOrder, int denomOrder) const
{
    std::map<Vector2i, IndexArray3>::const_iterator itr =
        mPermittivityIndices[i][j].find(Vector2i(numerOrder, denomOrder));
    if (itr != mPermittivityIndices[i][j].end())
        return itr->second;
    else
        return IndexArray3();
}

RLE::IndexArray3 MaterialIndices::
permeability(int i, int j, int numerOrder, int denomOrder) const
{
    std::map<Vector2i, IndexArray3>::const_iterator itr =
        mPermeabilityIndices[i][j].find(Vector2i(numerOrder, denomOrder));
    if (itr != mPermeabilityIndices[i][j].end())
        return itr->second;
    else
        return IndexArray3();
}

long MaterialIndices::
bytes() const
{
    std::map<Vector2i, IndexArray3>::const_iterator itr;
    
    long total = 0;
    
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    for (itr = mPermittivityIndices[ii][jj].begin();
        itr != mPermittivityIndices[ii][jj].end(); itr++)
    {
        total += itr->second.bytes();
//        if (itr->second.numRuns() > 0)
//            cerr << "Material " << itr->first << ", "
//            << itr->second.numRuns() << " runs.\n";
    }
    
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    for (itr = mPermeabilityIndices[ii][jj].begin();
        itr != mPermeabilityIndices[ii][jj].end(); itr++)
    {
        total += itr->second.bytes();
//        if (itr->second.numRuns() > 0)
//            cerr << "Material " << itr->first << ", "
//            << itr->second.numRuns() << " runs.\n";
    }
    
    return total;
}

void PMLIndices::
index(const NodeExtents & extents, 
    const ConstitutiveTensorField & inversePermittivity,
    const ConstitutiveTensorField & inversePermeability,
    const PMLField & pmlE,
    const PMLField & pmlH)
{
    for (int fieldXYZ = 0; fieldXYZ < 3; fieldXYZ++)
    for (int absorbXYZ = 0; absorbXYZ < 3; absorbXYZ++)
    if (extents.dimensions()[absorbXYZ] != false) // exclude unused dimensions explicitly
    {
        Rect3i eYee = halfToYee(extents.calcHalfCells(), octantE(fieldXYZ));
        Rect3i hYee = halfToYee(extents.calcHalfCells(), octantH(fieldXYZ));
        
        SupportRegion3 supportCalcE(extents.physicalYeeCells().num().asArray());
        supportCalcE.mark(eYee.p1[0], eYee.p1[1], eYee.p1[2],
            eYee.p2[0], eYee.p2[1], eYee.p2[2]);
        SupportRegion3 supportCalcH(extents.physicalYeeCells().num().asArray());
        supportCalcH.mark(hYee.p1[0], hYee.p1[1], hYee.p1[2],
            hYee.p2[0], hYee.p2[1], hYee.p2[2]);
        
        SupportRegion3 supportE = supportCalcE *
            SupportRegion3(pmlE.along(absorbXYZ)[fieldXYZ]);
        e(fieldXYZ, absorbXYZ, supportE);
        
        SupportRegion3 supportH = supportCalcH *
            SupportRegion3(pmlH.along(absorbXYZ)[fieldXYZ]);
        h(fieldXYZ, absorbXYZ, supportH);
    }
}

void PMLIndices::
e(int fieldXYZ, int absorptionXYZ, const RLE::SupportRegion3 & support)
{
//    cerr << "Set " << fieldXYZ << " " << absorptionXYZ << "\n";
    mPMLE[fieldXYZ][absorptionXYZ] = IndexArray3(support);
//    cerr << "\t" << mPMLE[fieldXYZ][absorptionXYZ].bytes() << " bytes\n";
}

void PMLIndices::
h(int fieldXYZ, int absorptionXYZ, const RLE::SupportRegion3 & support)
{
//    cerr << "Set " << fieldXYZ << " " << absorptionXYZ << "\n";
    mPMLH[fieldXYZ][absorptionXYZ] = IndexArray3(support);
//    cerr << "\t" << mPMLH[fieldXYZ][absorptionXYZ].bytes() << " bytes\n";
}

const RLE::IndexArray3 PMLIndices::
e(int fieldXYZ, int absorptionXYZ) const
{
    return mPMLE[fieldXYZ][absorptionXYZ];
}

const RLE::IndexArray3 PMLIndices::
h(int fieldXYZ, int absorptionXYZ) const
{
    return mPMLH[fieldXYZ][absorptionXYZ];
}

long PMLIndices::
bytes() const
{
    long total = 0;
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    {
        total += mPMLE[ii][jj].bytes();
        total += mPMLH[ii][jj].bytes();
//        cerr << "PML runs " << mPMLE[ii][jj].numRuns() << " and "
//            << mPMLH[ii][jj].numRuns() << "\n";
    }
    return total;
}


