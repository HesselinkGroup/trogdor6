/*
 *  TensorGridConstants.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 4/14/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "TensorGridConstants.h"

#include "Log.h"
#include "TimeWrapper.h"
#include "YeeUtilities.h"
#include "RLEOperations.h"
#include <map>
#include <fstream>
#include <stdexcept>
#include <stdint.h>
#include "WriteSomeRLE.h"
#include "rle/IndexArray3.h"
#include "rle/MultiApply.h"
#include "UserPreferences.h"

using namespace std;
using namespace YeeUtilities;
using namespace RLE;

static inline Precision::RationalFunction sTimes(const Precision::RationalFunction & r,
    Precision::Float d)
{
    return d*r;
}

// Make the denominator polynomial lead with a unit constant term.
static inline Precision::RationalFunction
sStandardForm(const Precision::RationalFunction & r)
{
    assert(r.denominator()[0] != 0);
    assert(r.denominator().order() == r.numerator().order());
    Precision::Float d0_inv = 1.0/r.denominator()[0];
    
    
    return Precision::RationalFunction(d0_inv*r.numerator(),
        d0_inv*r.denominator());
}

#pragma mark *** New Stuff ***


TensorGridConstants::
TensorGridConstants(int i, int j,
    DynamicRLE3<Precision::RationalFunction> & invConstits,
    Rect3i voxelBounds,
    Precision::RationalFunction backgroundMaterial) :
    mInverseConstitutives(invConstits),
    mBackground(backgroundMaterial),
    mVoxels(voxelBounds),
    m_i(i),
    m_j(j)
{
}

void TensorGridConstants::
calculateConstitutives(
    const std::vector<Precision::RationalFunction> & materials,
    const std::vector<RLE::DynamicRLE3<double> > & fillFactors,
    const RLE::DynamicRLE3<double> & orientations)
{
    if (UserPreferences::defines("noAveraging"))
        mInverseConstitutives = winnerTakesAll(materials, fillFactors);
    else if (m_i == m_j)
    {
        mInverseConstitutives = diagonalElement(materials, fillFactors,
            orientations);
    }
    else
    {
        mInverseConstitutives = offDiagonalElement(materials,
            fillFactors, orientations);
    }
}

void TensorGridConstants::
writeSensitivities(
    const std::vector<Precision::RationalFunction> & materials,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices,
    const std::vector<RLE::DynamicRLE3<double> > & fillFactors,
    const RLE::DynamicRLE3<double> & orientations,
    const std::vector<std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > > & fillFactorJacobians,
    const std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > & orientationJacobians,
    const SupportRegion3 & savedFields,
    std::ostream & file)
{
    uint32_t ii(m_i), jj(m_j);
    file.write("TENS", 4);
    file.write((char*)&ii, sizeof(uint32_t));
    file.write((char*)&jj, sizeof(uint32_t));
    
    uint32_t numV = controlVertices.size();
    file.write((char*)&numV, sizeof(uint32_t));
    
    for (uint32_t vv = 0; vv < controlVertices.size(); vv++)
    if (controlVertices.at(vv).freeDirections() != Vector3b(0,0,0))
    {
        file.write("VERT", 4);
        file.write((char*)&vv, sizeof(uint32_t));
        
        for (uint32_t freeXYZ = 0; freeXYZ < 3; freeXYZ++)
        {
            file.write((char*)"DIRE", 4);
            file.write((char*)&freeXYZ, sizeof(uint32_t));
            
            if (controlVertices.at(vv).freeDirections()[freeXYZ])
            {
                // extract the right bits of the sensitivities.
                
                vector<DynamicRLE3<double> > dFillFactors;
                DynamicRLE3<double>
                    dOrientation_unmasked(orientations.capacity()),
                    dOrientation(orientations.capacity());
                
                selectFillFactorSensitivity(fillFactorJacobians, vv,
                    freeXYZ, dFillFactors);
                selectOrientationSensitivity(orientationJacobians, vv,
                    freeXYZ, dOrientation_unmasked);
                
                // IMPORTANT: mask dOrientation!
                restriction(dOrientation_unmasked, orientations, dOrientation);
                
                Pointer<DynamicRLE3<Precision::RationalFunction> > sensitivity;
                
                if (m_i == m_j)
                {
                    sensitivity = diagonalSensitivity(materials, fillFactors,
                        orientations, dFillFactors, dOrientation);
                }
                else
                {
                    sensitivity = offDiagonalSensitivity(materials, fillFactors,
                        orientations, dFillFactors, dOrientation);
                }
                
                writeSensitivity(sensitivity, savedFields, file);
            }
            else
            {
                file.write("NENDDEND", 8);
            }
        }
    }
    file.write("TEND", 4);
}

static void sCheckNegative(const DynamicRLE3<Precision::RationalFunction> & rle)
{
    DynamicRLE3<Precision::RationalFunction>::ConstIterator itr;
    
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        const Precision::RationalFunction & rf = itr.mark();
        
        if (rf.numerator().order() == 0 && rf.denominator().order() == 0)
        {
            if (rf.numerator()[0] * rf.denominator()[0] < 0)
            {
                long x, y, z;
                cerr << "Negative value!  Where? ";
                rle.cartesianCoordinates(itr.position(), x, y, z);
                cerr << x << " , " << y << ", " << z << "\n";
                throw(std::logic_error("Found negative number."));
            }
        }
    }
}

static void sCheckOrder(const DynamicRLE3<Precision::RationalFunction> & rle)
{
    DynamicRLE3<Precision::RationalFunction>::ConstIterator itr;

    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        const Precision::RationalFunction & rf = itr.mark();
        
        if (rf.numerator().order() != rf.denominator().order())
        {
            std::cout << "Found material " << rf << ".\n";
        }
    }
}


struct WinnerTakesAll
{
    WinnerTakesAll(
        const vector<Precision::RationalFunction> & constitutives,
        DynamicRLE3<Precision::RationalFunction> & outArithmeticMean) :
        mConstitutives(constitutives),
        mOut(outArithmeticMean)
    {
    }
    
    void operator()(long first, long last, double const* fillFactor,
        int const* markedArrays, int numStreams)
    {
        
        Precision::Float largestFillFactor = 0.0;
        int indexLargest = -1;
        
        for (int nn = 0; nn < numStreams-1; nn++)
        if (markedArrays[nn] != false)
        {
            if (fillFactor[nn] > largestFillFactor)
            {
                largestFillFactor = fillFactor[nn];
                indexLargest = nn;
            }
        }
        
        if (indexLargest == -1)
            throw(std::logic_error("Looks like there's no material here"));
        
        mOut.RLEBase<DRLESegment<Precision::RationalFunction> >::mark(
            first, last, mConstitutives[indexLargest]);
    }
    
    const vector<Precision::RationalFunction> & mConstitutives;
    DynamicRLE3<Precision::RationalFunction> & mOut;
};

struct WinnerTakesAllWithBackground
{
    WinnerTakesAllWithBackground(
        const vector<Precision::RationalFunction> & constitutives,
        Precision::RationalFunction background,
        DynamicRLE3<Precision::RationalFunction> & outArithmeticMean) :
        mConstitutives(constitutives),
        mBackground(background),
        mOut(outArithmeticMean)
    {
    }
    
    void operator()(long first, long last, double const* fillFactor,
        int const* markedArrays, int numStreams)
    {
        Precision::RationalFunction sum;
        Precision::Float backgroundFill = 1.0;
        
        Precision::Float largestFillFactor = 0.0;
        int indexLargest = -1;
        
        for (int nn = 0; nn < numStreams; nn++)
        if (markedArrays[nn] != false)
        {
//            if (fillFactor[1] > 0.0)
//                int farb = 'a';
            if (fillFactor[nn] > largestFillFactor)
            {
                largestFillFactor = fillFactor[nn];
                indexLargest = nn;
            }
            backgroundFill -= fillFactor[nn];
        }
        
//        cerr << "from " << first << " to " << last << " fills: ";
//        for (int mm = 0; mm < numStreams; mm++)
//        {
//            cerr << "[";
//            if (markedArrays[mm])
//                cerr << fillFactor[mm] << ": marked] ";
//            else
//                cerr << "unmarked] ";
//        }
//        cerr << "\n";
//        
//        if (numStreams == 2 && markedArrays[0] && markedArrays[1])
//        {
//            if (fillFactor[1] > fillFactor[0])
//            {
//                cerr << "got one!\n";
//            }
//        }
        
        if (backgroundFill > largestFillFactor)
            mOut.RLEBase<DRLESegment<Precision::RationalFunction> >::mark(
                first, last, mBackground);
        else
        {
            if (indexLargest == -1)
                throw(std::logic_error("Looks like there's no material here"));
//            if (indexLargest != 0)
//                int flarg = 'b';
            mOut.RLEBase<DRLESegment<Precision::RationalFunction> >::mark(
                first, last, mConstitutives[indexLargest]);
        }
    }
    
    const vector<Precision::RationalFunction> & mConstitutives;
    Precision::RationalFunction mBackground;
    DynamicRLE3<Precision::RationalFunction> & mOut;
};


DynamicRLE3<Precision::RationalFunction> TensorGridConstants::
winnerTakesAll(
    const vector<Precision::RationalFunction> & materials,
    const vector<DynamicRLE3<double> > & fillFactors)
{
    vector<Precision::RationalFunction> inverseTensors;
    std::transform(materials.begin(), materials.end(),
        std::back_inserter(inverseTensors),
        reciprocal<Precision::Float>);
        
    vector<DynamicRLE3<double> const*> fillFactorPointers;
    for (int nn = 0; nn < fillFactors.size(); nn++)
        fillFactorPointers.push_back(&fillFactors[nn]);
    
    DynamicRLE3<Precision::RationalFunction> outInvPerm, monicDenom;
//    outInvPerm.dimensions(voxels().num().asArray());
    
    if (backgroundMaterial() != Precision::Float(0.0))
    {
        outInvPerm.mark(mVoxels.p1[0], mVoxels.p1[1], mVoxels.p1[2],
            mVoxels.p2[0], mVoxels.p2[1], mVoxels.p2[2],
            reciprocal(backgroundMaterial()));
        
        WinnerTakesAllWithBackground callback(inverseTensors, 
            reciprocal(backgroundMaterial()), outInvPerm);
        RLE::multiApply(&(fillFactorPointers[0]), fillFactors.size(), callback,
            RLE::ExcludeGaps);
    }
    else
    {
        WinnerTakesAll callback(inverseTensors, outInvPerm);
        RLE::multiApply(&(fillFactorPointers[0]), fillFactors.size(), callback,
            RLE::ExcludeGaps);
    }
    
    transform(outInvPerm, monicDenom, sStandardForm);
    
    return monicDenom;
}

DynamicRLE3<Precision::RationalFunction> TensorGridConstants::
diagonalElement(
    const vector<Precision::RationalFunction> & materials,
    const vector<DynamicRLE3<double> > & fillFactors,
    const DynamicRLE3<double> & orientation_ij)
{
    vector<Precision::RationalFunction> inverseTensors;
    std::transform(materials.begin(), materials.end(),
        std::back_inserter(inverseTensors),
        reciprocal<Precision::Float>);
    
//    for (int tt = 0; tt < inverseTensors.size(); tt++)
//    {
//        std::cout << "MATERIAL " << tt << ":" << inverseTensors[tt] << "\n";
//    }
    
    DynamicRLE3<Precision::RationalFunction> monicDenom(orientation_ij.nonSymmetricDimensions()),
        arithmeticAddend(orientation_ij.nonSymmetricDimensions()),
        harmonicAddend(orientation_ij.nonSymmetricDimensions());
    
    // Find the arithmetic average of inverse(constitutive tensor)
    // and the harmonic average of inverse(constitutive tensor).
    {
        DynamicRLE3<Precision::RationalFunction> arithmeticMean(orientation_ij.nonSymmetricDimensions());
        calculateArithmeticMean(fillFactors, inverseTensors, arithmeticMean);
        
        sCheckOrder(arithmeticMean);
        //cerr << "Arith " << arithmeticMean << "\n\n";
        transform(arithmeticMean, orientation_ij, arithmeticAddend, sTimes, Precision::RationalFunction(0), 0.0);
        
        sCheckOrder(arithmeticAddend);
        //checkNegative(arithmeticAddend);
    }
    
    {
        DynamicRLE3<Precision::RationalFunction> harmonicMean(orientation_ij.nonSymmetricDimensions());
        calculateHarmonicMean(fillFactors, inverseTensors, harmonicMean);
        sCheckOrder(harmonicMean);
//        cerr << "Harm " << harmonicMean << "\n\n";
        transform(harmonicMean, -orientation_ij, harmonicAddend, sTimes,
            Precision::RationalFunction(0), 0.0);
        harmonicAddend = harmonicMean + harmonicAddend;
        
        sCheckOrder(harmonicAddend);
        //checkNegative(harmonicMean);
        //checkNegative(harmonicAddend);
    }
    
    transform(harmonicAddend + arithmeticAddend, monicDenom, sStandardForm);
    
    sCheckOrder(monicDenom);
    
//    cerr << "Diagonal element: " << monicDenom << "\n";
//    for (int ff = 0; ff < fillFactors.size(); ff++)
//        cerr << "From fill " << ff << ": " << fillFactors[ff] << "\n";
//    cerr << "From o: " << orientation_ij << "\n";
    
//    std::cout << monicDenom << "\n";
    
    return monicDenom;
}


DynamicRLE3<Precision::RationalFunction> TensorGridConstants::
offDiagonalElement(
    const vector<Precision::RationalFunction> & materials,
    const vector<DynamicRLE3<double> > & fillFactors,
    const DynamicRLE3<double> & orientation_ij)
{
    vector<Precision::RationalFunction> inverseTensors;
    std::transform(materials.begin(), materials.end(),
        std::back_inserter(inverseTensors),
        reciprocal<Precision::Float>);
    
    DynamicRLE3<Precision::RationalFunction> monicDenom, arithmeticAddend,
        harmonicAddend, returnVal;
    
    // Find the arithmetic average of inverse(constitutive tensor)
    // and the harmonic average of inverse(constitutive tensor).
    {
        DynamicRLE3<Precision::RationalFunction> arithmeticMean;
        calculateArithmeticMean(fillFactors, inverseTensors, arithmeticMean);
        transform(arithmeticMean, orientation_ij, arithmeticAddend, sTimes,
            Precision::RationalFunction(0), 0.0);
    }
    
    {
        DynamicRLE3<Precision::RationalFunction> harmonicMean;
        calculateHarmonicMean(fillFactors, inverseTensors, harmonicMean);
        transform(harmonicMean, -orientation_ij, harmonicAddend, sTimes,
            Precision::RationalFunction(0), 0.0);
    }
    
    if (false == UserPreferences::defines("flipOffDiagonalSign"))
        transform(harmonicAddend + arithmeticAddend, monicDenom, sStandardForm);
    else
        transform(-harmonicAddend - arithmeticAddend, monicDenom, sStandardForm);
    
    filter(monicDenom, returnVal, IsNonzero());
    
    return returnVal;
}

static void checkNaN(const DynamicRLE3<double> & rle)
{
    DynamicRLE3<double>::ConstIterator itr;
    
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        if (isnan(itr.mark()))
        {
            long x,y,z;
            
            cerr << "Not a number!  Where? ";
            rle.cartesianCoordinates(itr.position(), x, y, z);
            cerr << x << ", " << y << ", " << z << "\n";
            throw(std::logic_error("Found NaN"));
        }
    }
}

Pointer<RLE::DynamicRLE3<Precision::RationalFunction> >
TensorGridConstants::
diagonalSensitivity(
    const vector<Precision::RationalFunction> & materials,
    const vector<DynamicRLE3<double> > & fillFactors,
    const DynamicRLE3<double> & orientations,
    const vector<DynamicRLE3<double> > & fillFactorSensitivities,
    const DynamicRLE3<double> & orientationSensitivities)
{
    assert(fillFactors.size() == materials.size());
    assert(fillFactorSensitivities.size() == materials.size());
    
    //const Precision::Float DELTA = sqrt(f::epsilon()); // for finite-differencing.
    //const Precision::Float DELTA = mDerivDelta; //1e-8;
    Precision::Float DELTA;
    if (UserPreferences::defines("derivDelta"))
        DELTA = UserPreferences::valueAs<double>("derivDelta");
    else if (sizeof(Precision::Float) == sizeof(double))
        DELTA = 1e-4;
    else
        DELTA = 1e-2;
    
    // Obtain the total region that can change (it's the total support of the
    // fill factors and orientation)
    SupportRegion3 vertexInfluence;
    if (fillFactorSensitivities.size() > 0)
        vertexInfluence = SupportRegion3(fillFactorSensitivities[0]);
    for (int nn = 1; nn < fillFactorSensitivities.size(); nn++)
        vertexInfluence += SupportRegion3(fillFactorSensitivities[nn]);
    if (vertexInfluence.numRuns() > 0)
        vertexInfluence += SupportRegion3(orientationSensitivities);
    else
        vertexInfluence = SupportRegion3(orientationSensitivities);
    
    // Declare temporary arrays
    DynamicRLE3<Precision::RationalFunction> constit1, constit2;
    vector<DynamicRLE3<double> > fillFactors0(fillFactors.size()),
        fillFactors1(fillFactors.size()),
        fillFactors2(fillFactors.size());
    DynamicRLE3<double> orientation0, orientation1, orientation2;
    
    // Obtain perturbed fill factors
    for (int ff = 0; ff < fillFactors.size(); ff++)
    {
        restriction(fillFactors[ff], vertexInfluence, fillFactors0[ff]);
        fillFactors1[ff] = fillFactors0[ff] - DELTA*fillFactorSensitivities[ff];
        fillFactors2[ff] = fillFactors0[ff] + DELTA*fillFactorSensitivities[ff];
//        cerr << "fill1[" << ff << "] = " << fillFactors1[ff] << "\n\n";
//        cerr << "fill2[" << ff << "] = " << fillFactors2[ff] << "\n\n";

        checkNaN(fillFactors[ff]);
        checkNaN(fillFactors1[ff]);
        checkNaN(fillFactors2[ff]);
    }
    
    // Obtain perturbed orientations
    restriction(orientations, vertexInfluence, orientation0);
    orientation2 = orientation0 + DELTA*orientationSensitivities;
    orientation1 = orientation0 - DELTA*orientationSensitivities;
    
    checkNaN(orientations);
    checkNaN(orientation2);
    checkNaN(orientation1);
    
//    cerr << "Orientation sensitivity " << orientationSensitivities << "\n\n";
//    cerr << "orientation0 " << orientation0 << "\n\n";
//    cerr << "orientation1 " << orientation1 << "\n\n";
//    cerr << "orientation2 " << orientation2 << "\n\n";
    
    // Obtain perturbed constitutives
    constit1 = diagonalElement(materials, fillFactors1, orientation1);
    constit2 = diagonalElement(materials, fillFactors2, orientation2);
    
//    cerr << "constit1 = " << constit1 << "\n\n";
//    cerr << "constit2 = " << constit2 << "\n\n";
    
    // Calculate the sensitivity.  This amounts to separately differencing the
    // numerators and denominators, term by term.
    Pointer<DynamicRLE3<Precision::RationalFunction> > sensitivity(
        new DynamicRLE3<Precision::RationalFunction>());
    
    RLE::transform(constit2, constit1, *sensitivity,
        Differential(2.0*DELTA), Precision::RationalFunction(0),
        Precision::RationalFunction(0));
    
    //cerr << "sensitivity = " << *sensitivity << "\n\n";
    
    return sensitivity;
}

Pointer<RLE::DynamicRLE3<Precision::RationalFunction> >
TensorGridConstants::
offDiagonalSensitivity(
    const vector<Precision::RationalFunction> & materials,
    const vector<DynamicRLE3<double> > & fillFactors,
    const DynamicRLE3<double> & orientations,
    const vector<DynamicRLE3<double> > & fillFactorSensitivities,
    const DynamicRLE3<double> & orientationSensitivities)
{
    assert(fillFactors.size() == materials.size());
    assert(fillFactorSensitivities.size() == materials.size());
    
    //const Precision::Float DELTA = sqrt(f::epsilon()); // for finite-differencing.
    Precision::Float DELTA;
    if (UserPreferences::defines("derivDelta"))
        DELTA = UserPreferences::valueAs<double>("derivDelta");
    else
        DELTA = 1e-4;
    
    // Obtain the total region that can change (it's the total support of the
    // fill factors and orientation)
    SupportRegion3 vertexInfluence;
    if (fillFactorSensitivities.size() > 0)
        vertexInfluence = SupportRegion3(fillFactorSensitivities[0]);
    for (int nn = 1; nn < fillFactorSensitivities.size(); nn++)
        vertexInfluence += SupportRegion3(fillFactorSensitivities[nn]);
    if (vertexInfluence.numRuns() > 0)
        vertexInfluence += SupportRegion3(orientationSensitivities);
    else
        vertexInfluence = SupportRegion3(orientationSensitivities);
    
    // Declare temporary arrays
    DynamicRLE3<Precision::RationalFunction> constit1, constit2;
    vector<DynamicRLE3<double> > fillFactors0(fillFactors.size()),
        fillFactors1(fillFactors.size()),
        fillFactors2(fillFactors.size());
    DynamicRLE3<double> orientation0, orientation1, orientation2;
    
    // Obtain perturbed fill factors
    for (int ff = 0; ff < fillFactors.size(); ff++)
    {
        restriction(fillFactors[ff], vertexInfluence, fillFactors0[ff]);
        fillFactors1[ff] = fillFactors0[ff] - DELTA*fillFactorSensitivities[ff];
        fillFactors2[ff] = fillFactors0[ff] + DELTA*fillFactorSensitivities[ff];
//        cerr << "fill1[" << ff << "] = " << fillFactors1[ff] << "\n\n";
//        cerr << "fill2[" << ff << "] = " << fillFactors2[ff] << "\n\n";

        checkNaN(fillFactors[ff]);
        checkNaN(fillFactors1[ff]);
        checkNaN(fillFactors2[ff]);
    }
    
    // Obtain perturbed orientations
    restriction(orientations, vertexInfluence, orientation0);
    orientation2 = orientation0 + DELTA*orientationSensitivities;
    orientation1 = orientation0 - DELTA*orientationSensitivities;
    
    checkNaN(orientations);
    checkNaN(orientation2);
    checkNaN(orientation1);
    
//    cerr << "Orientation sensitivity " << orientationSensitivities << "\n\n";
//    cerr << "orientation0 " << orientation0 << "\n\n";
//    cerr << "orientation1 " << orientation1 << "\n\n";
//    cerr << "orientation2 " << orientation2 << "\n\n";
    
    // Obtain perturbed constitutives
    //constit0 = offDiagonalElement(materials, fillFactors0, orientation0);
    constit1 = offDiagonalElement(materials, fillFactors1, orientation1);
    constit2 = offDiagonalElement(materials, fillFactors2, orientation2);
    
//    cerr << "constit0 = " << constit0 << "\n\n";
//    cerr << "constit1 = " << constit1 << "\n\n";
//    cerr << "constit2 = " << constit2 << "\n\n";
    
    // Calculate the sensitivity.  This amounts to separately differencing the
    // numerators and denominators, term by term.
    Pointer<DynamicRLE3<Precision::RationalFunction> > sensitivity(
        new DynamicRLE3<Precision::RationalFunction>());
    
    RLE::transform(constit2, constit1, *sensitivity,
        Differential(2.0*DELTA), Precision::RationalFunction(0),
        Precision::RationalFunction(0));
    
//    cerr << "sensitivity = " << *sensitivity << "\n\n";
    
    return sensitivity;
}

void TensorGridConstants::
selectFillFactorSensitivity(
    const vector<map<unsigned int, DynamicRLE3<Vector3d> > > & fillFactorJacobians,
    int vertexNumber, int xyz,
    vector<DynamicRLE3<double> > & out)
{
    out.resize(fillFactorJacobians.size());
    
    for (int ff = 0; ff < fillFactorJacobians.size(); ff++)
    {
        if (fillFactorJacobians[ff].count(vertexNumber))
        {
            RLE::transform(fillFactorJacobians[ff].find(vertexNumber)->second,
                out[ff], Vector3d::GetElement(xyz));
        }
    }
}

void TensorGridConstants::
selectOrientationSensitivity(const map<unsigned int, DynamicRLE3<Vector3d> > & 
        orientationJacobians,
    int vertexNumber, int xyz,
    DynamicRLE3<double> & out)
{
    if (orientationJacobians.count(vertexNumber))
    {
        RLE::transform(orientationJacobians.find(vertexNumber)->second,
            out, Vector3d::GetElement(xyz));
    }
}

static inline Precision::RationalFunction sScaleRationalFunction(double lhs,
    const Precision::RationalFunction & rhs)
{
    return Precision::Float(lhs)*rhs;
}

static inline Precision::RationalFunction sScaleRationalFunction2(
    const Precision::RationalFunction & lhs,
    double rhs)
{
    return lhs*Precision::Float(rhs);
}

struct HarmonicMean
{
    HarmonicMean(const vector<Precision::RationalFunction> & constitutives,
        DynamicRLE3<Precision::RationalFunction> & outHarmonicMean) :
        mReciprocals(constitutives.size()),
        mOut(outHarmonicMean)
    {
        for (int nn = 0; nn < constitutives.size(); nn++)
            mReciprocals[nn] = reciprocal(constitutives[nn]);
    }
    
    void operator()(long first, long last, double const* data,
        int const* markedArrays, int numStreams)
//    void operator()(long first, long last, Precision::Float const* data,
//        int const* markedArrays, int numStreams)
    {
        Precision::RationalFunction sum;
        for (int nn = 0; nn < numStreams; nn++)
        if (markedArrays[nn] != false)
        {
            sum += mReciprocals[nn] * Precision::Float(data[nn]);
        }
        
        mOut.RLEBase<DRLESegment<Precision::RationalFunction> >::mark(
            first, last, reciprocal(sum));
    }
    
    vector<Precision::RationalFunction> mReciprocals;
    DynamicRLE3<Precision::RationalFunction> & mOut;
};

// Interpret the last material as a background material, and use the last fill
// factor array as a mask.  Its numerical value is meaningless!
struct HarmonicMeanWithBackground
{
    HarmonicMeanWithBackground(const vector<Precision::RationalFunction> & constitutives,
        Precision::RationalFunction background,
        DynamicRLE3<Precision::RationalFunction> & outHarmonicMean) :
        mReciprocals(constitutives.size()),
        mBackground(reciprocal(background)),
        mOut(outHarmonicMean)
    {
        for (int nn = 0; nn < constitutives.size(); nn++)
            mReciprocals[nn] = reciprocal(constitutives[nn]);
    }
    
    void operator()(long first, long last, double const* data,
        int const* markedArrays, int numStreams)
    {
        Precision::RationalFunction sum;
        Precision::Float backgroundFill = 1.0;
        
        for (int nn = 0; nn < numStreams-1; nn++)
        if (markedArrays[nn] != false)
        {
            sum += mReciprocals[nn] * Precision::Float(data[nn]);
            backgroundFill -= data[nn];
        }
        
        sum += backgroundFill * mBackground;
        
        mOut.RLEBase<DRLESegment<Precision::RationalFunction> >::mark(
            first, last, reciprocal(sum));
    }
    
    vector<Precision::RationalFunction> mReciprocals;
    Precision::RationalFunction mBackground;
    DynamicRLE3<Precision::RationalFunction> & mOut;
};


void TensorGridConstants::
calculateHarmonicMean(
    const vector<DynamicRLE3<double> > & yeeFillFactors,
    vector<Precision::RationalFunction> inverseConstits,
    DynamicRLE3<Precision::RationalFunction> & outHarmonicMean)
{
    // HEY.  The reason I can't give the permittivity the "right" dimensions is
    // that the fill factors and orientations start out 3D and are downsampled.
    // This screws everything up.
//    outHarmonicMean.dimensions(voxels().numNonSingularDims()
    //if (yeeFillFactors.size() > 0)
    //    outHarmonicMean.dimensions(yeeFillFactors.begin()->capacity());
    //else
//    outHarmonicMean.dimensions(voxels().size(0), voxels().size(1), voxels().size(2));
    
    vector<DynamicRLE3<double> const*> fillFactors;
    for (int nn = 0; nn < yeeFillFactors.size(); nn++)
        fillFactors.push_back(&yeeFillFactors[nn]);
    
    if (backgroundMaterial() != Precision::Float(0.0))
    {
        DynamicRLE3<double> thisSignifiesBackground(voxels().num().asArray());
//        for (int ff = 0; ff < yeeFillFactors.size(); ff++)
//        {
//            thisSignifiesBackground = thisSignifiesBackground + 
//                (0.0*yeeFillFactors[ff] + 1.0);
//        }
        thisSignifiesBackground.mark(
            voxels().p1[0], voxels().p1[1], voxels().p1[2],
            voxels().p2[0], voxels().p2[1], voxels().p2[2], 1.0);
        
        fillFactors.push_back(&thisSignifiesBackground);
        inverseConstits.push_back(backgroundMaterial());
        
        HarmonicMeanWithBackground callback(inverseConstits,
            reciprocal(backgroundMaterial()), outHarmonicMean);
        RLE::multiApply(&(fillFactors[0]), fillFactors.size(), callback,
            RLE::ExcludeGaps);
    }
    else
    {
        HarmonicMean callback(inverseConstits, outHarmonicMean);
        RLE::multiApply(&(fillFactors[0]), fillFactors.size(), callback,
            RLE::ExcludeGaps);
    }
}

struct ArithmeticMean
{
    ArithmeticMean(const vector<Precision::RationalFunction> & constitutives,
        DynamicRLE3<Precision::RationalFunction> & outArithmeticMean) :
        mConstitutives(constitutives),
        mOut(outArithmeticMean)
    {
    }
    
//    void operator()(long first, long last, Precision::Float const* data,
//        int const* markedArrays, int numStreams)
    void operator()(long first, long last, double const* data,
        int const* markedArrays, int numStreams)
    {
        Precision::RationalFunction sum;
        for (int nn = 0; nn < numStreams; nn++)
        if (markedArrays[nn] != false)
        {
            sum += mConstitutives[nn] * Precision::Float(data[nn]);
        }
        
        mOut.RLEBase<DRLESegment<Precision::RationalFunction> >::mark(
            first, last, sum);
    }
    
    const vector<Precision::RationalFunction> & mConstitutives;
    DynamicRLE3<Precision::RationalFunction> & mOut;
};

struct ArithmeticMeanWithBackground
{
    ArithmeticMeanWithBackground(const vector<Precision::RationalFunction> & constitutives,
        Precision::RationalFunction background,
        DynamicRLE3<Precision::RationalFunction> & outArithmeticMean) :
        mConstitutives(constitutives),
        mBackground(background),
        mOut(outArithmeticMean)
    {
    }
    
    void operator()(long first, long last, double const* data,
        int const* markedArrays, int numStreams)
    {
        Precision::RationalFunction sum;
        Precision::Float backgroundFill = 1.0;
        
//        std::cout << "==== NEW:\n";
//        std::cout << sum << "\n";
        for (int nn = 0; nn < numStreams-1; nn++)
        if (markedArrays[nn] != false)
        {
//            std::cout << "Add " << mConstitutives[nn] * Precision::Float(data[nn]) << "\n";
            sum += mConstitutives[nn] * Precision::Float(data[nn]);
//            std::cout << sum << "\n";
            backgroundFill -= data[nn];
        }
//        std::cout << "Add " << backgroundFill * mBackground << "\n";
        sum += backgroundFill * mBackground;
//        std::cout << sum << "\n";
        
//        if (sum.numerator().order() != sum.denominator().order())
//        {
//            std::cout << backgroundFill*mBackground << " + ";
//            for (int nn = 0; nn < numStreams-1; nn++)
//            if (markedArrays[nn] != false)
//            {
//                std::cout << data[nn]*mConstitutives[nn];
//                if (nn < numStreams-2)
//                {
//                    std::cout << " + ";
//                }
//            }
//            std::cout << " = " << sum << "\n";
//        }
        
        mOut.RLEBase<DRLESegment<Precision::RationalFunction> >::mark(first, last, sum);
    }
    
    const vector<Precision::RationalFunction> & mConstitutives;
    Precision::RationalFunction mBackground;
    DynamicRLE3<Precision::RationalFunction> & mOut;
};

void TensorGridConstants::
calculateArithmeticMean(
    const vector<DynamicRLE3<double> > & yeeFillFactors,
    vector<Precision::RationalFunction> inverseConstits,
    DynamicRLE3<Precision::RationalFunction> & outArithmeticMean)
{
//    outArithmeticMean.dimensions(voxels().num().asArray());
//    outArithmeticMean.dimensions(true,true,true);
//    if (yeeFillFactors.size() > 0)
//        outArithmeticMean.dimensions(yeeFillFactors.begin()->capacity());
//    else
//    outArithmeticMean.dimensions(voxels().size(0), voxels().size(1), voxels().size(2));
    
    vector<DynamicRLE3<double> const*> fillFactors;
    for (int nn = 0; nn < yeeFillFactors.size(); nn++)
        fillFactors.push_back(&yeeFillFactors[nn]);
    
    if (backgroundMaterial() != Precision::Float(0.0))
    {
        DynamicRLE3<double> thisSignifiesBackground(voxels().num().asArray());
//        for (int ff = 0; ff < yeeFillFactors.size(); ff++)
//        {
//            thisSignifiesBackground = thisSignifiesBackground + 
//                (0.0*yeeFillFactors[ff] + 1.0);
//        }
        thisSignifiesBackground.mark(
            voxels().p1[0], voxels().p1[1], voxels().p1[2],
            voxels().p2[0], voxels().p2[1], voxels().p2[2], 1.0);
        
        fillFactors.push_back(&thisSignifiesBackground);
        inverseConstits.push_back(backgroundMaterial());
        
        ArithmeticMeanWithBackground callback(inverseConstits,
            reciprocal(backgroundMaterial()), outArithmeticMean);
        RLE::multiApply(&(fillFactors[0]), fillFactors.size(), callback,
            RLE::ExcludeGaps);
    }
    else
    {
        ArithmeticMean callback(inverseConstits, outArithmeticMean);
        RLE::multiApply(&(fillFactors[0]), fillFactors.size(), callback,
            RLE::ExcludeGaps);
    }
}

void TensorGridConstants::
writeSensitivity(
    Pointer<RLE::DynamicRLE3<Precision::RationalFunction> > sensitivity,
    const RLE::SupportRegion3 & savedFields,
    std::ostream & file) const
{
    if (false == savedFields.encloses(SupportRegion3(*sensitivity)))
    {
//        cerr << "WARNING: savedFields does not enclose all sensitive cells.\n";
//        cerr << "(This may occur when filtering boundary outputs to include "
//            "only axis-aligned boundaries.)\n";
        
//        cerr << "----- Saved:\n";
//        cerr << savedFields << "\n";
//        cerr << "----- Sensitive:\n";
//        cerr << SupportRegion3(*sensitivity) << "\n";
    }
    IndexArray3 savedIndices(savedFields);
    
    const uint32_t MAXORDER = 10;
    
    const int ANYORDER = -1;
    for (uint32_t numerOrder = 0; numerOrder < MAXORDER; numerOrder++)
    {
        DynamicRLE3<Precision::RationalFunction> sensitivityAtOrder;
        DynamicRLE3<Precision::RationalFunction> selectedSensitivityAtOrder;
        filter(*sensitivity, sensitivityAtOrder,
            OrderGreaterThan(numerOrder, ANYORDER));
        
        IndexArray3 indices;
        restriction(savedIndices, sensitivityAtOrder, indices);
        restriction(sensitivityAtOrder, indices, selectedSensitivityAtOrder);
        sensitivityAtOrder = selectedSensitivityAtOrder;
        
        if (sensitivityAtOrder.numRuns() == 0)
            continue;
            
        file.write("NUME", 4);
        file.write((char*)&numerOrder, sizeof(uint32_t));
        uint32_t len = sensitivityAtOrder.length();
        file.write((char*) &len, sizeof(uint32_t));
        
        DynamicRLE3<Precision::RationalFunction>::ConstIterator itr;
        for (itr = sensitivityAtOrder.begin(); itr != sensitivityAtOrder.end();
            itr.nextMarkedRun())
        {
            for (long ll = itr.runStart(); ll <= itr.runEnd(); ll++)
            {
                assert(!isnan(itr->markAt(ll).numerator()[numerOrder]));
                double val = itr->markAt(ll).numerator()[numerOrder];
                file.write((char*) &val, sizeof(double));
            }
        }
        
        IndexArray3::ConstIterator itr2;
        for (itr2 = indices.begin(); itr2 != indices.end();
            itr2.nextMarkedRun())
        {
            for (long ll = itr2.runStart(); ll <= itr2.runEnd();
                ll++)
            {
                uint32_t ind = itr2->markAt(ll)+1;
                file.write((char*) &ind, sizeof(uint32_t));
            }
        }
    }
    file.write("NEND", 4);
    
    for (uint32_t denomOrder = 0; denomOrder < MAXORDER; denomOrder++)
    {
        DynamicRLE3<Precision::RationalFunction> sensitivityAtOrder;
        DynamicRLE3<Precision::RationalFunction> selectedSensitivityAtOrder;
        
        filter(*sensitivity, sensitivityAtOrder,
            OrderGreaterThan(ANYORDER, denomOrder));
        
        IndexArray3 indices;
        restriction(savedIndices, sensitivityAtOrder, indices);
        restriction(sensitivityAtOrder, indices, selectedSensitivityAtOrder);
        sensitivityAtOrder = selectedSensitivityAtOrder;
        
        if (sensitivityAtOrder.numRuns() == 0)
            continue;
            
        file.write("DENO", 4);
        file.write((char*)&denomOrder, sizeof(uint32_t));
        uint32_t len = sensitivityAtOrder.length();
        file.write((char*) &len, sizeof(uint32_t));
        
        DynamicRLE3<Precision::RationalFunction>::ConstIterator itr;
        for (itr = sensitivityAtOrder.begin(); itr != sensitivityAtOrder.end();
            itr.nextMarkedRun())
        {
            for (long ll = itr.runStart(); ll <= itr.runEnd();
                ll++)
            {
                assert(!isnan(itr->markAt(ll).denominator()[denomOrder]));
                double val = itr->markAt(ll).denominator()[denomOrder];
                file.write((char*) &val, sizeof(double));
            }
        }
        
        IndexArray3::ConstIterator itr2;
        for (itr2 = indices.begin(); itr2 != indices.end();
            itr2.nextMarkedRun())
        {
            for (long ll = itr2.runStart(); ll <= itr2.runEnd();
                ll++)
            {
                uint32_t ind = itr2->markAt(ll)+1;
                file.write((char*) &ind, sizeof(uint32_t));
            }
        }
    }
    file.write("DEND", 4);
}

