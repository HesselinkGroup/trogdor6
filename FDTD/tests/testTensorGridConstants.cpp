/*
 *  testOrientedVoxels.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/29/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestVoxelizer

#include <boost/test/unit_test.hpp>

#include "TensorGridConstants.h"

#include <iomanip>

using namespace std;
using namespace RLE;
using namespace YeeUtilities;

BOOST_AUTO_TEST_CASE( EffectivePermittivity )
{
    double tolerance;
    if (sizeof(Precision::Float) == sizeof(float))
        tolerance = 1e-5;
    else
        tolerance = 1e-9;
    vector<Precision::RationalFunction> materials(2);
    materials[0] = 1.0;
    materials[1] = 2.0;
    
    vector<DynamicRLE3<double> > fillFactors(2);
    fillFactors[0].mark(0,0,0,0.5);
    fillFactors[1].mark(0,0,0,0.5);
    
    DynamicRLE3<Precision::RationalFunction> inversePermittivity;
    
    // Parameters just to fill space in the TensorGridConstants constructor
    const int WHATEVER = 0;
    DynamicRLE3<Precision::RationalFunction> WHO_CARES;
    TensorGridConstants tgc(WHATEVER, WHATEVER, WHO_CARES, Rect3i(0,0,0,0,0,0));
    
    // Test normal-to-surface diagonal component:
    //  permittivity is harmonic mean
    //  permittivity is of the form 1/a0, so the numerator of 1/eps is a0,
    //   or 1/(harmonic mean).  The harmonic mean of 1 and 2 is 4/3, so a0 is 
    //   0.75.
    {
        DynamicRLE3<double> orientation;
        orientation.mark(0,0,0,1.0);
        inversePermittivity = tgc.diagonalElement(materials, fillFactors, orientation);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).numerator().order(), 0);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).denominator().order(), 0);
        
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).numerator()[0], 0.75, tolerance);
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).denominator()[0], 1.0, tolerance);
    }
    
    // Test parallel-to-surface diagonal component:
    //  permittivity is arithmetic mean
    //  1/a0 = 1.5, so a0 = 2/3
    {
        DynamicRLE3<double> orientation;
        orientation.mark(0,0,0,0.0);
        inversePermittivity = tgc.diagonalElement(materials, fillFactors, orientation);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).numerator().order(), 0);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).denominator().order(), 0);
        
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).numerator()[0], 2.0/3, tolerance);
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).denominator()[0], 1.0, tolerance);
    }
    
    // Test off-diagonal component for diagonal orientation tensor:
    // should be zero.
    {
        DynamicRLE3<double> orientation;
        orientation.mark(0,0,0,0.0);
        inversePermittivity = tgc.offDiagonalElement(materials, fillFactors, orientation);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).numerator().order(), 0);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).denominator().order(), 0);
        
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).numerator()[0], 0.0, tolerance);
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).denominator()[0], 1.0, tolerance);
    }
    
    // Test a nonzero off-diagonal of the tensor.
    {
        double harmonicMean = 4.0/3;
        double arithmeticMean = 1.5;
        double orientation_ij = 0.5;
        
        double epsInv = orientation_ij*(1.0/harmonicMean - 1.0/arithmeticMean);
        
        DynamicRLE3<double> orientation;
        orientation.mark(0,0,0, orientation_ij);
        inversePermittivity = tgc.offDiagonalElement(materials,
            fillFactors, orientation);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).numerator().order(), 0);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).denominator().order(), 0);
        
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).numerator()[0], epsInv, tolerance);
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).denominator()[0], 1.0, tolerance);
    }
}

BOOST_AUTO_TEST_CASE( BackgroundMaterial )
{
    double tolerance;
    if (sizeof(Precision::Float) == sizeof(float))
        tolerance = 1e-5;
    else
        tolerance = 1e-9;
    vector<Precision::RationalFunction> materials(1);
    materials[0] = 1.0;
//    materials[1] = 2.0;
    
    vector<DynamicRLE3<double> > fillFactors(1);
    fillFactors[0].mark(0,0,0,0.5);
//    fillFactors[1].mark(0,0,0,0.5);
    
    DynamicRLE3<Precision::RationalFunction> inversePermittivity;
    
    // Parameters just to fill space in the TensorGridConstants constructor
    const int WHATEVER = 0;
    DynamicRLE3<Precision::RationalFunction> WHO_CARES;
    TensorGridConstants tgc(WHATEVER, WHATEVER, WHO_CARES,
        Rect3i(0,0,0,0,0,0),
        Precision::RationalFunction(2.0));
    
    // Test normal-to-surface diagonal component:
    //  permittivity is harmonic mean
    //  permittivity is of the form 1/a0, so the numerator of 1/eps is a0,
    //   or 1/(harmonic mean).  The harmonic mean of 1 and 2 is 4/3, so a0 is 
    //   0.75.
    {
        DynamicRLE3<double> orientation;
        orientation.mark(0,0,0,1.0);
        inversePermittivity = tgc.diagonalElement(materials,
            fillFactors, orientation);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).numerator().order(), 0);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).denominator().order(), 0);
        
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).numerator()[0], 0.75, tolerance);
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).denominator()[0], 1.0, tolerance);
    }
    
    // Test parallel-to-surface diagonal component:
    //  permittivity is arithmetic mean
    //  1/a0 = 1.5, so a0 = 2/3
    {
        DynamicRLE3<double> orientation;
        orientation.mark(0,0,0,0.0);
        inversePermittivity = tgc.diagonalElement(materials,
            fillFactors, orientation);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).numerator().order(), 0);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).denominator().order(), 0);
        
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).numerator()[0], 2.0/3, tolerance);
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).denominator()[0], 1.0, tolerance);
    }
    
    // Test off-diagonal component for diagonal orientation tensor:
    // should be zero.
    {
        DynamicRLE3<double> orientation;
        orientation.mark(0,0,0,0.0);
        inversePermittivity = tgc.offDiagonalElement(materials, fillFactors, orientation);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).numerator().order(), 0);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).denominator().order(), 0);
        
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).numerator()[0], 0.0, tolerance);
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).denominator()[0], 1.0, tolerance);
    }
    
    // Test a nonzero off-diagonal of the tensor.
    {
        double harmonicMean = 4.0/3;
        double arithmeticMean = 1.5;
        double orientation_ij = 0.5;
        
        double epsInv = orientation_ij*(1.0/harmonicMean - 1.0/arithmeticMean);
        
        DynamicRLE3<double> orientation;
        orientation.mark(0,0,0, orientation_ij);
        inversePermittivity = tgc.offDiagonalElement(materials, fillFactors, orientation);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).numerator().order(), 0);
        BOOST_CHECK_EQUAL(inversePermittivity.at(0,0,0).denominator().order(), 0);
        
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).numerator()[0], epsInv, tolerance);
        BOOST_CHECK_CLOSE(inversePermittivity.at(0,0,0).denominator()[0], 1.0, tolerance);
    }
}

BOOST_AUTO_TEST_CASE( Sensitivity_FillFactor_Analytical )
{
    double eps1 = 1.0, eps2 = 2.0;
    double fill1 = 0.5, fill2 = 1.0 - fill1;
    vector<Precision::RationalFunction> materials(2);
    materials[0] = eps1;
    materials[1] = eps2;
    
    vector<DynamicRLE3<double> > fillFactors(2);
    fillFactors[0].mark(0,0,0, fill1);
    fillFactors[1].mark(0,0,0, fill2);
    
    DynamicRLE3<Precision::RationalFunction> epsInv;
    DynamicRLE3<Precision::RationalFunction> DepsInv;
    
    vector<DynamicRLE3<double> > dFillFactors(2);
    dFillFactors[0].mark(0,0,0, 1.0);
    dFillFactors[1].mark(0,0,0, -1.0);
    DynamicRLE3<double> dOrientation; // all zeros!
    
    const int WHATEVER = 0;
    DynamicRLE3<Precision::RationalFunction> WHO_CARES;
    TensorGridConstants tgc(WHATEVER, WHATEVER, WHO_CARES, Rect3i(0,0,0,0,0,0));

    
    // Test normal-to-surface diagonal component:
    //  permittivity is harmonic mean
    //  permittivity is of the form 1/a0, so the numerator of 1/eps is a0,
    //   or 1/(harmonic mean).  The harmonic mean of 1 and 2 is 4/3, so a0 is 
    //   0.75.
    {
        DynamicRLE3<double> orientation;
        orientation.mark(0,0,0,1.0);
        epsInv = tgc.diagonalElement(materials,
            fillFactors, orientation);
        DepsInv = *tgc.diagonalSensitivity(materials,
            fillFactors, orientation, dFillFactors, dOrientation);
        double Da0 = 1/eps1 - 1/eps2; // derive it yerself!
        
        BOOST_CHECK_CLOSE(DepsInv.at(0,0,0).numerator()[0], Da0, 1); // to 1%
    }
    
    // Test parallel-to-surface diagonal component:
    //  permittivity is arithmetic mean
    //  1/a0 = 1.5, so a0 = 2/3
    {
        DynamicRLE3<double> orientation;
        orientation.mark(0,0,0,0.0);
        epsInv = tgc.diagonalElement(materials,
            fillFactors, orientation);
        DepsInv = *tgc.diagonalSensitivity(materials,
            fillFactors, orientation, dFillFactors, dOrientation);
        double meanEps = fill1*eps1 + fill2*eps2;
        double Da0 = (eps2-eps1)/meanEps/meanEps; // derive it yerself!
        
        BOOST_CHECK_CLOSE(DepsInv.at(0,0,0).numerator()[0], Da0, 1); // to 1%
    }
}



// Test sensitivity to variation of fill factors.
BOOST_AUTO_TEST_CASE( Sensitivity_FillFactor_FiniteDifference )
{
    vector<Precision::RationalFunction> materials(2);
    materials[0] = 1.0;
    materials[1] = Precision::RationalFunction(
        Precision::Polynomial(10.0, -20.0, 30.0),
        Precision::Polynomial(11.0, -22.0, 33.0) );
        
    const int WHATEVER = 0;
    DynamicRLE3<Precision::RationalFunction> WHO_CARES;
    TensorGridConstants tgc(WHATEVER, WHATEVER, WHO_CARES, Rect3i(0,0,0,0,0,0));
    
//    typedef numeric_limits<Precision::Float> f;
    const double DELTA = 0.01;
    //const double DELTA = sqrt(f::epsilon());
    const double TOLERANCE = 10*DELTA;
//    cerr << "Using DELTA = " << DELTA << ", TOLERANCE = " << TOLERANCE << "\n";
    
    vector<DynamicRLE3<double> > fillFactors1(2), fillFactors2(2);
    fillFactors1[0].mark(0,0,0,0.5);
    fillFactors1[1].mark(0,0,0,0.5);
    fillFactors2[0] = fillFactors1[0] + DELTA;
    fillFactors2[1] = fillFactors1[1] - DELTA;
    vector<DynamicRLE3<double> > dFillFactors(2);
    dFillFactors[0] = (fillFactors2[0] - fillFactors1[0])/DELTA/2;
    dFillFactors[1] = (fillFactors2[1] - fillFactors1[1])/DELTA/2;
    
    DynamicRLE3<double> orientation;
    DynamicRLE3<double> dOrientation; // empty
    orientation.mark(0,0,0, 0.4); // some arbitrary value.
    
    DynamicRLE3<Precision::RationalFunction> epsInv1, epsInv2;
    DynamicRLE3<Precision::RationalFunction> epsInvSensitivity;
    
    epsInv1 = tgc.diagonalElement(materials, fillFactors1, orientation);
    epsInv2 = tgc.diagonalElement(materials, fillFactors2, orientation);
    epsInvSensitivity = *tgc.diagonalSensitivity(materials, fillFactors1, orientation, dFillFactors, dOrientation);
    
    TensorGridConstants::Differential diff(2*DELTA);
    
    if (0)
    {
        cerr << "eps1:\n" << epsInv1 << "\n";
        cerr << "eps2:\n" << epsInv2 << "\n";
        cerr << "eps sensitivity: " << epsInvSensitivity << "\n";
        cerr << "by differential:\n" << diff(epsInv2.at(0,0,0), epsInv1.at(0,0,0))
            << "\n";
    }
    
    double deltaNumer = epsInv2.at(0,0,0).numerator()[0] -
        epsInv1.at(0,0,0).numerator()[0];
    double deltaDenom = epsInv2.at(0,0,0).denominator()[0] -
        epsInv1.at(0,0,0).denominator()[0];
    double dNumer = epsInvSensitivity.at(0,0,0).numerator()[0];
    double dDenom = epsInvSensitivity.at(0,0,0).denominator()[0];
    
    BOOST_CHECK_SMALL(dDenom, TOLERANCE*100); // technically zero
    BOOST_CHECK_SMALL(deltaDenom, TOLERANCE*100); // technically zero
    BOOST_CHECK_CLOSE(dNumer, deltaNumer/DELTA/2, 1.0); // to within 1% is ok
}



BOOST_AUTO_TEST_CASE( Sensitivity_Orientation_FiniteDifference )
{
    const int WHATEVER = 0;
    DynamicRLE3<Precision::RationalFunction> WHO_CARES;
    TensorGridConstants tgc(WHATEVER, WHATEVER, WHO_CARES, Rect3i(0,0,0,0,0,0));
    
    vector<Precision::RationalFunction> materials(2);
    materials[0] = 1.0;
    materials[1] = 2.0;
    
    typedef numeric_limits<Precision::Float> f;
    const double DELTA = sqrt(f::epsilon());
    const double TOLERANCE = 10*DELTA;
    BOOST_TEST_MESSAGE("Using DELTA = " << DELTA << ", TOLERANCE = " << TOLERANCE);
    
    vector<DynamicRLE3<double> > fillFactors(2);
    vector<DynamicRLE3<double> > deltaFillFactors(2);
    fillFactors[0].mark(0,0,0,0.5);
    fillFactors[1].mark(0,0,0,0.5);
    
    DynamicRLE3<double> orientation1, orientation2;
    orientation1.mark(0,0,0, 0.4); // some arbitrary value.
    orientation2 = orientation1 + DELTA;
    DynamicRLE3<double> dOrientation = (orientation2 - orientation1)/DELTA;
    
    DynamicRLE3<Precision::RationalFunction> epsInv1, epsInv2;
    DynamicRLE3<Precision::RationalFunction> epsInvSensitivity;
    
    epsInv1 = tgc.diagonalElement(materials, fillFactors,
        orientation1);
    epsInv2 = tgc.diagonalElement(materials, fillFactors,
        orientation2);
    epsInvSensitivity = *tgc.diagonalSensitivity(materials,
        fillFactors, orientation1, deltaFillFactors, dOrientation);
    
//    cout << "eps1:\n" << epsInv1 << "\n";
//    cout << "eps2:\n" << epsInv2 << "\n";
//    cout << "eps sensitivity: " << epsInvSensitivity << "\n";
//    cout << "diff:\n" << epsInv2 - epsInv1 << "\n";
    
    double deltaNumer = epsInv2.at(0,0,0).numerator()[0] -
        epsInv1.at(0,0,0).numerator()[0];
    double deltaDenom = epsInv2.at(0,0,0).denominator()[0] -
        epsInv1.at(0,0,0).denominator()[0];
    double dNumer = epsInvSensitivity.at(0,0,0).numerator()[0];
    double dDenom = epsInvSensitivity.at(0,0,0).denominator()[0];
    
    BOOST_CHECK_SMALL(dDenom, TOLERANCE*100); // technically zero; tolerance %
    BOOST_CHECK_SMALL(deltaDenom, TOLERANCE*100); // technically zero
    BOOST_CHECK_CLOSE(dNumer, deltaNumer/DELTA, TOLERANCE*100); // within 0.01%
}







