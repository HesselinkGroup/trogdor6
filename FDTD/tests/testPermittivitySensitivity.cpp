/*
 *  testPermittivitySensitivity.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 9/11/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestPermittivitySensitivity

#include <boost/test/unit_test.hpp>

#include "TensorGridConstants.h"

#include <iomanip>

using namespace std;
using namespace RLE;
using namespace YeeUtilities;

Matrix3d sBigDot(const Matrix3<Vector3d> & lhs, const Vector3d & rhs)
{
    return Matrix3d(
        dot(lhs[0], rhs), dot(lhs[1], rhs), dot(lhs[2], rhs),
        dot(lhs[3], rhs), dot(lhs[4], rhs), dot(lhs[5], rhs),
        dot(lhs[6], rhs), dot(lhs[7], rhs), dot(lhs[8], rhs));
}

BOOST_AUTO_TEST_CASE( TestDiagonalSensitivity )
{
    const int WHATEVER = 0;
    DynamicRLE3<Precision::RationalFunction> WHO_CARES;
    TensorGridConstants tgc(WHATEVER, WHATEVER, WHO_CARES, Rect3i(0,0,0,0,0,0));
    
    vector<Precision::RationalFunction> materials;
    materials.push_back(Precision::RationalFunction(1.0));
    materials.push_back(Precision::RationalFunction(2.0));
    
    vector<DynamicRLE3<double> > fillFactors0(2, DynamicRLE3<double>(0,0,0)),
        fillFactors1(2, DynamicRLE3<double>(0,0,0));
    DynamicRLE3<double> orientations0(0,0,0), orientations1(0,0,0);
    vector<DynamicRLE3<double> > fillFactorSensitivities(2, DynamicRLE3<double>(0,0,0));
    DynamicRLE3<double> orientationSensitivities(0,0,0);
    
    fillFactors0[0].mark(0,0,0,0.75);
    fillFactors0[1].mark(0,0,0,0.25);
    fillFactors1[0].mark(0,0,0,0.76);
    fillFactors1[1].mark(0,0,0,0.24);
    
    orientations0.mark(0,0,0,0.5);
    orientations1.mark(0,0,0,0.6);
    
    fillFactorSensitivities[0].mark(0,0,0,0.01); // durr.  just to be sure.
    fillFactorSensitivities[1].mark(0,0,0,-0.01);
    
    orientationSensitivities.mark(0,0,0,0.1);
    
    Pointer<DynamicRLE3<Precision::RationalFunction> > dPermittivity =
        tgc.diagonalSensitivity(materials, fillFactors0,
        orientations0, fillFactorSensitivities, orientationSensitivities);
    
    Precision::RationalFunction xsi0 =
        tgc.diagonalElement(materials, fillFactors0,
        orientations0).at(0,0,0);
    Precision::RationalFunction xsi1 =
        tgc.diagonalElement(materials, fillFactors1,
        orientations1).at(0,0,0);
    Precision::RationalFunction dXsi = dPermittivity->at(0,0,0);
    
    if (0)
    {
        cerr << "Mat 1: " << xsi0 << "\nMat 2: " << xsi1 << "\n";
        cerr << "Deriv: " << dXsi << "\n";
        cerr << "DeltNum: " << xsi1.numerator()[0] - xsi0.numerator()[0] << "\n";
        cerr << "DerivNum: " << dXsi.numerator()[0] << "\n";
    }
    
    BOOST_CHECK(fabs(xsi1.numerator()[0] - xsi0.numerator()[0] - 
        dXsi.numerator()[0]) < 0.001);
}


