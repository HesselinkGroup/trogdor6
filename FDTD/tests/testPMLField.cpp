/*
 *  testPMLField.h
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 2/12/12.
 *  Copyright 2012 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestPMLField

#include <boost/test/unit_test.hpp>
#include <cmath>

#include "PMLField.h"
#include "PhysicalConstants.h"

using namespace std;

// Test the static function to return depth into the PML
BOOST_AUTO_TEST_CASE( TestPMLDepth_Left )
{
    Rect3i yeeCells(0, 0, 0, 99, 0, 0);
    Rect3i nonPMLHalfCells(20, 0, 0, 199, 1, 1);
    PhysicalExtents extents(yeeCells, nonPMLHalfCells);
    BOOST_CHECK(extents.hasPML(0));
    BOOST_CHECK(false == extents.hasPML(1));
    
    RLE::DynamicRLE3<double> depth = pmlDepth(extents, 0);
    
    BOOST_CHECK(depth.hasDimension(0));
    BOOST_CHECK(false == depth.hasDimension(1));
    BOOST_CHECK(false == depth.hasDimension(2));
    
    Rect3i leftPML = extents.pmlHalfCells(0);
    BOOST_CHECK_EQUAL(leftPML.count(), 4*20);
    
    for (int nn = 0; nn < 20; nn++)
    {
        double d = (20-nn)/20.0;
        BOOST_CHECK_CLOSE(depth.at(nn,0,0), d, 1e-6);
    }
}

// The way PML works is:
// - the outermost non-PML cell has depth = 0
// - the outermost PML cell has depth = 1
BOOST_AUTO_TEST_CASE( TestPMLDepth_Both )
{
    Rect3i yeeCells(0, 0, 0, 99, 0, 0);
    Rect3i nonPMLHalfCells(20, 0, 0, 179, 1, 1);
    PhysicalExtents extents(yeeCells, nonPMLHalfCells);
    BOOST_CHECK(extents.hasPML(0));
    BOOST_CHECK(extents.hasPML(1));
    
    RLE::DynamicRLE3<double> depth = pmlDepth(extents, 0);
    
    BOOST_CHECK(depth.hasDimension(0));
    BOOST_CHECK(false == depth.hasDimension(1));
    BOOST_CHECK(false == depth.hasDimension(2));
    
    Rect3i leftPML = extents.pmlHalfCells(0);
    BOOST_CHECK_EQUAL(leftPML.count(), 4*20);
    Rect3i rightPML = extents.pmlHalfCells(1);
    BOOST_CHECK_EQUAL(rightPML.count(), 4*20);
    
//    cerr << "PML rects " << leftPML << " and " << rightPML << "\n";
    
    // Test the LEFT side depths
    for (int nn = 0; nn < 20; nn++)
    {
        double d = (20-nn)/20.0;
        BOOST_CHECK_CLOSE(depth.at(nn,0,0), d, 1e-6);
    }
    // In particular:
    BOOST_CHECK_CLOSE(depth.at(leftPML.p1[0],0,0), 1.0, 1e-6);
    BOOST_CHECK_CLOSE(depth.at(leftPML.p2[0],0,0), 1.0/20.0, 1e-6);
    
    // Test the RIGHT side depths
    for (int nn = 180; nn < 200; nn++)
    {
        double d = (nn-179)/20.0;
        BOOST_CHECK_CLOSE(depth.at(nn,0,0), d, 1e-6);
    }
    for (int nn = 20; nn < 180; nn++)
        BOOST_CHECK_EQUAL(depth.at(nn,0,0), 0.0);
    // In particular:
    BOOST_CHECK_CLOSE(depth.at(rightPML.p1[0],0,0), 1.0/20.0, 1e-6);
    BOOST_CHECK_CLOSE(depth.at(rightPML.p2[0],0,0), 1.0, 1e-6);
}


BOOST_AUTO_TEST_CASE( TestSigma )
{
    Precision::Vec3 dxyz(1, 1, 1);
    Rect3i yeeCells(0, 0, 0, 99, 0, 0);
    Rect3i nonPMLHalfCells(20, 0, 0, 179, 1, 1);
    PhysicalExtents extents(yeeCells, nonPMLHalfCells);
    BOOST_CHECK(extents.hasPML(0));
    BOOST_CHECK(extents.hasPML(1));
    
    Rect3i leftPML = extents.pmlHalfCells(0);
    BOOST_CHECK_EQUAL(leftPML.count(), 4*20);
    Rect3i rightPML = extents.pmlHalfCells(1);
    BOOST_CHECK_EQUAL(rightPML.count(), 4*20);
    
//    cerr << "PML rects " << leftPML << " and " << rightPML << "\n";
    
    PMLField eField(YeeUtilities::Octant::electric());
    
    eField.init(extents, dxyz, Map<std::string, std::string>());
    
    const RLE::DynamicRLE3<PMLParameters>* xParams = eField.along(0);
    
    double eps0 = Constants::eps0;
    double mu0 = Constants::mu0;
    double c = 1.0/sqrt(eps0*mu0);
    
    double gridDiagonal = norm(dxyz*extents.physicalYeeCells().size());
    
    // Test the left side PML.
    Rect3i pmlYee;
    
    // Ex field: no attenuation should be carried out, but I guess I still 
    // calculate the PML parameters.  Seems a little lazy/misleading, but as
    // long as the actual PML implementation doesn't do anything to Ex it's
    // probably ok, right?
    
    pmlYee = YeeUtilities::halfToYee(leftPML, YeeUtilities::octantE(0));
    
    for (int ii = pmlYee.p1[0]; ii <= pmlYee.p2[0]; ii++)
    {
        double x = (ii+0.5)*dxyz[0];
        double d = (20*dxyz[0]*0.5-x)/(20.0*dxyz[0]*0.5);
        
        double sigma = d*d*d*0.8*4/sqrt(mu0/eps0)/dxyz[0];
        BOOST_CHECK_CLOSE(xParams[0].at(ii,0,0).sigma(), sigma, 1e-6);
        
        double kappa = 1 + (5-1)*d*d*d;
        BOOST_CHECK_CLOSE(xParams[0].at(ii,0,0).kappa(), kappa, 1e-6);
        
        double alpha = (1-d)*c*eps0/gridDiagonal;
        BOOST_CHECK_CLOSE(xParams[0].at(ii,0,0).alpha(), alpha, 1e-6);
        
//        cerr << "k(" << x << ") = " << kappa << "\n";
    }
    
    // Ey field: 
    
    for (int fieldXYZ = 0; fieldXYZ < 3; fieldXYZ++)
    {
        Rect3i r;
    }
}



BOOST_AUTO_TEST_CASE( Test1 )
{
    float dx = 5e-9;
    float gridDiagonalSize = 1000e-9;
    float greaterThanOne = 1.001;
    
    float sigmaBounds[] = { 0.0, greaterThanOne*0.8f*4/float(sqrt(Constants::eta0))/dx };
    float alphaBounds[] = { 0.0, greaterThanOne*float(Constants::c)*float(Constants::eps0)/50/dx };
    float kappaBounds[] = { 1.0, greaterThanOne*5.0f };
    
    BOOST_TEST_MESSAGE("Cell size " << dx << " meters.");
    BOOST_TEST_MESSAGE("Sigma bounds [" << sigmaBounds[0] << ", " <<
        sigmaBounds[1] << "].");
    BOOST_TEST_MESSAGE("Kappa bounds [" << kappaBounds[0] << ", " <<
        kappaBounds[1] << "].");
    BOOST_TEST_MESSAGE("Alpha bounds [" << alphaBounds[0] << ", " <<
        alphaBounds[1] << "].");
    
    CalculatePMLParameters callback(dx, gridDiagonalSize,
        Constants::eps0, Constants::mu0);
    
    for (double depth = -1.0; depth <= 2.0; depth += 0.1)
    {
        PMLParameters p = callback(depth);
        
        // depths out of range should yield default PML parameters
        if (depth <= 0.0 || depth > 1.0)
        {
            BOOST_CHECK(p == PMLParameters());
        }
        else
        {
            BOOST_CHECK(p.sigma() > 0.0);
            BOOST_CHECK(p.kappa() > 0.0);
            BOOST_CHECK(p.alpha() >= 0.0);
            
            BOOST_CHECK(p.sigma() >= sigmaBounds[0]);
            BOOST_CHECK(p.sigma() <= sigmaBounds[1]);
            BOOST_CHECK(p.kappa() >= kappaBounds[0]);
            BOOST_CHECK(p.kappa() <= kappaBounds[1]);
            BOOST_CHECK(p.alpha() >= alphaBounds[0]);
            BOOST_CHECK(p.alpha() <= alphaBounds[1]);
            
            // these are sanity tests really, aren't they?  well they helped
            // me find a bug.
        }
    }
}
