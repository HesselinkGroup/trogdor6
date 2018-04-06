/*
 *  testPMLParameters.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 4/1/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test3D

#include <boost/test/unit_test.hpp>
#include <cmath>

#include "PMLParameters.h"
#include "PhysicalConstants.h"

using namespace std;

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


