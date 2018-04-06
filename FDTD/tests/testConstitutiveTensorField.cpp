/*
 *  testConstitutiveTensorField.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 1/11/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestConstitutiveTensorField

#include <boost/test/unit_test.hpp>

#include "ConstitutiveTensorField.h"

using namespace std;
using namespace RLE;
using namespace YeeUtilities;

BOOST_AUTO_TEST_CASE( ConstructorEtCetera )
{
    Rect3i yeeCells(0,0,0,10,10,0);
    
    ConstitutiveTensorField c1;
//    BOOST_CHECK(c1.fieldOctant() == Octant::electric());
//    BOOST_CHECK(c1.fieldOctant(2,1) == Octant::electric()(2,1));
    for (int xyz = 0; xyz < 3; xyz++)
    {
//        BOOST_CHECK(c1.halfCellOffset(xyz,xyz) == eFieldOffset(xyz));
//        for (int abc = 0; abc < 3; abc++)
//        if (abc != xyz).
//            BOOST_CHECK(c1.halfCellOffset(xyz,abc) == Vector3i(0,0,0));
    }
    
    ConstitutiveTensorField c2(Octant::magnetic(), yeeCells);
//    BOOST_CHECK(c2.fieldOctant() == Octant::magnetic());
//    BOOST_CHECK(c2.fieldOctant(1,0) == Octant::magnetic()(1,0));
    for (int xyz = 0; xyz < 3; xyz++)
    {
//        BOOST_CHECK(c2.halfCellOffset(xyz,xyz) == hFieldOffset(xyz));
//        for (int abc = 0; abc < 3; abc++)
//        if (abc != xyz)
//            BOOST_CHECK(c2.halfCellOffset(xyz,abc) == Vector3i(1,1,1));
    }
    
    for (int xyz = 0; xyz < 3; xyz++)
    for (int abc = 0; abc < 3; abc++)
    {
        BOOST_CHECK_EQUAL(c2(xyz,abc).hasDimension(0), true);
        BOOST_CHECK_EQUAL(c2(xyz,abc).hasDimension(1), true);
        BOOST_CHECK_EQUAL(c2(xyz,abc).hasDimension(2), false);
    }
}

static Precision::RationalFunction ofOrder(int numerOrder, int denomOrder)
{
    Precision::RationalFunction r;
    for (int nn = 0; nn <= numerOrder; nn++)
        r.numerator().coefficient(nn, nn+1);
    for (int dd = 0; dd <= denomOrder; dd++)
        r.denominator().coefficient(dd, dd+1);
    return r;
}

BOOST_AUTO_TEST_CASE( SupportAndFilter )
{
    // Test each octant separately
    for (int xyz = 0; xyz < 3; xyz++)
    for (int abc = 0; abc < 3; abc++)
    {
        ConstitutiveTensorField c;
        c(xyz,abc).mark(0,0,0, ofOrder(1,0));
        c(xyz,abc).mark(0,1,0, ofOrder(1,1));
        c(xyz,abc).mark(0,2,0, ofOrder(0,1));
        
        DynamicRLE3<Precision::RationalFunction> order10, order11, order01;
        SupportRegion3 support10, support11, support01;
        
        order10.mark(0,0,0, ofOrder(1,0));
        order11.mark(0,1,0, ofOrder(1,1));
        order01.mark(0,2,0, ofOrder(0,1));
        support10 = SupportRegion3(order10);
        support11 = SupportRegion3(order11);
        support01 = SupportRegion3(order01);
        
        BOOST_CHECK(c.support(xyz,abc,1,0) == support10);
        BOOST_CHECK(c.support(xyz,abc,0,1) == support01);
        BOOST_CHECK(c.support(xyz,abc,1,1) == support11);
        
        for (int yy = 0; yy < 3; yy++)
        {
            BOOST_CHECK(c.filter(xyz,abc,1,0).at(0,yy,0) == order10.at(0,yy,0));
            BOOST_CHECK(c.filter(xyz,abc,0,1).at(0,yy,0) == order01.at(0,yy,0));
            BOOST_CHECK(c.filter(xyz,abc,1,1).at(0,yy,0) == order11.at(0,yy,0));
        }
        
        // Wildcards too.
        // Find all the cells with numerator order 1, and ANY denominator.
        // Find all the cells with denominator order 1, and ANY numerator.
        DynamicRLE3<Precision::RationalFunction> order1_, order_1;
        SupportRegion3 support1_, support_1;
        
        order1_.mark(0,0,0, ofOrder(1,0));
        order1_.mark(0,1,0, ofOrder(1,1));
        order_1.mark(0,1,0, ofOrder(1,1));
        order_1.mark(0,2,0, ofOrder(0,1));
        support1_ = SupportRegion3(order1_);
        support_1 = SupportRegion3(order_1);
        
        BOOST_CHECK(c.support(xyz,abc,1,-1) == support1_);
        BOOST_CHECK(c.support(xyz,abc,-1,1) == support_1);
        
        for (int yy = 0; yy < 3; yy++)
        {
            BOOST_CHECK(c.filter(xyz,abc,1,-1).at(0,yy,0) == order1_.at(0,yy,0));
            BOOST_CHECK(c.filter(xyz,abc,-1,1).at(0,yy,0) == order_1.at(0,yy,0));
        }
    }
}

BOOST_AUTO_TEST_CASE( OffDiagonalSupport_Simple )
{
//    Rect3i physicalYeeCells(-10,-10,0,10,10,0);
//    NodeExtents extents(physicalYeeCells, yeeToHalf(physicalYeeCells),
//        yeeToHalf(physicalYeeCells), yeeToHalf(physicalYeeCells));
//    
//    ConstitutiveTensorField eps(YeeUtilities::Octant::electric(),
//        physicalYeeCells);
//    ConstitutiveTensorField mu(YeeUtilities::Octant::magnetic(),
//        physicalYeeCells);
//    
//    DynamicRLE3<Precision::RationalFunction> gridMu(true,true,false),
//        gridEps(true,true,false);
//    
//    gridMu.mark(-10,-10,0,10,10,0, 1.0);
//    gridEps.mark(-10,-10,0,10,10,0, 1.0);
//    
//    for (int ii = 0; ii < 3; ii++)
//    for (int jj = 0; jj < 3; jj++)
//        eps(ii,jj) = gridEps;
//    
//    for (int xyz = 0; xyz < 3; xyz++)
//        mu(xyz,xyz) = gridMu;
//    
//    SupportRegion3 everywhere(gridMu);
//    
//    BOOST_CHECK(eps.offDiagonalSupport() == everywhere);
//    BOOST_CHECK(mu.offDiagonalSupport().numRuns() == 0);
}

BOOST_AUTO_TEST_CASE( Extrude )
{
    BOOST_CHECK(true);
    // This functionality should actually move to OrientedVoxels, I think.
}

BOOST_AUTO_TEST_CASE( Copy )
{
    BOOST_CHECK(true);
    // Low priority now (Jan 11 2011).
}









