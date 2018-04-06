/*
 *  testLine.cpp
 *  SimpleMesh
 *
 *  Created by Paul C Hansen on 6/26/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 *
 */


#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestLine

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "SimpleMesh.h"

using namespace std;
using namespace SimpleMesh;

BOOST_AUTO_TEST_CASE( Construction )
{
    SimpleMesh::Line l1(Vector3d::unit(0), Vector3d::unit(1));
    BOOST_CHECK_EQUAL(l1.numControlVertices(), 2);
    
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    Triangle3d xy(o, x, y);
    Triangle3d yz(o, y, z);
    SimpleMesh::Line l2(xy, yz);
    BOOST_CHECK_EQUAL(l2.numControlVertices(), 6);
    
    SimpleMesh::Line l3(xy, yz.plane());
    BOOST_CHECK_EQUAL(l3.numControlVertices(), 3);
}

BOOST_AUTO_TEST_CASE( LineSensitivity )
{
    SimpleMesh::Line l(Vector3d(0,0,0), Vector3d(3,4,5));
    
    BOOST_CHECK(l.x() == Vector3d(0,0,0));
    BOOST_CHECK(l.v() == Vector3d(3,4,5) - Vector3d(0,0,0));
    
    BOOST_CHECK(l.Dx(0) == Matrix3d::eye());
    BOOST_CHECK(l.Dx(1) == Matrix3d::zero());
    
    BOOST_CHECK(l.Dv(0) == -Matrix3d::eye());
    BOOST_CHECK(l.Dv(1) == Matrix3d::eye());
}

BOOST_AUTO_TEST_CASE( TestTriangleTriangleIntersection )
{
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    Triangle3d xy(o, x, y);
    Triangle3d yz(o, y, z);
    
    SimpleMesh::Line l(xy, yz);
    
    // The line is represented by a point and a direction.  Strictly speaking,
    // all I can test in terms of intersection is where the line strikes
    // various planes.
    
    Plane3d yEquals0(Vector3d(0,1,0), -0.0);
    Plane3d yEquals1(Vector3d(0,1,0), -1.0);
    Plane3d yEquals2(Vector3d(0,1,0), -2.0);
    
    Vector3d p0 = intersection(yEquals0, l.x(), l.v());
    Vector3d p1 = intersection(yEquals1, l.x(), l.v());
    Vector3d p2 = intersection(yEquals2, l.x(), l.v());
    
    BOOST_CHECK_SMALL(p0[0], 1e-9);
    BOOST_CHECK_SMALL(p0[1], 1e-9);
    BOOST_CHECK_SMALL(p0[2], 1e-9);
    
    BOOST_CHECK_SMALL(p1[0], 1e-9);
    BOOST_CHECK_CLOSE(p1[1], 1.0, 1e-9);
    BOOST_CHECK_SMALL(p1[2], 1e-9);
    
    BOOST_CHECK_SMALL(p2[0], 1e-9);
    BOOST_CHECK_CLOSE(p2[1], 2.0, 1e-9);
    BOOST_CHECK_SMALL(p2[2], 1e-9);
}

// This is a few by-hand tests to check for silly errors; it's not exhaustive.
BOOST_AUTO_TEST_CASE( TestIntersectionSensitivityFunction )
{
    using SimpleMesh::Line;
    
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    Triangle3d xy(o, x, y);
    Triangle3d yz(o, y, z);
    
    // Walk through the (valid) line-plane intersections and check them all.
    // I can't check the intersection of edges coplanar with the other triangle.
    // This does test for all cases, at least.
    
    // tri0 plane intersect tri1 plane
    BOOST_CHECK_SMALL(norm(Line::intersection(xy, yz, 0) - Vector3d(0,0,0)),
        1e-6);
    BOOST_CHECK_SMALL(norm(Line::intersection(xy, yz, 1) - Vector3d(0,1,0)),
        1e-6);
    BOOST_CHECK_SMALL(dot(yz.normal(), xy.edgeDirection(2)), 1e-6); // coplanar
    
    // tri1 plane intersect tri0 plane
    BOOST_CHECK_SMALL(dot(xy.normal(), yz.edgeDirection(0)), 1e-6); // coplanar
    BOOST_CHECK_SMALL(norm(Line::intersection(yz, xy, 1) - Vector3d(0,1,0)),
        1e-6);
    BOOST_CHECK_SMALL(norm(Line::intersection(yz, xy, 2) - Vector3d(0,0,0)),
        1e-6);
    
    // Now test a few sensitivities.
    
    Triangle3d zedTri(Vector3d(0,0,0), Vector3d(0,0,0), Vector3d(0,0,0));
    
    // Sensitivity of tri0 plane intersect tri1 plane to perturbations in tri0;
    // only looking at edge 0 of tri0.  tri0 contributes the line and tri1
    // contributes the plane, so only perturbations that change the direction
    // of the line will change the intersection point.
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, Test::perturbTri(0,0), zedTri, 0) - Vector3d(0,0,0)),
        1e-6);
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, Test::perturbTri(0,1), zedTri, 0) - Vector3d(0,1,0)),
        1e-6);
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, Test::perturbTri(0,2), zedTri, 0) - Vector3d(0,0,1)),
        1e-6);
    
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, Test::perturbTri(1,0), zedTri, 0) - Vector3d(0,0,0)),
        1e-6);
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, Test::perturbTri(1,1), zedTri, 0) - Vector3d(0,0,0)),
        1e-6);
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, Test::perturbTri(1,2), zedTri, 0) - Vector3d(0,0,0)),
        1e-6);
    
    // Same idea, perturbing tri1 and looking at edge 1.  tri0 contributes the
    // line, and tri1 contributes the plane, so perturbations that do not
    // change the plane will not change the intersection point.
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, zedTri, Test::perturbTri(0,0), 1) - Vector3d(0,0,0)),
        1e-6);
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, zedTri, Test::perturbTri(0,1), 1) - Vector3d(0,0,0)),
        1e-6);
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, zedTri, Test::perturbTri(0,2), 1) - Vector3d(0,0,0)),
        1e-6);
    
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, zedTri, Test::perturbTri(1,0), 1) - Vector3d(1,-1,0)),
        1e-6);
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, zedTri, Test::perturbTri(1,1), 1) - Vector3d(0,0,0)),
        1e-6);
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, zedTri, Test::perturbTri(1,2), 1) - Vector3d(0,0,0)),
        1e-6);
    
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, zedTri, Test::perturbTri(2,0), 1) - Vector3d(0,0,0)),
        1e-6);
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, zedTri, Test::perturbTri(2,1), 1) - Vector3d(0,0,0)),
        1e-6);
    BOOST_CHECK_SMALL(norm(Line::
        dIntersection(xy, yz, zedTri, Test::perturbTri(2,2), 1) - Vector3d(0,0,0)),
        1e-6);
    
}

// Automatic brute-force testing of Dx() and Dv().  Failure here would be harder
// to track down by hand than failure in the simpler test above.
BOOST_AUTO_TEST_CASE( TestDxDv )
{
    // Test the sensitivity of Line::x() and Line::v().
    
    // roughly an xy triangle
    Triangle3d xy(Vector3d(0.2, 0.1, -0.1),
        Vector3d(1.2, -0.1, 0.0),
        Vector3d(-0.01, 1.0, 0.1));
    
    // roughly a yz triangle
    Triangle3d yz(Vector3d(-0.1, 0.0, 0.0),
        Vector3d(0.0, 1.1, -0.03),
        Vector3d(0.0, 0.0, 1.0));
    
    SimpleMesh::Line l(xy, yz);
    const double EPSILON = 1e-8;
    const double BASICALLY_ZERO = EPSILON*1e-6;
    const double MAX_RELATIVE_ERROR = 1e-6;
    const double MAX_ABSOLUTE_ERROR = 1e-6;
    
    // PERTURB XY TRIANGLE
    
    for (int vertex = 0; vertex < 3; vertex++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Triangle3d xyPerturbed = xy;
        Vector3d dVert = EPSILON*Vector3d::unit(xyz);
        xyPerturbed[vertex] += dVert;
        
        SimpleMesh::Line ll(xyPerturbed, yz);
        
        Vector3d deltaX = ll.x() - l.x();
        Vector3d deltaV = ll.v() - l.v();
        
        Vector3d deltaXExpected = l.Dx(vertex)*dVert;
        Vector3d deltaVExpected = l.Dv(vertex)*dVert;
        
        double errX = norm(deltaXExpected - deltaX);
        double errV = norm(deltaVExpected - deltaV);
        
        double relErrX = errX / norm(deltaXExpected);
        double relErrV = errV / norm(deltaVExpected);
        
        if (norm(deltaXExpected) > BASICALLY_ZERO)
            BOOST_CHECK_SMALL(relErrX, MAX_RELATIVE_ERROR);
        else
            BOOST_CHECK_SMALL(errX, MAX_ABSOLUTE_ERROR);
        
        if (norm(deltaVExpected) > BASICALLY_ZERO)
            BOOST_CHECK_SMALL(relErrV, MAX_RELATIVE_ERROR);
        else
            BOOST_CHECK_SMALL(errV, MAX_ABSOLUTE_ERROR);
    }
    
    // PERTURB YZ TRIANGLE
    
    for (int vertex = 0; vertex < 3; vertex++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Triangle3d yzPerturbed = yz;
        Vector3d dVert = EPSILON*Vector3d::unit(xyz);
        yzPerturbed[vertex] += dVert;
        
        SimpleMesh::Line ll(xy, yzPerturbed);
        
        Vector3d deltaX = ll.x() - l.x();
        Vector3d deltaV = ll.v() - l.v();
        
        Vector3d deltaXExpected = l.Dx(vertex+3)*dVert;
        Vector3d deltaVExpected = l.Dv(vertex+3)*dVert;
        
        double errX = norm(deltaXExpected - deltaX);
        double errV = norm(deltaVExpected - deltaV);
        
        double relErrX = errX / norm(deltaXExpected);
        double relErrV = errV / norm(deltaVExpected);
        
        if (norm(deltaXExpected) > BASICALLY_ZERO)
            BOOST_CHECK_SMALL(relErrX, MAX_RELATIVE_ERROR);
        else
            BOOST_CHECK_SMALL(errX, MAX_ABSOLUTE_ERROR);
        
        if (norm(deltaVExpected) > BASICALLY_ZERO)
            BOOST_CHECK_SMALL(relErrV, MAX_RELATIVE_ERROR);
        else
            BOOST_CHECK_SMALL(errV, MAX_ABSOLUTE_ERROR);
    }
}

BOOST_AUTO_TEST_CASE( TestTriPlaneIntersectionSimple )
{
    using SimpleMesh::Line;
    
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    Triangle3d xy(o, x, y);
    Triangle3d yz(o, y, z);
    
    Line ly(xy, yz.plane());
    
    // Check that (0,0,0) and (0,1,0) are on that line
    BOOST_CHECK_SMALL(pointLineDistance(ly.x(), ly.x() + ly.v(), o), 1e-9);
    BOOST_CHECK_SMALL(pointLineDistance(ly.x(), ly.x() + ly.v(), y), 1e-9);
    
    // Check that (1 0 0) and (0 0 1) are at the right distance
    BOOST_CHECK_CLOSE(pointLineDistance(ly.x(), ly.x() + ly.v(), x), 1.0, 1e-9);
    BOOST_CHECK_CLOSE(pointLineDistance(ly.x(), ly.x() + ly.v(), z), 1.0, 1e-9);
}


// Automatic brute-force testing of Dx() and Dv().  Failure here would be harder
// to track down by hand than failure in the simpler test above.
BOOST_AUTO_TEST_CASE( TestDxDv_TriPlane )
{
    // Test the sensitivity of Line::x() and Line::v().
    
    // roughly an xy triangle
    Triangle3d xy(Vector3d(0.2, 0.1, -0.1),
        Vector3d(1.2, -0.1, 0.0),
        Vector3d(-0.01, 1.0, 0.1));
    
    // roughly a yz triangle
    Triangle3d yz(Vector3d(-0.1, 0.0, 0.0),
        Vector3d(0.0, 1.1, -0.03),
        Vector3d(0.0, 0.0, 1.0));
    
//    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
//    Triangle3d xy(o, x, y);
//    Triangle3d yz(o, y, z);
    
    SimpleMesh::Line l(xy, yz.plane());
    const double EPSILON = 1e-8;
    const double BASICALLY_ZERO = EPSILON;
    const double MAX_RELATIVE_ERROR = 1e-6;
    const double MAX_ABSOLUTE_ERROR = 1e-6;
    
    // PERTURB XY TRIANGLE
    
    for (int vertex = 0; vertex < 3; vertex++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Triangle3d xyPerturbed = xy;
        Vector3d dVert = EPSILON*Vector3d::unit(xyz);
        xyPerturbed[vertex] += dVert;
        
//        cerr << "Perturbation of v" << vertex << "_" << xyz << "\n";
        //cerr << xy << " becomes " << xyPerturbed << "\n";
        
        SimpleMesh::Line ll(xyPerturbed, yz.plane());
        
//        cerr << " x from " << l.x() << " to " << ll.x() << "\n";
        
        Vector3d deltaX = ll.x() - l.x();
        Vector3d deltaV = ll.v() - l.v();
        
        Vector3d deltaXExpected = l.Dx(vertex)*dVert;
        Vector3d deltaVExpected = l.Dv(vertex)*dVert;
        
//        cerr << " dx expected " << deltaXExpected << " got " << deltaX << "\n";
        
        double errX = norm(deltaXExpected - deltaX);
        double errV = norm(deltaVExpected - deltaV);
        
        double relErrX = errX / norm(deltaXExpected);
        double relErrV = errV / norm(deltaVExpected);
        
        if (norm(deltaXExpected) > BASICALLY_ZERO)
            BOOST_CHECK_SMALL(relErrX, MAX_RELATIVE_ERROR);
        else
            BOOST_CHECK_SMALL(errX, MAX_ABSOLUTE_ERROR);
        
        if (norm(deltaVExpected) > BASICALLY_ZERO)
            BOOST_CHECK_SMALL(relErrV, MAX_RELATIVE_ERROR);
        else
            BOOST_CHECK_SMALL(errV, MAX_ABSOLUTE_ERROR);
    }
}


//
// Line( Triangle3d, Plane3d )
//
// Test the sensitivity of the Line of intersection of a triangle and a plane.
// (Supporting plane of triangle with plane, rather.)
BOOST_AUTO_TEST_CASE( TrianglePlaneIntersectionSensitivity )
{
    cerr << "==== TrianglePlaneIntersectionSensitivity\n\n";
    
    Vector3d x(1,0.01,-0.01), y(0.015,1,-0.005), z(0.001,-0.001,1);
    Triangle3d tri0(x,y,z);
    Plane3d zEqualsHalf(z, -0.5);
    
    Line line0(tri0, zEqualsHalf);
    
    //cerr << "line0:\n" << line0 << "\n";
    
    for (int vv = 0; vv < 3; vv++)
    for (int dxyz = 0; dxyz < 3; dxyz++)
    {
        const double DELTA = 1e-6;
        
        //cerr << "Perturbation of v" << vv << "_" << dxyz << "\n";
        
        Triangle3d tri1(tri0);
        tri1[vv][dxyz] += DELTA;
        
        Line line1(tri1, zEqualsHalf);
        //cerr << "line1:\n" << line1 << "\n";
        
        Vector3d deltax = line1.x() - line0.x();
        Vector3d deltav = line1.v() - line0.v();
        
        Vector3d deltax_expected = line0.Dx(vv) * Vector3d::unit(dxyz) * DELTA;
        Vector3d deltav_expected = line0.Dv(vv) * Vector3d::unit(dxyz) * DELTA;
        
        for (int nn = 0; nn < 3; nn++)
        if (deltax_expected[nn] == 0)
            BOOST_CHECK_SMALL(deltax[nn], 1e-6);
        else
            BOOST_CHECK_CLOSE(deltax[nn], deltax_expected[nn], 1e-2);
        
        for (int nn = 0; nn < 3; nn++)
        if (deltav_expected[nn] == 0)
            BOOST_CHECK_SMALL(deltav[nn], 1e-6);
        else
            BOOST_CHECK_CLOSE(deltav[nn], deltav_expected[nn], 1e-2);
    }
}




