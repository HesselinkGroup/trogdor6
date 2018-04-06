/*
 *  testAssembly.cpp
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 6/21/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 *
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestTriangleFacet

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "SimpleMesh.h"

using namespace std;
using namespace SimpleMesh;


BOOST_AUTO_TEST_CASE( TriangleConstruction )
{
    BOOST_CHECK(true);
    
    Triangle emptyTri;
    BOOST_CHECK(emptyTri.triangle() == Triangle3d());
    BOOST_CHECK(emptyTri.neighbors() == Vector3i(-1,-1,-1));
    
    for (int ee = 0; ee < 3; ee++)
        BOOST_CHECK_EQUAL(emptyTri.controlVertices()[ee].id(), -1);
    
    for (int ee = 0; ee < 3; ee++)
    {
        BOOST_CHECK_EQUAL(emptyTri.edgeControlVertices(ee).size(), 0);
    }
    
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0);
    Triangle3d tri(o, x, y);
    
    Triangle meshTriangle(tri);
    BOOST_CHECK(meshTriangle.triangle() == tri);
    BOOST_CHECK(meshTriangle.unitNormal() == unit(tri.normal()));
    BOOST_CHECK(meshTriangle.neighbors() == Vector3i(-1,-1,-1));
    
    for (int ee = 0; ee < 3; ee++)
        BOOST_CHECK_EQUAL(meshTriangle.controlVertices()[ee].id(), -1);
    
    for (int ee = 0; ee < 3; ee++)
    {
        BOOST_CHECK_EQUAL(meshTriangle.edgeControlVertices(ee).size(), 0);
    }
    
    Triangle meshTriangle2;
    meshTriangle2.triangle(tri);
    BOOST_CHECK(meshTriangle2.triangle() == tri);
    BOOST_CHECK(meshTriangle2.unitNormal() == unit(tri.normal()));
    for (int ee = 0; ee < 3; ee++)
    {
        BOOST_CHECK_EQUAL(meshTriangle2.edgeControlVertices(ee).size(), 0);
    }
}

// TrackedPolyhedron will typically work like this:
//
// 1. Triangulate
// 2. Find neighbor triangles
// 3. Inherit control vertices
// 4. Cache sensitivities
//
// I should follow these steps here.
//
// Assemble a simplex (convex hull of 0, unit(x), unit(y), and unit(z)).
// Build the whole mesh, including control vertices.
// Test the sensitivities.  Yeah.
//
BOOST_AUTO_TEST_CASE( CacheSensitivities )
{
    cerr << "Warning: CacheSensitivities test does not really do anything yet.\n";
    
    Vector3b threeTrues(true, true, true);
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    Triangle3d tri(o, x, y);
    
    vector<ControlVertex> cv;
    cv.push_back(ControlVertex(o, threeTrues, 0));
    cv.push_back(ControlVertex(x, threeTrues, 1));
    cv.push_back(ControlVertex(y, threeTrues, 2));
    cv.push_back(ControlVertex(z, threeTrues, 3));
    
    vector<Triangle> triangles;
    triangles.push_back(Triangle(Triangle3d(o, y, x))); // xy plane
    triangles.push_back(Triangle(Triangle3d(o, z, y))); // yz plane
    triangles.push_back(Triangle(Triangle3d(o, x, z))); // zx plane
    triangles.push_back(Triangle(Triangle3d(x, y, z))); // the diagonal surface
    
    triangles[0].controlVertices(cv[0], cv[2], cv[1]);
    triangles[0].edgeControlVertices(cv[0], cv[2], cv[1]);
    
    triangles[1].controlVertices(cv[0], cv[3], cv[2]);
    triangles[1].edgeControlVertices(cv[0], cv[3], cv[2]);
    
    triangles[2].controlVertices(cv[0], cv[1], cv[3]);
    triangles[2].edgeControlVertices(cv[0], cv[1], cv[3]);
    
    triangles[3].controlVertices(cv[1], cv[2], cv[3]);
    triangles[3].edgeControlVertices(cv[1], cv[2], cv[3]);
    
//    triangles[0].controlVertices(Vector3i(0, 2, 1));
//    triangles[0].edgeControlVertices(Vector3i(0, 2, 1));
//    
//    triangles[1].controlVertices(Vector3i(0, 3, 2));
//    triangles[1].edgeControlVertices(Vector3i(0, 3, 2));
//    
//    triangles[2].controlVertices(Vector3i(0, 1, 3));
//    triangles[2].edgeControlVertices(Vector3i(0, 1, 3));
//    
//    triangles[3].controlVertices(Vector3i(1, 2, 3));
//    triangles[3].edgeControlVertices(Vector3i(1, 2, 3));
    
    for (int tt = 0; tt < triangles.size(); tt++)
    {
        triangles[tt].cacheSensitivity();
        
        for (int ee = 0; ee < 3; ee++)
            BOOST_CHECK_EQUAL(triangles[tt].edgeControlVertices(ee).size(), 2);
    }
}

BOOST_AUTO_TEST_CASE( OrientationSensitivity_Positive )
{
//    cerr << "\n==== OrientationSensitivity_Positive\n\n";
    Vector3b threeTrues(true, true, true);
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    Triangle3d tri(o + 0.1*x, x+0.5*z, y-0.3*z); // Use simple triangle to debug
    
    Triangle meshTri0(Test::oneTri(tri[0], tri[1], tri[2]));
    
//    cerr << "Ori = \n" << meshTri0.nnTdA() << "\n";
    
    Matrix3<Vector3d> DNNT[3] = { meshTri0.Dorientation(0),
        meshTri0.Dorientation(1),
        meshTri0.Dorientation(2) };
    
    Matrix3<Vector3d> deltaNNT[3];
    
//    cerr << "DOri0 = \n" << DNNT[0] << "\n";
    
    double epsilon = 1e-8;
    for (int vv = 0; vv < 3; vv++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Triangle3d tri1(tri);
        tri1[vv][xyz] += epsilon;
        
        Triangle meshTri1(Test::oneTri(tri1[0], tri1[1], tri1[2]));
        
        Matrix3d deltaOr = (meshTri1.nnTdA() - meshTri0.nnTdA()) / epsilon;
        
        for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
        {
            deltaNNT[vv](ii,jj)[xyz] = deltaOr(ii,jj);
            
            if (DNNT[vv](ii,jj)[xyz] == 0.0)
                BOOST_CHECK_SMALL(deltaOr(ii,jj), 1e-6);
            else
                BOOST_CHECK_CLOSE(DNNT[vv](ii,jj)[xyz], deltaOr(ii,jj), 1e-4);
        }
    }
    
//    cerr << "DeltaOri0 = \n" << deltaNNT[0] << "\n";
}

BOOST_AUTO_TEST_CASE( OrientationSensitivity_Negative )
{
//    cerr << "\n==== OrientationSensitivity_Negative\n\n";
    Vector3b threeTrues(true, true, true);
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    Triangle3d tri(o + 0.1*x, x+0.5*z, y-0.3*z); // Use simple triangle to debug
    
    Triangle meshTri0(Test::oneFlippedTri(tri[0], tri[1], tri[2]));
    
//    cerr << "Ori = \n" << meshTri0.nnTdA() << "\n";
    
    Matrix3<Vector3d> DNNT[3] = { meshTri0.Dorientation(0),
        meshTri0.Dorientation(1),
        meshTri0.Dorientation(2) };
    
    Matrix3<Vector3d> deltaNNT[3];
    
//    cerr << "DOri0 = \n" << DNNT[0] << "\n";
    
    double epsilon = 1e-8;
    for (int vv = 0; vv < 3; vv++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Triangle3d tri1(tri);
        tri1[vv][xyz] += epsilon;
        
        Triangle meshTri1(Test::oneFlippedTri(tri1[0], tri1[1], tri1[2]));
        
        Matrix3d deltaOr = (meshTri1.nnTdA() - meshTri0.nnTdA()) / epsilon;
        
        for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
        {
            deltaNNT[vv](ii,jj)[xyz] = deltaOr(ii,jj);
            
            if (DNNT[vv](ii,jj)[xyz] == 0.0)
                BOOST_CHECK_SMALL(deltaOr(ii,jj), 1e-6);
            else
                BOOST_CHECK_CLOSE(DNNT[vv](ii,jj)[xyz], deltaOr(ii,jj), 1e-4);
        }
    }
    
//    cerr << "DeltaOri0 = \n" << deltaNNT[0] << "\n";
}



BOOST_AUTO_TEST_CASE( PrintSizes )
{
    BOOST_CHECK(true);
    return;
    cerr << "Line " << sizeof(Line) << " bytes\n";
    cerr << "SensitiveEdge " << sizeof(SensitiveEdge) << " bytes\n";
    cerr << "Triangle " << sizeof(Triangle) << " bytes\n";
    cerr << "ControlVertex " << sizeof(ControlVertex) << " bytes\n";
    cerr << "10k triangle mesh: roughly "
        << 1e-6*(10000*sizeof(Triangle) + 30000*sizeof(ControlVertex))
        << " MB\n";
    
    BOOST_CHECK(true);
}










