/*
 *  testClipTriangle.cpp
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 7/2/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */



#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testClipTriangle

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Voxelize/ClipTriangle.h"
#include "SimpleMesh/Search.h"
#include "TimeWrapper.h"

using namespace std;
using namespace TimeWrapper;

//#define CHECK_CLOSE(p0, p1, tol) BOOST_CHECK_SMALL(norm((p0)-(p1)), tol)
#define CHECK_CLOSE(p0, p1) BOOST_CHECK_SMALL(norm((p0)-(p1)), 1e-9)

const void* NULLPTR = 0L;
const int X = 0, Y = 1, Z = 2;

BOOST_AUTO_TEST_CASE( Constructor )
{
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    SimpleMesh::Triangle tri(SimpleMesh::Test::oneTri(o, x, y));
    ClippedTriangle ctri(tri);
    
    BOOST_CHECK_EQUAL(ctri.vertices().size(), 3);
    BOOST_CHECK_EQUAL(ctri.edges().size(), 3);
    BOOST_CHECK_EQUAL(ctri.triangle(), &tri);
    
    for (int edge = 0; edge < 3; edge++)
        BOOST_CHECK_EQUAL(ctri.edges()[edge], &tri.edgeSensitivity(edge));
    
    Rect3d bounds = ctri.bounds();
    BOOST_CHECK(bounds.p1 == o);
    BOOST_CHECK(bounds.p2 == x + y);
}


BOOST_AUTO_TEST_CASE( ClipInPlane )
{
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    SimpleMesh::Triangle tri(SimpleMesh::Test::oneTri(o, x, y));
    ClippedTriangle ctri(tri);
    ClippedTriangle clipped;
    
    const void* NULLPTR = 0L;
    
    // Clip below x = 0.5
    clipped = ctri.partBelow(X, 0.5);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 4);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 4);
    CHECK_CLOSE(clipped.vertices()[0], Vector3d(0,0,0));
    CHECK_CLOSE(clipped.vertices()[1], Vector3d(0.5, 0, 0));
    CHECK_CLOSE(clipped.vertices()[2], Vector3d(0.5, 0.5, 0));
    CHECK_CLOSE(clipped.vertices()[3], Vector3d(0, 1, 0));
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], &tri.edgeSensitivity(0));
    BOOST_CHECK_EQUAL(clipped.edges()[1], NULLPTR);
    BOOST_CHECK_EQUAL(clipped.edges()[2], &tri.edgeSensitivity(1));
    BOOST_CHECK_EQUAL(clipped.edges()[3], &tri.edgeSensitivity(2));
    
    // Clip above x = 0.5
    clipped = ctri.partAbove(X, 0.5);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 3);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 3);
    CHECK_CLOSE(clipped.vertices()[0], Vector3d(0.5,0,0));
    CHECK_CLOSE(clipped.vertices()[1], Vector3d(1,0,0));
    CHECK_CLOSE(clipped.vertices()[2], Vector3d(0.5,0.5,0));
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], &tri.edgeSensitivity(0));
    BOOST_CHECK_EQUAL(clipped.edges()[1], &tri.edgeSensitivity(1));
    BOOST_CHECK_EQUAL(clipped.edges()[2], NULLPTR);
    
    // Clip below y = 0.5
    clipped = ctri.partBelow(Y, 0.5);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 4);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 4);
    CHECK_CLOSE(clipped.vertices()[0], Vector3d(0,0,0));
    CHECK_CLOSE(clipped.vertices()[1], Vector3d(1,0,0));
    CHECK_CLOSE(clipped.vertices()[2], Vector3d(0.5,0.5,0));
    CHECK_CLOSE(clipped.vertices()[3], Vector3d(0,0.5,0));
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], &tri.edgeSensitivity(0));
    BOOST_CHECK_EQUAL(clipped.edges()[1], &tri.edgeSensitivity(1));
    BOOST_CHECK_EQUAL(clipped.edges()[2], NULLPTR);
    BOOST_CHECK_EQUAL(clipped.edges()[3], &tri.edgeSensitivity(2));
    
    // Clip above y = 0.5
    clipped = ctri.partAbove(Y, 0.5);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 3);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 3);
    CHECK_CLOSE(clipped.vertices()[0], Vector3d(0.5,0.5,0));
    CHECK_CLOSE(clipped.vertices()[1], Vector3d(0,1,0));
    CHECK_CLOSE(clipped.vertices()[2], Vector3d(0,0.5,0));
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], &tri.edgeSensitivity(1));
    BOOST_CHECK_EQUAL(clipped.edges()[1], &tri.edgeSensitivity(2));
    BOOST_CHECK_EQUAL(clipped.edges()[2], NULLPTR);
}

// Vertices on the boundary are considered to be outside.  This should prevent
// the existence of degenerate polygons (just boundary segments).
//
// This triangle has an edge at x = 0, an edge at y = 1, and a diagonal.
//
// +------+
// |     /
// |    /
// |   /
// |  /
// | /
// |/
// +
// 
// My clipping convention is that partBelow keeps stuff with, e.g., x < 1.0,
// and partAbove keeps stuff with, e.g., x >= 0.0.  All necessary vertices
// will be preserved, but the edge sensitivity will saved for the x = 0 line
// and not for the y = 1 line.  This way the edge of a triangle along a voxel
// boundary will contribute to fill factor sensitivities only on one side at
// a time.
BOOST_AUTO_TEST_CASE( VerticesOnBoundary )
{
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1), xy(1,1,0);
    SimpleMesh::Triangle tri(SimpleMesh::Test::oneTri(o, xy, y));
    ClippedTriangle ctri(tri);
    ClippedTriangle clipped;
    
    // Clip below x = 0.0.  The vertices on the edge are included but the edge
    // is explicitly not, so the edge pointers are null.
    clipped = ctri.partBelow(X, 0.0);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 2);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 2);
    
    CHECK_CLOSE(clipped.vertices()[0], o);
    CHECK_CLOSE(clipped.vertices()[1], y);
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], NULLPTR);
    BOOST_CHECK_EQUAL(clipped.edges()[1], NULLPTR);
    
    // Clip above x = 0.0.  All points with x >= 0.0 are included (all of them).
    clipped = ctri.partAbove(X,0.0);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 3);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 3);
    
    CHECK_CLOSE(clipped.vertices()[0], o);
    CHECK_CLOSE(clipped.vertices()[1], xy);
    CHECK_CLOSE(clipped.vertices()[2], y);
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], &tri.edgeSensitivity(0));
    BOOST_CHECK_EQUAL(clipped.edges()[1], &tri.edgeSensitivity(1));
    BOOST_CHECK_EQUAL(clipped.edges()[2], &tri.edgeSensitivity(2));
    
    // Clip below x = 1.0.  All points with x < 1.0 are included.
    clipped = ctri.partBelow(X, 1.0);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 3);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 3);
    
    CHECK_CLOSE(clipped.vertices()[0], o);
    CHECK_CLOSE(clipped.vertices()[1], xy);
    CHECK_CLOSE(clipped.vertices()[2], y);
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], &tri.edgeSensitivity(0));
    BOOST_CHECK_EQUAL(clipped.edges()[1], &tri.edgeSensitivity(1));
    BOOST_CHECK_EQUAL(clipped.edges()[2], &tri.edgeSensitivity(2));
    
    // Clip above x = 1.0.  All points with x >= 1.0 are included.  This will
    // initially toss the point at (1,0,0) into the polygon by itself, but the 
    // finaly polygon will be empty because degenerate polygons get BALEETED.
    clipped = ctri.partAbove(X, 1.0);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 0);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 0);
    
    // Clip below y = 0.0.  Should be empty.
    clipped = ctri.partBelow(Y, 0.0);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 0);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 0);
    
    // Clip above y = 0.0.  Should be the whole triangle.
    clipped = ctri.partAbove(Y, 0.0);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 3);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 3);
    
    CHECK_CLOSE(clipped.vertices()[0], o);
    CHECK_CLOSE(clipped.vertices()[1], xy);
    CHECK_CLOSE(clipped.vertices()[2], y);
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], &tri.edgeSensitivity(0));
    BOOST_CHECK_EQUAL(clipped.edges()[1], &tri.edgeSensitivity(1));
    BOOST_CHECK_EQUAL(clipped.edges()[2], &tri.edgeSensitivity(2));
    
    // Clip below y = 1.0.  Should be the whole triangle but exclude the
    // sensitivity of the top edge.
    clipped = ctri.partBelow(Y, 1.0);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 3);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 3);
    
    CHECK_CLOSE(clipped.vertices()[0], o);
    CHECK_CLOSE(clipped.vertices()[1], xy);
    CHECK_CLOSE(clipped.vertices()[2], y);
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], &tri.edgeSensitivity(0));
    BOOST_CHECK_EQUAL(clipped.edges()[1], NULLPTR);
    BOOST_CHECK_EQUAL(clipped.edges()[2], &tri.edgeSensitivity(2));
    
    // Clip above y = 1.0.  Should be an isolated edge (i.e. a forward edge
    // and a backward edge).  The forward edge has sensitivity, and the
    // backward edge does not.
    clipped = ctri.partAbove(Y, 1.0);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 2);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 2);
    
    CHECK_CLOSE(clipped.vertices()[0], xy);
    CHECK_CLOSE(clipped.vertices()[1], y);
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], &tri.edgeSensitivity(1));
    BOOST_CHECK_EQUAL(clipped.edges()[1], NULLPTR);
}

BOOST_AUTO_TEST_CASE( TriangleParallelToClippingPlane )
{
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1), xy(1,1,0);
    SimpleMesh::Triangle tri(SimpleMesh::Test::oneTri(o, xy, y));
    ClippedTriangle ctri(tri);
    ClippedTriangle clipped;
    
    // Clip below z = 0.0.  The vertices are preserved and the edges do NOT
    // save sensitivity information.  I'm not sure if this behavior is good or
    // not—mostly it's a side-effect of the clipping behavior for vertices
    // on clipping planes (it's not by design).  But is it harmful as such?
    // Probably not?  Anyway, this IS what it does:
    clipped = ctri.partBelow(Z, 0.0);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 3);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 3);
    
    CHECK_CLOSE(clipped.vertices()[0], o);
    CHECK_CLOSE(clipped.vertices()[1], xy);
    CHECK_CLOSE(clipped.vertices()[2], y);
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], NULLPTR);
    BOOST_CHECK_EQUAL(clipped.edges()[1], NULLPTR);
    BOOST_CHECK_EQUAL(clipped.edges()[2], NULLPTR);
    
    // Clip above z = 0.0.  This should be the whole triangle.
    clipped = ctri.partAbove(Z,0.0);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 3);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 3);
    
    CHECK_CLOSE(clipped.vertices()[0], o);
    CHECK_CLOSE(clipped.vertices()[1], xy);
    CHECK_CLOSE(clipped.vertices()[2], y);
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], &tri.edgeSensitivity(0));
    BOOST_CHECK_EQUAL(clipped.edges()[1], &tri.edgeSensitivity(1));
    BOOST_CHECK_EQUAL(clipped.edges()[2], &tri.edgeSensitivity(2));
}

BOOST_AUTO_TEST_CASE( ClipTwice )
{
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1), xy(1,1,0);
    SimpleMesh::Triangle triTwo(SimpleMesh::Test::oneTri(o, 2*xy, 2*y));
    ClippedTriangle ctri(triTwo);
    ClippedTriangle clipped;
    
    // CLIP X < 1.0
    clipped = ctri.partBelow(X, 1.0);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 4);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 4);
    CHECK_CLOSE(clipped.vertices()[0], o);
    CHECK_CLOSE(clipped.vertices()[1], xy);
    CHECK_CLOSE(clipped.vertices()[2], xy + y);
    CHECK_CLOSE(clipped.vertices()[3], 2*y);
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], &triTwo.edgeSensitivity(0));
    BOOST_CHECK_EQUAL(clipped.edges()[1], NULLPTR);
    BOOST_CHECK_EQUAL(clipped.edges()[2], &triTwo.edgeSensitivity(1));
    BOOST_CHECK_EQUAL(clipped.edges()[3], &triTwo.edgeSensitivity(2));
    
    // CLIP Y < 1.0
    clipped = clipped.partBelow(Y, 1.0);
    BOOST_CHECK_EQUAL(clipped.vertices().size(), 3);
    BOOST_CHECK_EQUAL(clipped.edges().size(), 3);
    CHECK_CLOSE(clipped.vertices()[0], o);
    CHECK_CLOSE(clipped.vertices()[1], xy);
    CHECK_CLOSE(clipped.vertices()[2], y);
    
    BOOST_CHECK_EQUAL(clipped.edges()[0], &triTwo.edgeSensitivity(0));
    BOOST_CHECK_EQUAL(clipped.edges()[1], NULLPTR);
    BOOST_CHECK_EQUAL(clipped.edges()[2], &triTwo.edgeSensitivity(2));
}

BOOST_AUTO_TEST_CASE( Speed )
{
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1), xy(1,1,0);
    SimpleMesh::Triangle tri(SimpleMesh::Test::oneTri(o, xy, y));
    ClippedTriangle ctri(tri);
    ClippedTriangle clipped;
    
    int numRepetitions = 100;
    int numPositions = 101;
    double delta = 1.0/(numPositions-1);
    
    TimePoint t0 = now();
    for (int repetition = 0; repetition < numRepetitions; repetition++)
    for (int axis = X; axis <= Z; axis++)
    for (double planePos = 0.0; planePos <= 1.0; planePos += delta)
    {
        clipped = ctri.partAbove(axis, planePos);
        clipped = ctri.partBelow(axis, planePos);
    }
    double sec = elapsedSeconds(t0, now());
    
    cout << "Total clips: " << 3*numRepetitions*numPositions*2
        << " in " << sec << " seconds.\n";
    cout << 3*numRepetitions*numPositions*2/sec << " calls/sec.\n";
    cout << sec/(3*numRepetitions*numPositions*2) << " sec/call.\n";
}

BOOST_AUTO_TEST_CASE( AreaAndVolume )
{
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1), xy(1,1,0);
    SimpleMesh::Triangle tri(SimpleMesh::Test::oneTri(o + z/2, 2*x + z/2, 2*y + z/2));
    ClippedTriangle ctri(tri);
    ctri = ctri.partBelow(X, 1.0);
    ctri = ctri.partBelow(Y, 1.0);
    
    // ctri should now be the square [0, 1]x[0, 1], elevated 0.5 along z
    
    BOOST_CHECK_CLOSE(ctri.area2d(Z), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(ctri.volume(Z), 0.5, 1e-10);
}










