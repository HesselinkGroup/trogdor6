/*
 *  testRowVoxelizer.cpp
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 7/6/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testRowVoxelizer

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Voxelize/RowVoxelizer.h"
#include "TimeWrapper.h"
#include "rle/DynamicRLE3.h"

using namespace std;
using namespace RLE;

#define CHECK_CLOSE(p0, p1) BOOST_CHECK_SMALL(norm((p0)-(p1)), 1e-9)
#define CHECK_MATRIX_CLOSE(m0, m1) BOOST_CHECK_SMALL(normInf((m0)-(m1)), 1e-9);

const void* NULLPTR = 0L;
static const int X = 0, Y = 1, Z = 2;
    
BOOST_AUTO_TEST_CASE( Constructor )
{
    Rect3d bounds(-3,-2,-1,3,2,1);
    Rect3i voxels(-1,-1,-1,1,1,1);
    Vector3b nonSymmetricDimensions(true,true,true);
    RowVoxelizer rv(Y, bounds, voxels, nonSymmetricDimensions);
    
    BOOST_CHECK_EQUAL(rv.rowAxis(), Y);
    BOOST_CHECK(rv.bounds() == bounds);
    BOOST_CHECK(rv.voxelBounds() == voxels);
    BOOST_CHECK(rv.voxelSize() == Vector3d(2.0, 4.0/3.0, 2.0/3.0));
    
    // Voxel extents:
    // X pixels: [-3 -1] [-1 1] [1 3]
    // Y pixels: [-2 -2/3] [-2/3 2/3] [2/3 2]
    // Z pixels: [-1 -1/3] [-1/3 1/3] [1/3 1]
    
    CHECK_CLOSE(rv.pixelBounds(0,0).p1, Vector2d(-1.0/3, -1));
    CHECK_CLOSE(rv.pixelBounds(0,0).p2, Vector2d(1.0/3, 1));
    CHECK_CLOSE(rv.pixelBounds(1,0).p1, Vector2d(1.0/3, -1));
    CHECK_CLOSE(rv.pixelBounds(1,0).p2, Vector2d(1.0, 1));
    CHECK_CLOSE(rv.pixelBounds(0,-1).p1, Vector2d(-1.0/3, -3));
    CHECK_CLOSE(rv.pixelBounds(0,-1).p2, Vector2d(1.0/3, -1));
    
    BOOST_CHECK_EQUAL(rv.voxel(0.5), 0);
    BOOST_CHECK_EQUAL(rv.voxel(1.5), 1);
    BOOST_CHECK_EQUAL(rv.voxel(-0.5), 0);
    
    BOOST_CHECK_CLOSE(rv.voxelMin(1), 2.0/3, 1e-9);
    BOOST_CHECK_CLOSE(rv.voxelMin(-1), -2.0, 1e-9);
}

BOOST_AUTO_TEST_CASE( PixelBounds )
{
    Rect3d bounds(-3,-2,-1,3,2,1);
    Rect3i voxels(-1,-1,-1,1,1,1);
    Vector3b nonSymmetricDimensions(true,true,true);
    RowVoxelizer rv(Z, bounds, voxels, nonSymmetricDimensions);
    
    // Test the lower-front-bottom pixel.
    
    Rect2d pixExpected(-3,-2,-1,-2.0/3.0);
    
    CHECK_CLOSE(rv.pixelBounds(-1,-1).p1, pixExpected.p1);
    CHECK_CLOSE(rv.pixelBounds(-1,-1).p2, pixExpected.p2);
}

BOOST_AUTO_TEST_CASE( ClipToRow )
{
    Rect3d bounds(0,0,0,1,1,1);
    Rect3d voxels(0,0,0, 9, 99, 9);
    Vector3b nonSymmetricDimensions(true,true,true);
    RowVoxelizer rv(Z, bounds, voxels, nonSymmetricDimensions);
    
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    SimpleMesh::Triangle triangleTwo(SimpleMesh::Test::oneTri(o, 2*x, 2*y));
    
    ClippedTriangle clipped = rv.clipToRow(triangleTwo, 0, 0);
    
    CHECK_CLOSE(clipped.bounds().p1, Vector3d(0,0,0));
    CHECK_CLOSE(clipped.bounds().p2, Vector3d(0.1,0.01,0));
}

//BOOST_AUTO_TEST_CASE( FillVoxelizeOneTri )
//{
//    Rect3d bounds(0,0,0,1,1,1);
//    Rect3d voxels(0,0,0,0,0,0);
//    DynamicRLE3<double> fillFactors(true, true, true);
//    FillFactorRowVoxelizer rv(Z, bounds, voxels, fillFactors);
//    
//    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
//    SimpleMesh::Triangle tri(SimpleMesh::Test::oneTri(o + z/2, x + z/2, y + z/2));
//    
//    vector<FillRecord> fillRecords;
//    rv.voxelizeOneTri(&tri, 0, 0, fillRecords);
//    
//    BOOST_CHECK_EQUAL(fillRecords.size(), 1);
//    BOOST_CHECK_EQUAL(fillRecords[0].voxel, 0);
//    BOOST_CHECK_CLOSE(fillRecords[0].area, 0.5, 1e-10);
//    BOOST_CHECK_CLOSE(fillRecords[0].volume, 0.25, 1e-10);
//    fillRecords.clear();
//    
//    SimpleMesh::Triangle tippedTri(SimpleMesh::Test::oneTri(o + 0.000001*z,
//        x + 0.999999*z, y + 0.999999*z) );
//    rv.voxelizeOneTri(&tippedTri, 0, 0, fillRecords);
//    
//    BOOST_CHECK_EQUAL(fillRecords.size(), 1);
//    BOOST_CHECK_EQUAL(fillRecords[0].voxel, 0);
//    BOOST_CHECK_CLOSE(fillRecords[0].area, 0.5, 1e-2);
//    BOOST_CHECK_CLOSE(fillRecords[0].volume, 1.0/3, 1e-2); // 1/3 base*height
//    fillRecords.clear();
//    
//    SimpleMesh::Triangle tallTri(SimpleMesh::Test::oneTri(o - 0.999999*z,
//        x + 0.999999*z, y + 0.99999*z));
//    rv.voxelizeOneTri(&tallTri, 0, 0, fillRecords);
//    
//    BOOST_CHECK_EQUAL(fillRecords.size(), 2);
//    BOOST_CHECK_EQUAL(fillRecords[0].voxel, -1);
//    BOOST_CHECK_CLOSE(fillRecords[0].area, 1.0/8, 1e-2);
//    BOOST_CHECK_CLOSE(fillRecords[0].volume, 2.0/8/3, 1e-2);
//    BOOST_CHECK_EQUAL(fillRecords[1].voxel, 0);
//    BOOST_CHECK_CLOSE(fillRecords[1].area, 0.5 - 0.125, 1e-2);
//    BOOST_CHECK_CLOSE(fillRecords[1].volume, 5.0/24, 1e-2);
//}

BOOST_AUTO_TEST_CASE( FillVoxelizeOneTri )
{
    Rect3d bounds(0,0,0,1,1,1);
    Rect3d voxels(-1,-1,-1,-1,-1,-1);
    Vector3b nonSymmetricDimensions(true,true,true);
    DynamicRLE3<double> fillFactors(nonSymmetricDimensions.asArray());
    FillFactorRowVoxelizer rv(Z, bounds, voxels, nonSymmetricDimensions, fillFactors);
    
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    SimpleMesh::Triangle tri(SimpleMesh::Test::oneTri(o + z/2, x + z/2, y + z/2));
    
    vector<FillRecord> fillRecords;
    rv.voxelizeOneTri(&tri, -1, -1, fillRecords);
    
    BOOST_CHECK_EQUAL(fillRecords.size(), 1);
    BOOST_CHECK_EQUAL(fillRecords[0].voxel, -1);
    BOOST_CHECK_CLOSE(fillRecords[0].area, 0.5, 1e-10);
    BOOST_CHECK_CLOSE(fillRecords[0].volume, 0.25, 1e-10);
    fillRecords.clear();
    
    SimpleMesh::Triangle tippedTri(SimpleMesh::Test::oneTri(o + 0.000001*z,
        x + 0.999999*z, y + 0.999999*z) );
    rv.voxelizeOneTri(&tippedTri, -1, -1, fillRecords);
    
    BOOST_CHECK_EQUAL(fillRecords.size(), 1);
    BOOST_CHECK_EQUAL(fillRecords[0].voxel, -1);
    BOOST_CHECK_CLOSE(fillRecords[0].area, 0.5, 1e-2);
    BOOST_CHECK_CLOSE(fillRecords[0].volume, 1.0/3, 1e-2); // 1/3 base*height
    fillRecords.clear();
    
    SimpleMesh::Triangle tallTri(SimpleMesh::Test::oneTri(o - 0.999999*z,
        x + 0.999999*z, y + 0.99999*z));
    rv.voxelizeOneTri(&tallTri, -1, -1, fillRecords);
    
    BOOST_CHECK_EQUAL(fillRecords.size(), 2);
    BOOST_CHECK_EQUAL(fillRecords[0].voxel, -2);
    BOOST_CHECK_CLOSE(fillRecords[0].area, 1.0/8, 1e-2);
    BOOST_CHECK_CLOSE(fillRecords[0].volume, 2.0/8/3, 1e-2);
    BOOST_CHECK_EQUAL(fillRecords[1].voxel, -1);
    BOOST_CHECK_CLOSE(fillRecords[1].area, 0.5 - 0.125, 1e-2);
    BOOST_CHECK_CLOSE(fillRecords[1].volume, 5.0/24, 1e-2);
}

BOOST_AUTO_TEST_CASE( HandleSortedFillRecords )
{
    Rect3d bounds(0,0,0,4,0,0);
    Rect3d voxels(0,0,0,3,0,0);
    Vector3b nonSymmetricDimensions(true,true,true);
    DynamicRLE3<double> fillFactors(nonSymmetricDimensions.asArray());
    FillFactorRowVoxelizer rv(X, bounds, voxels, nonSymmetricDimensions, fillFactors);
    
    vector<FillRecord> fillRecords;
    fillRecords.push_back(FillRecord(0, -1, -0.5)); // e.g. downward-facing tri
    fillRecords.push_back(FillRecord(2, 1, 0.5)); // e.g. upward-facing tri
    
    rv.handleSortedFillRecords(0, 0, fillRecords);
    
    BOOST_CHECK_CLOSE(fillFactors(0,0,0), 0.5, 1e-4);
    BOOST_CHECK_CLOSE(fillFactors(1,0,0), 1.0, 1e-4);
    BOOST_CHECK_CLOSE(fillFactors(2,0,0), 0.5, 1e-4);
    BOOST_CHECK_SMALL(fillFactors.at(3,0,0), 1e-7);
}



BOOST_AUTO_TEST_CASE( OrientationVoxelizeOneTri )
{
    Rect3d bounds(0,0,0,1,1,1);
    Rect3d voxels(0,0,0,0,0,0);
    Vector3b nonSymmetricDimensions(true,true,true);
    DynamicRLE3<Matrix3d> fillFactors(nonSymmetricDimensions.asArray());
    OrientationRowVoxelizer rv(Z, bounds, voxels, nonSymmetricDimensions, fillFactors);
    
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    SimpleMesh::Triangle tri(SimpleMesh::Test::oneTri(o + z/2, x + z/2, y + z/2));
    ClippedTriangle clipt = rv.clipToRow(tri, 0, 0);
    
    vector<OrientationRecord> orientationRecords;
    rv.voxelizeOneTri(&clipt, 0, 0, orientationRecords);
    
    BOOST_CHECK_EQUAL(orientationRecords.size(), 1);
    BOOST_CHECK_EQUAL(orientationRecords[0].voxel, 0);
    BOOST_CHECK_EQUAL(orientationRecords[0].polygon, &clipt);
    orientationRecords.clear();
    
    SimpleMesh::Triangle tippedTri(SimpleMesh::Test::oneTri(o + 0.000000001*z,
        x + 0.999999999*z, y + 0.999999999*z) );
    clipt = rv.clipToRow(tippedTri, 0, 0);
    
    rv.voxelizeOneTri(&clipt, 0, 0, orientationRecords);
    
    BOOST_CHECK_EQUAL(orientationRecords.size(), 1);
    BOOST_CHECK_EQUAL(orientationRecords[0].voxel, 0);
    BOOST_CHECK_EQUAL(orientationRecords[0].polygon, &clipt);
    orientationRecords.clear();
    
//    SimpleMesh::Triangle tallTri(SimpleMesh::Test::oneTri(o - 0.999999*z,
//        x + 0.999999*z, y + 0.99999*z));
//    rv.voxelizeOneTri(&tallTri, 0, 0, orientationRecords);
//    
//    BOOST_CHECK_EQUAL(orientationRecords.size(), 2);
//    BOOST_CHECK_EQUAL(orientationRecords[0].voxel, -1);
//    BOOST_CHECK_CLOSE(orientationRecords[0].area, 1.0/8, 1e-2);
//    BOOST_CHECK_CLOSE(orientationRecords[0].volume, 2.0/8/3, 1e-2);
//    BOOST_CHECK_EQUAL(orientationRecords[1].voxel, 0);
//    BOOST_CHECK_CLOSE(orientationRecords[1].area, 0.5 - 0.125, 1e-2);
//    BOOST_CHECK_CLOSE(orientationRecords[1].volume, 5.0/24, 1e-2);
}




