/*
 *  testOrientation.cpp
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 7/22/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testTrackedPolyhedronSensitivity

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Voxelize.h"
#include "SimpleMesh/Search.h"

using namespace std;

static Matrix3d sFaceMatrix(int xyz)
{
    return outerProduct(Vector3d::unit(xyz), Vector3d::unit(xyz));
}

#define CHECK_CLOSE(mat1, mat2) BOOST_CHECK_SMALL(normInf((mat1)-(mat2)), 1e-10);

const int X = 0, Y = 1, Z = 2;

BOOST_AUTO_TEST_CASE( FaceMatrix )
{
    Matrix3d Mat00 = sFaceMatrix(0);
    Matrix3d Mat00_ref(Matrix3d::zero());
    Mat00_ref(0,0) = 1.0;
    
    BOOST_CHECK(Mat00 == Mat00_ref);
}


BOOST_AUTO_TEST_CASE( StandardBlock )
{
    Rect3d bounds(-3,-2,-1,3,2,1);
    Rect3i voxels(-1,-1,-1,1,1,1);
    Rect3d cubeBounds(-2, -4.0/3.0, -2.0/3.0, 2, 4.0/3, 2.0/3);
    
    std::vector<SimpleMesh::ControlVertex> controlVertices1, controlVertices2;
    std::vector<SimpleMesh::Triangle> meshTris1, meshTris2;
    SimpleMesh::Test::oneBoxMesh(controlVertices1, meshTris1, cubeBounds);
    
    Voxelize vox(bounds, voxels);
    
    RLE::DynamicRLE3<Matrix3d> ori = vox.orientations(meshTris1);
    
    // These matrices are the sum of polygon areas times orientation matrices.
    // I had to draw some diagrams to get this straight.  :-)
    
    Matrix3d cornerOrientation = 2./9*sFaceMatrix(X) +
        1./3*sFaceMatrix(Y) + 2./3*sFaceMatrix(Z);
    
    Matrix3d xyOrientation = 4./9*sFaceMatrix(X) + 2./3*sFaceMatrix(Y);
    Matrix3d yzOrientation = 2./3*sFaceMatrix(Y) + 4./3*sFaceMatrix(Z);
    Matrix3d zxOrientation = 4./3*sFaceMatrix(Z) + 4./9*sFaceMatrix(X);
    
    Matrix3d xOrientation = 8./9*sFaceMatrix(X);
    Matrix3d yOrientation = 4./3*sFaceMatrix(Y);
    Matrix3d zOrientation = 8./3*sFaceMatrix(Z);
    
    Matrix3d zed = Matrix3d::zero();
    
    for (int xx = -1; xx <= 1; xx += 2)
    for (int yy = -1; yy <= 1; yy += 2)
    for (int zz = -1; zz <= 1; zz += 2)
    {
        CHECK_CLOSE(ori.at(xx,yy,zz), cornerOrientation);
    }
    
    for (int xx = -1; xx <= 1; xx += 2)
    for (int yy = -1; yy <= 1; yy += 2)
    {
        CHECK_CLOSE(ori.at(xx,yy,0), xyOrientation);
    }
    
    for (int yy = -1; yy <= 1; yy += 2)
    for (int zz = -1; zz <= 1; zz += 2)
    {
        CHECK_CLOSE(ori.at(0,yy,zz), yzOrientation);
    }
    
    for (int zz = -1; zz <= 1; zz += 2)
    for (int xx = -1; xx <= 1; xx += 2)
    {
        CHECK_CLOSE(ori.at(xx,0,zz), zxOrientation);
    }
    
    for (int ww = -1; ww <= 1; ww += 2)
    {
        CHECK_CLOSE(ori.at(ww,0,0), xOrientation);
        CHECK_CLOSE(ori.at(0,ww,0), yOrientation);
        CHECK_CLOSE(ori.at(0,0,ww), zOrientation);
    }
    
    CHECK_CLOSE(ori.at(0,0,0), zed);
}





























