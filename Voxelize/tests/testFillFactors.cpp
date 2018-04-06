/*
 *  testTrackedPolyhedronSensitivity.cpp
 *
 *  Created by Paul C Hansen on 6/27/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testFillFactors

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Voxelize.h"
#include "SimpleMesh/Search.h"

using namespace std;


BOOST_AUTO_TEST_CASE( StandardCube )
{
    Rect3d rect1 = Rect3d(-1,-1,-1,1,1,1);
    std::vector<SimpleMesh::ControlVertex> controlVertices1, controlVertices2;
    std::vector<SimpleMesh::Triangle> meshTris1, meshTris2;
    SimpleMesh::Test::oneBoxMesh(controlVertices1, meshTris1, rect1);
    
    Rect3d bounds(-3,-2,-1,3,2,1);
    Rect3i voxels(-1,-1,-1,1,1,1);
    Voxelize vox(bounds, voxels);
    
    RLE::DynamicRLE3<double> ff = vox.fillFactors(meshTris1);
    
    Vector3d voxSize = vox.voxelSize();
    double cellVolume = voxSize[0]*voxSize[1]*voxSize[2];
    
    for (int xx = -1; xx <= 1; xx++)
    for (int yy = -1; yy <= 1; yy++)
    for (int zz = -1; zz <= 1; zz++)
    {
        double volume = cellVolume;
        
        if (xx == -1 || xx == 1)
            volume *= 0;
        
        if (yy == -1 || yy == 1)
            volume *= 0.25;
        
        if (volume != 0.0)
            BOOST_CHECK_CLOSE(ff.at(xx,yy,zz), volume, 1e-10);
        else
            BOOST_CHECK_SMALL(ff.at(xx,yy,zz), 1e-10);
    }
}




BOOST_AUTO_TEST_CASE( Cube2x2x2 )
{
    Rect3d rect1 = Rect3d(-1,-1,-1,1,1,1);
    std::vector<SimpleMesh::ControlVertex> controlVertices1, controlVertices2;
    std::vector<SimpleMesh::Triangle> meshTris1, meshTris2;
    SimpleMesh::Test::oneBoxMesh(controlVertices1, meshTris1, rect1);
    
    Rect3d bounds(-2,-2,-2,2,2,2);
    Rect3i voxels(0,0,0,3,3,3);
    Voxelize vox(bounds, voxels);
    
    RLE::DynamicRLE3<double> ff = vox.fillFactors(meshTris1);
    
    for (int xx = 0; xx < 4; xx++)
    for (int yy = 0; yy < 4; yy++)
    for (int zz = 0; zz < 4; zz++)
    if (xx >= 1 && xx <= 2 && yy >= 1 && yy <= 2 && zz >= 1 && zz <= 2)
    {
        BOOST_CHECK_CLOSE(ff.at(xx,yy,zz), 1.0, 1e-10);
    }
    else
    {
        BOOST_CHECK_SMALL(ff.at(xx,yy,zz), 1e-10);
    }
}

BOOST_AUTO_TEST_CASE( Cube1x1x1 )
{
    Rect3d rect1 = 0.5*Rect3d(-1,-1,-1,1,1,1);
    std::vector<SimpleMesh::ControlVertex> controlVertices1, controlVertices2;
    std::vector<SimpleMesh::Triangle> meshTris1, meshTris2;
    SimpleMesh::Test::oneBoxMesh(controlVertices1, meshTris1, rect1);
    
    Rect3d bounds(-2,-2,-2,2,2,2);
    Rect3i voxels(0,0,0,3,3,3);
    Voxelize vox(bounds, voxels);
    
    RLE::DynamicRLE3<double> ff = vox.fillFactors(meshTris1);
    
    for (int xx = 0; xx < 4; xx++)
    for (int yy = 0; yy < 4; yy++)
    for (int zz = 0; zz < 4; zz++)
    if (xx >= 1 && xx <= 2 && yy >= 1 && yy <= 2 && zz >= 1 && zz <= 2)
    {
        BOOST_CHECK_CLOSE(ff.at(xx,yy,zz), 0.125, 1e-10);
    }
    else
    {
        BOOST_CHECK_SMALL(ff.at(xx,yy,zz), 1e-10);
    }
}

























