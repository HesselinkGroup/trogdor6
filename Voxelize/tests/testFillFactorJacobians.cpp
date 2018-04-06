/*
 *  testFillFactorJacobians.cpp
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 8/13/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testFillFactorJacobians

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Voxelize.h"
#include "SimpleMesh/Search.h"

using namespace std;
using namespace RLE;
using namespace SimpleMesh;

static vector<ControlVertex> sBlockVertices(Rect3d r);

static vector<Vector3d> sSubtract(const vector<ControlVertex> & cv0,
    const vector<ControlVertex> & cv1, double dh = 1.0)
{
    vector<Vector3d> out(cv0.size());
    assert(cv0.size() == cv1.size());
    
    for (int nn = 0; nn < cv0.size(); nn++)
        out[nn] = (cv0[nn].point() - cv1[nn].point()) / dh;
    
    return out;
}

// used in an RLE::transform
static double sDot(Vector3d lhs, Vector3d rhs)
{
//    if (isnan(dot(lhs, rhs)))
//        int blarg = 54;
    return dot(lhs, rhs);
}

static void printJacobians(const map<unsigned int, DynamicRLE3<Vector3d> > & m);

BOOST_AUTO_TEST_CASE( StandardBlock )
{
    Rect3d bounds(-3,-2,-1,3,2,1);
    Rect3i voxels(-1,-1,-1,1,1,1);
    Voxelize vox(bounds, voxels);
    
    Rect3d rect1(-2, -4.0/3.0, -2.0/3.0, 2, 4.0/3, 2.0/3);
    
    std::vector<SimpleMesh::ControlVertex> controlVertices1, controlVertices2;
    std::vector<SimpleMesh::Triangle> meshTris1, meshTris2;
    SimpleMesh::Test::oneBoxMesh(controlVertices1, meshTris1, rect1);
    
    DynamicRLE3<double> fillFactors0 = vox.fillFactors(meshTris1);
    map<unsigned int, DynamicRLE3<Vector3d> > DF = vox.fillFactorJacobians(meshTris1);
        
//    cerr << "Fill factors #0:\n" << fillFactors0 << "\n";
    
    double epsilon = 1e-6;
    for (int translateXYZ = 0; translateXYZ < 3; translateXYZ++)
    {
        Rect3d rect2 = rect1 + epsilon*Vector3d::unit(translateXYZ);
        SimpleMesh::Test::oneBoxMesh(controlVertices2, meshTris2, rect2);
    
        DynamicRLE3<double> fillFactors1 = vox.fillFactors(meshTris2);
        DynamicRLE3<double> deltaFF = (fillFactors1 - fillFactors0);
        vector<Vector3d> deltaCV = SimpleMesh::Test::subtractPositions(controlVertices2, controlVertices1);
        
        // Test the jacobian now.
        DynamicRLE3<double> ffPredicted(fillFactors0);
        
        for (int vv = 0; vv < deltaCV.size(); vv++)
        {
            DynamicRLE3<double> dotProd;
            RLE::scalarTransform(DF[vv], deltaCV[vv], dotProd, sDot);
            ffPredicted = ffPredicted + dotProd;
            
//            cerr << "v" << vv << ": dot =\n" << dotProd << "\n";
        }
        
//        cerr << "Prediction:\n" << ffPredicted << "\n";
        
//        cerr << "Derivative measured:\n" << (fillFactors1 - fillFactors0)/epsilon;
//        cerr << "\nDerivative expected:\n" << (ffPredicted - fillFactors0)/epsilon
//            << "\n";
        
        DynamicRLE3<double> dfPredicted = (ffPredicted - fillFactors0);
        
        for (int xx = voxels.p1[0]; xx <= voxels.p2[0]; xx++)
        for (int yy = voxels.p1[1]; yy <= voxels.p2[1]; yy++)
        for (int zz = voxels.p1[2]; zz <= voxels.p2[2]; zz++)
        {
            if (dfPredicted.markAt(xx,yy,zz,0.0) != 0.0)
            {
                BOOST_CHECK_CLOSE(dfPredicted.markAt(xx,yy,zz,0.0), deltaFF.markAt(xx,yy,zz,0.0), 1e-4);
            }
            else
            {
                BOOST_CHECK_SMALL(deltaFF.markAt(xx,yy,zz,0.0), 1e-6);
            }
        }
    }
}


BOOST_AUTO_TEST_CASE( HalfFilled )
{
    Rect3d bounds(0,0,0,1,1,1);
//    Rect3i voxels(0,0,0,0,0,0); // one voxel
//    Rect3i voxels(0,0,0,1,0,0); // split X
//    Rect3i voxels(0,0,0,0,1,0); // split Y
//    Rect3i voxels(0,0,-1,0,0,1); // split Z
    Rect3i voxels(-1,-1,-1,1,1,1);
//    Rect3i voxels(0,0,0,1,1,1); 
    Voxelize vox(bounds, voxels);
    
    Rect3d rect0(0.1, 0.1, 0.1, 0.9, 0.9, 0.9); // inset XYZ
//    Rect3d rect0(0.1, 0.1, 0.1, 0.9, 0.9, 0.9);


    std::vector<SimpleMesh::ControlVertex> controlVertices1, controlVertices2;
    std::vector<SimpleMesh::Triangle> meshTris1, meshTris2;
    SimpleMesh::Test::oneBoxMesh(controlVertices1, meshTris1, rect0);
//
//    TrackedPolyhedron block0 = sBlock(rect0);
    
    DynamicRLE3<double> fillFactors0 = vox.fillFactors(meshTris1);
    map<unsigned int, DynamicRLE3<Vector3d> > DF = vox.fillFactorJacobians(meshTris1);
    
    double epsilon = 1e-3;
    for (int translateXYZ = 0; translateXYZ < 3; translateXYZ++)
    {
//        cerr << "Translate " << char('x'+translateXYZ) << "\n";
        Rect3d rect1 = rect0 + epsilon*Vector3d::unit(translateXYZ);
        
//        TrackedPolyhedron block1 = sBlock(rect1);
        SimpleMesh::Test::oneBoxMesh(controlVertices2, meshTris2, rect1);
        DynamicRLE3<double> fillFactors1 = vox.fillFactors(meshTris2);
        vox.fillFactorJacobians(meshTris2);
        
        //cerr << "Fill factors #1:\n" << fillFactors1 << "\n";
        
        DynamicRLE3<double> deltaFF = (fillFactors1 - fillFactors0);
        
        vector<Vector3d> deltaCV = SimpleMesh::Test::subtractPositions(controlVertices2, controlVertices1);
        
        // Test the jacobian now.
        DynamicRLE3<double> ffPredicted(fillFactors0);
        
        for (int vv = 0; vv < deltaCV.size(); vv++)
        {
            DynamicRLE3<double> dotProd;
            RLE::scalarTransform(DF[vv], deltaCV[vv], dotProd, sDot);
            ffPredicted = ffPredicted + dotProd;
//            cerr << "v" << vv << ": dot =\n" << dotProd << "\n";
        }
            
        DynamicRLE3<double> dfPredicted = (ffPredicted - fillFactors0);
        
        for (int xx = voxels.p1[0]; xx <= voxels.p2[0]; xx++)
        for (int yy = voxels.p1[1]; yy <= voxels.p2[1]; yy++)
        for (int zz = voxels.p1[2]; zz <= voxels.p2[2]; zz++)
        {
            if (dfPredicted.markAt(xx,yy,zz,0.0) != 0.0)
            {
                BOOST_CHECK_CLOSE(dfPredicted.markAt(xx,yy,zz,0.0), deltaFF.markAt(xx,yy,zz,0.0), 1e-4);
            }
            else
            {
                BOOST_CHECK_SMALL(deltaFF.markAt(xx,yy,zz,0.0), 1e-6);
            }
        }
    }
    
}



static void printJacobians(const map<unsigned int, DynamicRLE3<Vector3d> > & m)
{
    map<unsigned int, DynamicRLE3<Vector3d> >::const_iterator itr;
    
    for (itr = m.begin(); itr != m.end(); itr++)
    {
        cout << "CV" << itr->first << ": \n";
        cout << itr->second << "\n";
    }
}



