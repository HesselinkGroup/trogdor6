/*
 *  testOrientationJacobians.cpp
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 8/14/11.
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

static vector<ControlVertex> sBlockVertices(Rect3d r, double angle = 0.0);

static vector<Vector3d> sSubtract(const vector<ControlVertex> & cv0,
    const vector<ControlVertex> & cv1, double dh = 1.0)
{
    vector<Vector3d> out(cv0.size());
    assert(cv0.size() == cv1.size());
    
    for (int nn = 0; nn < cv0.size(); nn++)
        out[nn] = (cv0[nn].point() - cv1[nn].point()) / dh;
    
    return out;
}

static Matrix3d sDot(Matrix3<Vector3d> lhs, Vector3d rhs)
{
    Matrix3d mat;
    for (int nn = 0; nn < 9; nn++)
        mat[nn] = dot(lhs[nn], rhs);
    
    return mat;
}


static void printJacobians(const map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > > & m);

static void sStandardBlock(int chosenTriangle);

BOOST_AUTO_TEST_CASE( StandardBlock )
{
    for (int tt = 0; tt < 12; tt++)
    {
//        cerr << "Chosen triangle " << tt << "\n";
        sStandardBlock(tt);
    }
    
    sStandardBlock(-1);
}

static void sStandardBlock(int chosenTriangle)
{
    Rect3d bounds(-3,-2,-1,3,2,1);
    Rect3i voxels(-1,-1,-1,1,1,1);
    Voxelize vox(bounds, voxels);
    
    Rect3d rect0(-2, -4.0/3.0, -2.0/3.0, 2, 4.0/3, 2.0/3);
    
    std::vector<SimpleMesh::ControlVertex> controlVertices1, controlVertices2;
    std::vector<SimpleMesh::Triangle> meshTris1, meshTris2;
    SimpleMesh::Test::oneBoxMesh(controlVertices1, meshTris1, rect0);
    
    if (chosenTriangle != -1)
        meshTris1 = std::vector<SimpleMesh::Triangle>(1, meshTris1.at(chosenTriangle));
    
    DynamicRLE3<Matrix3d> orientations0 = vox.orientations(meshTris1);
    map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > > DNNT =
        vox.orientationJacobians(meshTris1);
    
//    cerr << block0 << "\n";
//    cerr << "Orientations #0:\n" << orientations0 << "\n";
//    
//    printJacobians(DNNT);
    
//    vector<ControlVertex> cvs = sBlockVertices(rect0);
//    copy(cvs.begin(), cvs.end(), ostream_iterator<ControlVertex>(cerr, "\n"));
    
    double epsilon = 1e-7;
    for (int translateXYZ = 0; translateXYZ < 3; translateXYZ++)
    {
//        cerr << "Translate " << char('x'+translateXYZ) << "\n";
        Rect3d rect1 = rect0 + epsilon*Vector3d::unit(translateXYZ);
        
        SimpleMesh::Test::oneBoxMesh(controlVertices2, meshTris2, rect1);
        if (chosenTriangle != -1)
            meshTris2 = std::vector<SimpleMesh::Triangle>(1, meshTris2.at(chosenTriangle));
        
        DynamicRLE3<Matrix3d> orientations1 = vox.orientations(meshTris2);
//        cerr << "Orientations #1:\n" << orientations1 << "\n";
        
        DynamicRLE3<Matrix3d> deltaOr = (orientations1 - orientations0);
        vector<Vector3d> deltaCV = SimpleMesh::Test::subtractPositions(controlVertices2, controlVertices1);
        
        // Test the jacobian now.
        DynamicRLE3<Matrix3d> orPredicted(orientations0);
        
        for (int vv = 0; vv < deltaCV.size(); vv++)
        {
            DynamicRLE3<Matrix3d> dotProd;
            RLE::scalarTransform(DNNT[vv], deltaCV[vv], dotProd, sDot);
            orPredicted = orPredicted + dotProd;
//            cerr << "v" << vv << ": dot =\n" << dotProd << "\n";
//            cerr << "orPredicted is " << orPredicted << "\n";
        }
        
        //cerr << "Prediction:\n" << ffPredicted << "\n";
        
        DynamicRLE3<Matrix3d> deltaOrMeas = (orientations1-orientations0);
        DynamicRLE3<Matrix3d> deltaOrPred = (orPredicted - orientations0);
        
//        cerr << "Derivative measured:\n" << dOrMeas;
//        cerr << "\nDerivative expected:\n" << dOrPred
//            << "\n";
        
        bool failedAnywhere = false;
        for (int zz = voxels.p1[2]; zz <= voxels.p2[2]; zz++)
        for (int yy = voxels.p1[1]; yy <= voxels.p2[1]; yy++)
        for (int xx = voxels.p1[0]; xx <= voxels.p2[0]; xx++)
        {
//            cerr << "xyz " << xx << ", " << yy << ", " << zz << "\n";
            bool failed = 0;
            for (int ii = 0; ii < 3; ii++)
            for (int jj = 0; jj < 3; jj++)
            {
//                BOOST_CHECK_SMALL(normInf(
//                    dOrPred.markAt(xx,yy,zz, Matrix3d()) -
//                    dOrMeas.markAt(xx,yy,zz, Matrix3d())), 1e-5);
                
                if (normInf(
                    deltaOrPred.markAt(xx,yy,zz, Matrix3d()) -
                    deltaOrMeas.markAt(xx,yy,zz, Matrix3d())) > 1e-5)
                {
                    failed = true;
                    failedAnywhere = true;
                }
            }
            
            if (failed)
            {
//                cerr << "FAIL: tri " << chosenTriangle << " @" << Vector3i(xx,yy,zz) << " translate "
//                    << char('x'+translateXYZ) << "\n";
//
//                cerr << "quickPatch(" << block0.triangles()[chosenTriangle]
//                    << ", 'g');\n";
//                cerr << "hold on;\n";
//
                Rect3d thisVoxel;
                
                thisVoxel = vox.voxel(xx, yy ,zz);
//                cerr << "quickRect(" << thisVoxel << ", 'r');\n";
//
//                cerr << "Predicted derivative:\n"
//                    << dOrPred.markAt(xx,yy,zz, Matrix3d()) << "\n"
//                    << "Measured derivative:\n"
//                    << dOrMeas.markAt(xx,yy,zz, Matrix3d()) << "\n";
            }
        }
        BOOST_CHECK(failedAnywhere == false);
    }
}

BOOST_AUTO_TEST_CASE( RotateBlock )
{
//    Rect3d bounds(-2, -2, -2, 2, 2, 2);
//    Rect3i voxels(-1,-1,-1,1,1,1);
//    Rect3i voxels(0,0,0,1,1,0);
    Rect3d bounds(0,-2,-2,2,0,0);
    Rect3i voxels(0,0,0,0,0,0);
    Voxelize vox(bounds, voxels);
    
    Rect3d rect0(-1, -1, -1, 1, 1, 1);
    
    std::vector<SimpleMesh::ControlVertex> controlVertices1, controlVertices2;
    std::vector<SimpleMesh::Triangle> meshTris1, meshTris2;
    SimpleMesh::Test::oneBoxMesh(controlVertices1, meshTris1, rect0);
    
//    TrackedPolyhedron block0 = sBlock(rect0);
    
    DynamicRLE3<Matrix3d> orientations0 = vox.orientations(meshTris1);
    map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > > DNNT =
        vox.orientationJacobians(meshTris1);
        
//    cerr << "Orientations #0:\n" << orientations0 << "\n";
//    
//    printJacobians(DNNT);
    
//    vector<ControlVertex> cvs = sBlockVertices(rect0);
//    copy(cvs.begin(), cvs.end(), ostream_iterator<ControlVertex>(cerr, "\n"));
    
    double epsilon = 1e-8;
    
    {
        double theta = epsilon*2*M_PI;
        
        SimpleMesh::Test::oneBoxMesh(controlVertices2, meshTris2, rect0, Matrix3d::rotation(theta, Vector3d(0,0,1)));
//        TrackedPolyhedron block1 = sBlock(rect0, theta);
        BOOST_CHECK(true);
        DynamicRLE3<Matrix3d> orientations1 = vox.orientations(meshTris2);

//        cerr << "Orientations #1:\n" << orientations1 << "\n";
        BOOST_CHECK(true);
        DynamicRLE3<Matrix3d> deltaOr = (orientations1 - orientations0);
        
        vector<Vector3d> deltaCV = SimpleMesh::Test::subtractPositions(controlVertices2, controlVertices1);
        
//        std::copy(deltaCV.begin(), deltaCV.end(), ostream_iterator<Vector3d>(cerr, "\n"));
        
        // Test the jacobian now.
        DynamicRLE3<Matrix3d> orPredicted(orientations0);
        
//        cerr << "orPredicted BEGINS as " << orPredicted << "\n";
        
        for (int vv = 0; vv < deltaCV.size(); vv++)
        if (DNNT.count(vv))
        {
            BOOST_CHECK(true);
            DynamicRLE3<Matrix3d> dotProd;
            BOOST_CHECK(true);
            RLE::scalarTransform(DNNT[vv], deltaCV[vv], dotProd, sDot);
            BOOST_CHECK(true);
            orPredicted = orPredicted + dotProd;
            BOOST_CHECK(true);
//            cerr << "v" << vv << ": dot =\n" << dotProd << "\n";
//            cerr << "orPredicted is " << orPredicted << "\n";
        }
        
        //cerr << "Prediction:\n" << ffPredicted << "\n";
        BOOST_CHECK(true);
        DynamicRLE3<Matrix3d> dOrMeas = (orientations1-orientations0)/epsilon;
        BOOST_CHECK(true);
        DynamicRLE3<Matrix3d> dOrPred = (orPredicted - orientations0)/epsilon;
        
//        cerr << "Derivative measured:\n" << dOrMeas;
//        cerr << "\nDerivative expected:\n" << dOrPred
//            << "\n";
        
        for (int zz = voxels.p1[2]; zz <= voxels.p2[2]; zz++)
        for (int yy = voxels.p1[1]; yy <= voxels.p2[1]; yy++)
        for (int xx = voxels.p1[0]; xx <= voxels.p2[0]; xx++)
        {
            for (int ii = 0; ii < 3; ii++)
            for (int jj = 0; jj < 3; jj++)
            {
                BOOST_CHECK_SMALL(normInf(
                    dOrPred.markAt(xx,yy,zz, Matrix3d()) -
                    dOrMeas.markAt(xx,yy,zz, Matrix3d())), 1e-5);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( RotateTriangle )
{
//    cerr << "\n\n*** One tri ***\n\n";
//    Rect3d bounds(-2, -2, -2, 2, 2, 2);
//    Rect3i voxels(-1,-1,-1,1,1,1);
//    Rect3i voxels(0,0,0,1,1,0);
//    Rect3d bounds(0,0,0,2,2,2);
//    Rect3i voxels(0,0,0,0,0,0);
    
    Rect3d bounds(0,-2,-2,2,0,0);
    Rect3i voxels(0,0,0,0,0,0);
    Voxelize vox(bounds, voxels);
    
//    Triangle3d tri0(Vector3d(1, 1, -1), Vector3d(1, 1, 3),
//        Vector3d(1, -3, -1));
    
//    Triangle3d tri0(Vector3d(1,3,-1), Vector3d(1,-1,3), Vector3d(1,-1,-1));
    Triangle3d tri0(Vector3d(1,3,-3), Vector3d(1,-1,-3), Vector3d(1,-1,1));
    
    vector<Triangle> tris0;
    tris0.push_back(SimpleMesh::Test::oneTri(tri0[0], tri0[1], tri0[2]));
    
    BOOST_CHECK(true);
    DynamicRLE3<Matrix3d> orientations0 = vox.orientations(tris0);
    map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > > DNNT =
        vox.orientationJacobians(tris0);
        
//    cerr << "Orientations #0:\n" << orientations0 << "\n";
    
//    printJacobians(DNNT);
    
//    vector<ControlVertex> cvs = sBlockVertices(rect0);
//    copy(cvs.begin(), cvs.end(), ostream_iterator<ControlVertex>(cerr, "\n"));
    
    double epsilon = 1e-8;
    
    {
        double theta = epsilon*2*M_PI;
        Matrix3d rot(cos(theta), -sin(theta), 0.0,
            sin(theta), cos(theta), 0.0,
            0.0, 0.0, 1.0);
        
        Triangle3d tri1 = rot*tri0;
        vector<Triangle> tris1;
        tris1.push_back(SimpleMesh::Test::oneTri(tri1[0], tri1[1], tri1[2]));
        
        DynamicRLE3<Matrix3d> orientations1 = vox.orientations(tris1);

//        cerr << "Orientations #1:\n" << orientations1 << "\n";
        
        DynamicRLE3<Matrix3d> deltaOr = (orientations1 - orientations0)/epsilon;
        
        BOOST_CHECK(true);
        vector<Vector3d> deltaCV(3);
        for (int nn = 0; nn < 3; nn++)
            deltaCV[nn] = tri1[nn] - tri0[nn];
        
//        cerr << "Change in CVs:\n";
//        std::copy(deltaCV.begin(), deltaCV.end(), ostream_iterator<Vector3d>(cerr, "\n"));
        
        // Test the jacobian now.
        DynamicRLE3<Matrix3d> orPredicted(orientations0);
        
        BOOST_CHECK(true);
//        cerr << "orPredicted BEGINS as " << orPredicted << "\n";
        
        for (int vv = 0; vv < deltaCV.size(); vv++)
        {
            BOOST_CHECK(true);
            DynamicRLE3<Matrix3d> dotProd;
            BOOST_CHECK(true);
            RLE::scalarTransform(DNNT[vv], deltaCV[vv], dotProd, sDot);
            BOOST_CHECK(true);
            orPredicted = orPredicted + dotProd;
            BOOST_CHECK(true);
//            cerr << "v" << vv << ": dot =\n" << dotProd << "\n";
//            cerr << "orPredicted is " << orPredicted << "\n";
        }
        
        //cerr << "Prediction:\n" << ffPredicted << "\n";
        
        BOOST_CHECK(true);
        DynamicRLE3<Matrix3d> dOrMeas = (orientations1-orientations0)/epsilon;
        
        BOOST_CHECK(true);
        DynamicRLE3<Matrix3d> dOrPred = (orPredicted - orientations0)/epsilon;
        
//        cerr << "Derivative measured:\n" << dOrMeas;
//        cerr << "\nDerivative expected:\n" << dOrPred
//            << "\n";
        
        for (int zz = voxels.p1[2]; zz <= voxels.p2[2]; zz++)
        for (int yy = voxels.p1[1]; yy <= voxels.p2[1]; yy++)
        for (int xx = voxels.p1[0]; xx <= voxels.p2[0]; xx++)
        {
            for (int ii = 0; ii < 3; ii++)
            for (int jj = 0; jj < 3; jj++)
            {
                BOOST_CHECK_SMALL(normInf(
                    dOrPred.markAt(xx,yy,zz, Matrix3d()) -
                    dOrMeas.markAt(xx,yy,zz, Matrix3d())), 1e-5);
            }
        }
    }
}


BOOST_AUTO_TEST_CASE( OneTriangle )
{
//    return;
//    cerr << "OneTriangle\n";
    
    Rect3d bounds(0,0,0,1,1,1);
//    Rect3i voxels(0,0,0,0,0,0); // one voxel
    Rect3i voxels(0,0,0,1,0,0); // split X
//    Rect3i voxels(0,0,0,0,1,0); // split Y
//    Rect3i voxels(0,0,0,0,0,1); // split Z
    Voxelize vox(bounds, voxels);
    
//    Triangle3d t0(Vector3d(-0.5, 0.5, 0.5),
//        Vector3d(0.5, 0.5, 0.5), Vector3d(-0.5, 1.5, 0.5));
//    Triangle3d t0(Vector3d(0.1, 0.1, 0.1),
//        Vector3d(0.1, 0.9, 0.1), Vector3d(0.9, 0.1, 0.1)); // downward
    Triangle3d t0(Vector3d(0.1, 0.9, 0.1),
        Vector3d(0.9, 0.9, 0.1), Vector3d(0.9, 0.1, 0.1));
    
    vector<Triangle> triangles0;
    triangles0.push_back(SimpleMesh::Test::oneTri(t0[0], t0[1], t0[2]));
    
    DynamicRLE3<Matrix3d> orientations0 = vox.orientations(triangles0);
    map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > > DNNT =
        vox.orientationJacobians(triangles0);
        
//    cerr << "Orientations #0:\n" << orientations0 << "\n";
    
//    printJacobians(DNNT);
    
    double epsilon = 1e-7;
    for (int translateXYZ = 0; translateXYZ < 3; translateXYZ++)
    {
//        cerr << "Translate " << char('x'+translateXYZ) << "\n";
        
        Triangle3d t1 = t0;
        for (int mm = 0; mm < 3; mm++)
            t1[mm] += epsilon*Vector3d::unit(translateXYZ);
        
        vector<Triangle> triangles1;
        triangles1.push_back(SimpleMesh::Test::oneTri(t1[0], t1[1], t1[2]));
        
        DynamicRLE3<Matrix3d> orientations1 = vox.orientations(triangles1);
        //vox.fillFactorJacobians(block1.triangles());
        
//        cerr << "Orientations #1:\n" << orientations1 << "\n";
        
        DynamicRLE3<Matrix3d> deltaOr = (orientations1 - orientations0)/epsilon;
        
        vector<Vector3d> deltaCV(3);
        for (int nn = 0; nn < 3; nn++)
            deltaCV[nn] = t1[nn] - t0[nn];
        
        //std::copy(deltaCV.begin(), deltaCV.end(), ostream_iterator<Vector3d>(cerr, "\n"));
        
        // Test the jacobian now.
        DynamicRLE3<Matrix3d> orPredicted(orientations0);
        
//        cerr << "orPredicted BEGINS as " << orPredicted << "\n";
        
        for (int vv = 0; vv < deltaCV.size(); vv++)
        {
            DynamicRLE3<Matrix3d> dotProd;
            RLE::scalarTransform(DNNT[vv], deltaCV[vv], dotProd, sDot);
            orPredicted = orPredicted + dotProd;
//            cerr << "v" << vv << ": dot =\n" << dotProd << "\n";
//            cerr << "orPredicted is " << orPredicted << "\n";
        }
        
        //cerr << "Prediction:\n" << ffPredicted << "\n";
        
        DynamicRLE3<Matrix3d> dOrMeas = (orientations1-orientations0)/epsilon;
        DynamicRLE3<Matrix3d> dOrPred = (orPredicted - orientations0)/epsilon;
        
//        cerr << "Derivative measured:\n" << dOrMeas;
//        cerr << "\nDerivative expected:\n" << dOrPred
//            << "\n";
        

        for (int zz = voxels.p1[2]; zz <= voxels.p2[2]; zz++)
        for (int yy = voxels.p1[1]; yy <= voxels.p2[1]; yy++)
        for (int xx = voxels.p1[0]; xx <= voxels.p2[0]; xx++)
        {
            bool problem = 0;
            
            for (int ii = 0; ii < 3; ii++)
            for (int jj = 0; jj < 3; jj++)
            {
//                BOOST_CHECK_SMALL(normInf(
//                    dOrPred.markAt(xx,yy,zz, Matrix3d()) -
//                    dOrMeas.markAt(xx,yy,zz, Matrix3d())), 1e-5);
                
                if (normInf(dOrPred.markAt(xx,yy,zz, Matrix3d()) -
                    dOrMeas.markAt(xx,yy,zz, Matrix3d())) > 1e-5)
                    problem = true;
            }
            if (problem)
            {
                cerr << "Triangle (" << triangles0[0] << " normal "
                    << triangles0[0].unitNormal() << ") perturbed along "
                    << char('x'+translateXYZ) << " is wrong in cell "
                    << Vector3i(xx,yy,zz) << ".\n";
            }
        }
    }
}



//static TrackedPolyhedron sBlock(Rect3d r, double angle)
//{
//    Matrix3d rot(cos(angle), -sin(angle), 0.0,
//        sin(angle), cos(angle), 0.0,
//        0.0, 0.0, 1.0);
//
//    NefBuilder nef = NefBuilder::block(r, 0, rot);
//
//    TrackedPolyhedron tp(nef.polyhedron(), nef.getControlVertexMap());
////    tp.cacheSensitivity();
//
//    return tp;
//}
//
//static vector<ControlVertex> sBlockVertices(Rect3d r, double angle)
//{
//    Matrix3d rot(cos(angle), -sin(angle), 0.0,
//        sin(angle), cos(angle), 0.0,
//        0.0, 0.0, 1.0);
//    NefBuilder nef = NefBuilder::block(r, 0, rot);
//
//    return nef.controlVertices();
//}

static void printJacobians(const map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > > & m)
{
    map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > >::const_iterator itr;
    
    for (itr = m.begin(); itr != m.end(); itr++)
    {
        cout << "CV" << itr->first << ": \n";
        cout << itr->second << "\n";
    }
}
