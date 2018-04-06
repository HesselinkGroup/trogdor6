/*
 *  testClipTriangleJacobians.cpp
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 8/11/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */


#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testClipTriangleJacobians

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Voxelize/ClipTriangle.h"
#include "SimpleMesh/Search.h"
#include "TimeWrapper.h"

using namespace std;
using namespace SimpleMesh;

//static ClippedTriangle sClip(const ClippedTriangle & ctri, Rect2d xyBounds);
static Matrix3d sDelta(const map<unsigned int, Matrix3<Vector3d> > & jacobian,
    Vector3d delta);

//#define CHECK_CLOSE(p0, p1, tol) BOOST_CHECK_SMALL(norm((p0)-(p1)), tol)
#define CHECK_CLOSE(p0, p1) BOOST_CHECK_SMALL(norm((p0)-(p1)), 1e-9)
#define CHECK_MATRIX_CLOSE(m0, m1) BOOST_CHECK_SMALL(vectorNorm((m0)-(m1)), 1e-9);

const void* NULLPTR = 0L;
const int X = 0, Y = 1, Z = 2;
//
//void printJacobian(const map<unsigned int, Vector3d> & dFdv)
//{
//    map<unsigned int, Vector3d>::const_iterator itr;
//    for (itr = dFdv.begin(); itr != dFdv.end(); itr++)
//    {
//        cerr << itr->first << ": " << itr->second << "\n";
//    }
//}
//void printJacobian(const map<unsigned int, Matrix3<Vector3d> > & dFdv)
//{
//    map<unsigned int, Matrix3<Vector3d> >::const_iterator itr;
//    for (itr = dFdv.begin(); itr != dFdv.end(); itr++)
//    {
//        cerr << itr->first << ": " << itr->second << "\n";
//    }
//}

BOOST_AUTO_TEST_CASE( IntegrateDz_Triangle )
{
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1), xy(1,1,0);
    Triangle xyTriangle0(SimpleMesh::Test::oneTri(o + z/2, x + z/2, y + z/2));
    ClippedTriangle ctri0(xyTriangle0);
    
//    double vol0 = ctri0.volume(Z);
    double area0 = ctri0.area2d(Z);
    
    map<unsigned int, Vector3d> jacobian;
    ctri0.intDz(jacobian);
    
    for (int vv = 0; vv < 3; vv++)
    {
        Vector3d dFdv = jacobian[vv];
        BOOST_CHECK_EQUAL(dFdv[0], 0.0);
        BOOST_CHECK_EQUAL(dFdv[1], 0.0);
        BOOST_CHECK_CLOSE(dFdv[2], area0/3.0, 1e-7);
    }
    
    /*
    cerr << "Initial volume is " << vol0 << "\n";
    
    Vector3d dz = 1e-6*z;
    Triangle xyTriangle1(SimpleMesh::Test::oneTri(o + z/2 + dz, x + z/2, y + z/2));
    ClippedTriangle ctri1(xyTriangle1);
    
    double vol1 = ctri1.volume(Z);
    
    cerr << "New volume is " << vol1 << "\n";
    cerr << "deltaV = " << vol1 - vol0 << "\n";
    
    cerr << "Jacobian has " << jacobian.size() << " elements.\n";
    printJacobian(jacobian);
    */
}

BOOST_AUTO_TEST_CASE( IntegrateDz_ClippedTriangle2D )
{
//    cerr << "SKIPPING\n"; return;
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1), xy(1,1,0);
    Triangle xyTriangle0(SimpleMesh::Test::oneTri(o + z/2, x + z/2, y + z/2));
    
    Rect2d clipBounds(0.2, 0.1, 0.9, 0.8);
//    ClippedTriangle ctri0 = sClip(ClippedTriangle(xyTriangle0), clipBounds);
    ClippedTriangle ctri0 = ClippedTriangle(xyTriangle0).partInside(clipBounds);
    
    double vol0 = ctri0.volume(Z);
//    double area0 = ctri0.area2d(Z);
    
    map<unsigned int, Vector3d> jacobian;
    ctri0.intDz(jacobian);
    
//    cerr << "Initial volume is " << vol0 << "\n";
    
    Vector3d dz = 1e-6*z;
    for (int vv = 0; vv < 3; vv++)
    {
        Triangle3d tri(o + z/2, x + z/2, y + z/2);
        tri[vv] += dz;
        Triangle xyTriangle1(SimpleMesh::Test::oneTri(tri[0], tri[1], tri[2]));
        ClippedTriangle ctri1 = ClippedTriangle(xyTriangle1).partInside(clipBounds);
        
        double vol1 = ctri1.volume(Z);
        
        double dVdv = (vol1-vol0)/norm(dz);
        
//        cerr << "dVdv = " << dVdv << "\n";
        BOOST_CHECK_EQUAL(jacobian[vv][X], 0.0);
        BOOST_CHECK_EQUAL(jacobian[vv][Y], 0.0);
        BOOST_CHECK_CLOSE(jacobian[vv][Z], dVdv, 1e-7);
    }
    
//    printJacobian(jacobian);
}

BOOST_AUTO_TEST_CASE( IntegrateZDhTriangle )
{
//    cerr << "SKIPPING\n"; return;
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1), xy(1,1,0);
    Triangle xyTriangle0(SimpleMesh::Test::oneTri(o + z/2, x + z/2, y + z/2));
    ClippedTriangle ctri0(xyTriangle0);
    
//    double vol0 = ctri0.volume(Z);
//    double area0 = ctri0.area2d(Z);
    
    map<unsigned int, Vector3d> jacobian;
    ctri0.intZDh(jacobian, Z);
    
    //printJacobian(jacobian);
    
    // It's pretty easy to work out these answers on paper.
    
    Vector3d dFdv0(-0.25, -0.25, 0.0);
    Vector3d dFdv1(0.25, 0, 0);
    Vector3d dFdv2(0.0, 0.25, 0.0);
    
    CHECK_CLOSE(jacobian[0], dFdv0);
    CHECK_CLOSE(jacobian[1], dFdv1);
    CHECK_CLOSE(jacobian[2], dFdv2);
}

BOOST_AUTO_TEST_CASE( IntegrateZDhTriangle_ClipXY )
{
//    cerr << "SKIPPING\n"; return;
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1), xy(1,1,0);
    Triangle3d tri0(o+z/2, x+z/2, y+z/2);
    Triangle xyTriangle0(SimpleMesh::Test::oneTri(tri0[0], tri0[1], tri0[2]));
    
    Rect2d clipBounds(0.2, 0.1, 0.9, 0.8);
    ClippedTriangle ctri0 = ClippedTriangle(xyTriangle0).partInside(clipBounds);
    
    double vol0 = ctri0.volume(Z);
    double area0 = ctri0.area2d(Z);
    
    map<unsigned int, Vector3d> jacobian;
    ctri0.intZDh(jacobian, Z);
    
    
    double eps = 1e-8; // roughly the most accurate point?  eh.
    for (int vv = 0; vv < 3; vv++)
    for (int xy = 0; xy < 2; xy++)
    {
        Triangle3d tri1(tri0);
        tri1[vv][xy] += eps;
        
        Triangle xyTriangle1(SimpleMesh::Test::oneTri(tri1[0], tri1[1], tri1[2]));
        ClippedTriangle ctri1 = ClippedTriangle(xyTriangle1).partInside(clipBounds);
        
        double vol1 = ctri1.volume(Z);
        
        double dFdv = (vol1-vol0)/eps;
        
        if (jacobian[vv][xy] == 0.0)
            BOOST_CHECK_SMALL(dFdv, 1e-3);
        else
            BOOST_CHECK_CLOSE(jacobian[vv][xy], dFdv, 1e-3);
    }
}

BOOST_AUTO_TEST_CASE( ClipZ )
{
//    cerr << "SKIPPING\n"; return;
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1), xy(1,1,0);
    Triangle3d tri0(o, x+z, y);
    Triangle xyTriangle0(SimpleMesh::Test::oneTri(tri0[0], tri0[1], tri0[2]));
    
    Plane3d planeZ(Vector3d::unit(Z), -0.5); // z = 0.5
    
    SensitiveEdge newEdge(xyTriangle0, planeZ);
    
    ClippedTriangle ctri0(xyTriangle0);
    ctri0 = ctri0.partAbove(Z, 0.5, &newEdge);
    
    double vol0 = ctri0.volume(Z);
    double area0 = ctri0.area2d(Z);
    
    map<unsigned int, Vector3d> interiorJacobian, totalJacobian;
//    ctri0.intZDh(jacobian);
    ctri0.intDz(interiorJacobian);
    
    // Yes, these are magic numbers.  It took me two hours or so to derive them.
    // Just be careful and you'll be ok.  :-D
    BOOST_CHECK_CLOSE(interiorJacobian[0][0], -1.0/48, 1e-6);
    BOOST_CHECK_SMALL(interiorJacobian[0][1], 1e-6);
    BOOST_CHECK_CLOSE(interiorJacobian[0][2], 1.0/48, 1e-6);
    BOOST_CHECK_CLOSE(interiorJacobian[1][0], -1.0/12, 1e-6);
    BOOST_CHECK_SMALL(interiorJacobian[1][1], 1e-6);
    BOOST_CHECK_CLOSE(interiorJacobian[1][2], 1.0/12, 1e-6);
    BOOST_CHECK_CLOSE(interiorJacobian[2][0], -1.0/48, 1e-6);
    BOOST_CHECK_SMALL(interiorJacobian[2][1], 1e-6);
    BOOST_CHECK_CLOSE(interiorJacobian[2][2], 1.0/48, 1e-6);
    
    ctri0.intZDh(totalJacobian, Z);
    ctri0.intDz(totalJacobian);
    
//    printJacobian(totalJacobian);
    
    double eps = 1e-8; // roughly the most accurate point?  eh.
    for (int vv = 0; vv < 3; vv++)
    {
        Vector3d dFdv(0,0,0);
        for (int xyz = 0; xyz < 3; xyz++)
        {
            Triangle3d tri1(tri0);
            tri1[vv][xyz] += eps;
            
            Triangle xyTriangle1(SimpleMesh::Test::oneTri(tri1[0], tri1[1], tri1[2]));
            ClippedTriangle ctri1 = ClippedTriangle(xyTriangle1).partAbove(Z, 0.5);
            
            double vol1 = ctri1.volume(Z);
            
            dFdv[xyz] = (vol1-vol0)/eps;
            
            if (fabs(totalJacobian[vv][xyz]) < 1e-10)
                BOOST_CHECK_SMALL(dFdv[xyz], 1e-6);
            else
                BOOST_CHECK_CLOSE(totalJacobian[vv][xyz], dFdv[xyz], 1e-5);
        }
//        cerr << vv << ": " << dFdv << "\n";
    }
}

static void sIntegrateZDhTriangle_ClipZ(Triangle3d tri0)
{
//    cerr << "SKIPPING\n"; return;
    Triangle xyTriangle0(SimpleMesh::Test::oneTri(tri0[0], tri0[1], tri0[2]));
    
    Rect2d clipBounds(0.3, 0.2, 0.8, 0.7);
    Plane3d zHalf(Vector3d::unit(Z), -0.5);
    
    SensitiveEdge zEdge(xyTriangle0, zHalf);
    
    ClippedTriangle ctri0 = ClippedTriangle(xyTriangle0).partInside(clipBounds);
    ctri0 = ctri0.partAbove(Z, 0.5, &zEdge);
    
    double vol0 = ctri0.volume(Z);
    double area0 = ctri0.area2d(Z);
    
    map<unsigned int, Vector3d> jacobian;
    ctri0.intZDh(jacobian, Z);
    ctri0.intDz(jacobian);
    
    double eps = 1e-8; // roughly the most accurate point?  eh.
    for (int vv = 0; vv < 3; vv++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Triangle3d tri1(tri0);
        tri1[vv][xyz] += eps;
        
        Triangle xyTriangle1(SimpleMesh::Test::oneTri(tri1[0], tri1[1], tri1[2]));
        ClippedTriangle ctri1 = ClippedTriangle(xyTriangle1).partInside(clipBounds);
        ctri1 = ctri1.partAbove(Z, 0.5);
        
        double vol1 = ctri1.volume(Z);
        
        double dFdv = (vol1-vol0)/eps;
        
        if (fabs(jacobian[vv][xyz]) < 1e-10)
            BOOST_CHECK_SMALL(dFdv, 1e-6);
        else
            BOOST_CHECK_CLOSE(jacobian[vv][xyz], dFdv, 1e-4);
    }
}

BOOST_AUTO_TEST_CASE( IntegrateZDhTriangle_ClipZ_Upfacing )
{
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1), xy(1,1,0);
    Triangle3d tri0(o, x+z, y);
    sIntegrateZDhTriangle_ClipZ(tri0);
}

BOOST_AUTO_TEST_CASE( IntegrateZDhTriangle_ClipZ_Downfacing )
{
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1), xy(1,1,0);
    Triangle3d tri0(o, y, x+z);
    sIntegrateZDhTriangle_ClipZ(tri0);
}

void sIntDNNT(Triangle3d tri0)
{
    Triangle xyTriangle0(SimpleMesh::Test::oneTri(tri0[0], tri0[1], tri0[2]));
    
    Rect2d clipBounds(0.1, 0.1, 0.4, 0.4);
    ClippedTriangle ctri0 = ClippedTriangle(xyTriangle0).partInside(clipBounds);
    
    Matrix3d nnT = outerProduct(ctri0.arealNormal(), unit(ctri0.arealNormal()));
    
    map<unsigned int, Matrix3<Vector3d> > jacobian;
    ctri0.intDnnT(jacobian);
    
//    cerr << "Orientation:\n" << nnT << "\n";
//    cerr << "Jacobian:\n";
//    printJacobian(jacobian);
    
    double eps = 1e-8; // roughly the most accurate point?  eh.
    for (int vv = 0; vv < 3; vv++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Triangle3d tri1(tri0);
        tri1[vv][xyz] += eps;
        
        Triangle xyTriangle1(SimpleMesh::Test::oneTri(tri1[0], tri1[1], tri1[2]));
        ClippedTriangle ctri1 = ClippedTriangle(xyTriangle1).partInside(clipBounds);
        
        Matrix3d nnT1 = outerProduct(ctri1.arealNormal(),
            unit(ctri1.arealNormal()));
        
        Matrix3d deltaNNT = (nnT1 - nnT)/eps;
        
        for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
        if (jacobian[vv](ii,jj)[xyz] == 0.0)
            BOOST_CHECK_SMALL(deltaNNT(ii,jj), 1e-6);
        else
            BOOST_CHECK_CLOSE(jacobian[vv](ii,jj)[xyz], deltaNNT(ii,jj), 1e-3);
    }
}

// To test this we can clip the triangle down to have no boundaries left.
BOOST_AUTO_TEST_CASE( IntDNNT_Upfacing )
{
//    cerr << "SKIPPING\n"; return;
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
//    Triangle3d tri0(o+z/2, x+z/2, y+z/2);
    Triangle3d tri0(o + 0.1*z, x - 0.2*z, y + 0.1*x + 0.3*z);
    
    sIntDNNT(tri0);
}

BOOST_AUTO_TEST_CASE( IntDNNT_Downfacing )
{
//    cerr << "SKIPPING\n"; return;
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
//    Triangle3d tri0(o+z/2, x+z/2, y+z/2);
    Triangle3d tri0(o + 0.1*z, y + 0.1*x + 0.3*z, x - 0.2*z);
    
    sIntDNNT(tri0);
}


BOOST_AUTO_TEST_CASE( IntNNTDh_ClipRight_Unflipped )
{
//    cerr << "SKIPPING\n"; return;
    cerr << "\n==== IntNNTDh_CLipRight_Unflipped\n\n";
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    Triangle3d tri0(o+z/2, x+z/2, y+z/2);
    Triangle xyTriangle0(SimpleMesh::Test::oneTri(tri0[0], tri0[1], tri0[2]));
    
    ClippedTriangle ctri0 = ClippedTriangle(xyTriangle0).partBelow(X, 0.5);
//    cerr << "Tri: \n" << ctri0 << "\n";
    
    Matrix3d nnTA = ctri0.nnTA(); //outerProduct(ctri0.arealNormal(), unit(ctri0.arealNormal()));
    
    map<unsigned int, Matrix3<Vector3d> > jacobian;
    ctri0.intNNTDh(jacobian);
    
//    cerr << "Orientation:\n" << nnTA << "\n";
//    cerr << "Jacobian:\n";
//    printJacobian(jacobian);
    
    // The integral contributions along each edge can be worked out in yer head.
    // Here they are:
    
    // Edge 0:
    // dA/dp0 = [0 -3/8 0]
    // dA/dp1 = [0 -1/8 0]
    // dA/dp2 = 0
    
    // Edge 1:
    // dA/dp0 = 0
    // dA/dp1 = [1/8 1/8 0]
    // dA/dp2 = [3/8 3/8 0]
    
    // Edge 2:
    // dA/dp0 = [-1/2 0 0]
    // dA/dp1 = 0
    // dA/dp2 = [-1/2 0 0]
    
    // Total:
    // dA/dp0 = [-1/2 -3/8 0]
    // dA/dp1 = [1/8 0 0]
    // dA/dp2 = [-1/8 3/8 0]
    
    // The sensitivity in question is just the orientation factor times the
    // areal sensitivity.  When the normal vector is unit(z), the answer is:
    CHECK_CLOSE(jacobian[0](Z,Z), Vector3d(-0.5, -0.375, 0));
    CHECK_CLOSE(jacobian[1](Z,Z), Vector3d(0.125, 0, 0));
    CHECK_CLOSE(jacobian[2](Z,Z), Vector3d(-0.125, 0.375, 0));
    
    // make sure all the parts that should be zero ARE zero.
    for (int vert = 0; vert < 3; vert++)
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    if (ii != Z || jj != Z)
    {
        BOOST_CHECK_SMALL(norm(jacobian[vert](ii,jj)), 1e-6);
    }
    
    Vector3d delta(1e-6, 0, 0);
    Triangle3d tri1 = tri0 + delta;
    Triangle xyTriangle1(SimpleMesh::Test::oneTri(tri1[0], tri1[1], tri1[2]));
    ClippedTriangle ctri1 = ClippedTriangle(xyTriangle1).partBelow(X, 0.5);
//    cerr << "Tri perturbed:\n" << ctri1 << "\n";
    
    Matrix3d nnTA1 = ctri1.nnTA(); //outerProduct(ctri1.arealNormal(), unit(ctri1.arealNormal()));
    
    Matrix3d deltaNNT_expected = sDelta(jacobian, delta);
    
//    cerr << "Delta:\n" << nnTA1 - nnTA << "\n";
//    cerr << "Delta expected:\n" << deltaNNT_expected << "\n";
}

BOOST_AUTO_TEST_CASE( IntNNTDh_ClipRight_Tilty )
{
    cerr << "SKIPPING\n"; return;
    cerr << "\n\n==== IntNNTDh_CLipRight_Unflipped\n\n";
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    Triangle3d tri0(o+z/2, x+z/2, y+0.49*z);
    Triangle xyTriangle0(SimpleMesh::Test::oneTri(tri0[0], tri0[1], tri0[2]));
    
    ClippedTriangle ctri0 = ClippedTriangle(xyTriangle0).partBelow(X, 0.5);
//    cerr << "Tri: \n" << ctri0 << "\n";
    
//    Matrix3d nnT = outerProduct(ctri0.arealNormal(), unit(ctri0.arealNormal()))
//        / unit(ctri0.arealNormal())[Z];
    
    Matrix3d nnTdA = ctri0.triangle()->nnTdA();
    
    map<unsigned int, Matrix3<Vector3d> > jacobian;
    ctri0.intNNTDh(jacobian);
    
//    cerr << "Orientation:\n" << nnTdA << "\n";
//    cerr << "Jacobian:\n";
//    printJacobian(jacobian);
    
    // The integral contributions along each edge can be worked out in yer head.
    // Here they are:
    
    // Edge 0:
    // dA/dp0 = [0 -3/8 0]
    // dA/dp1 = [0 -1/8 0]
    // dA/dp2 = 0
    
    // Edge 1:
    // dA/dp0 = 0
    // dA/dp1 = [1/8 1/8 0]
    // dA/dp2 = [3/8 3/8 0]
    
    // Edge 2:
    // dA/dp0 = [-1/2 0 0]
    // dA/dp1 = 0
    // dA/dp2 = [-1/2 0 0]
    
    // Total:
    // dA/dp0 = [-1/2 -3/8 0]
    // dA/dp1 = [1/8 0 0]
    // dA/dp2 = [-1/8 3/8 0]
    
    // The sensitivity in question is just the orientation factor times the
    // areal sensitivity.
    
//    for (int ii = 0; ii < 3; ii++)
//    for (int jj = 0; jj < 3; jj++)
//    {
//        Vector3d jacobian0 = nnTdA(ii,jj) * Vector3d(-0.5, -0.375, 0);
//        
//        Vector3d jacobian1 = nnTdA(ii,jj) * Vector3d(0.125, 0, 0);
//        
//        Vector3d jacobian2 = nnTdA(ii,jj) * Vector3d(-0.125, 0.375, 0.0);
//        
//        cerr << "jacobian0 " << jacobian0 << " vs " << jacobian[0](ii,jj) << "\n";
//        cerr << "jacobian1 " << jacobian1 << " vs " << jacobian[1](ii,jj) << "\n";
//        cerr << "jacobian2 " << jacobian2 << " vs " << jacobian[2](ii,jj) << "\n";
//    }
    
    Vector3d delta(1e-6, 0, 0);
    Triangle3d tri1 = tri0 + delta;
    Triangle xyTriangle1(SimpleMesh::Test::oneTri(tri1[0], tri1[1], tri1[2]));
    ClippedTriangle ctri1 = ClippedTriangle(xyTriangle1).partBelow(X, 0.5);
//    cerr << "Tri perturbed:\n" << ctri1 << "\n";
    
    Matrix3d nnT1dA = ctri1.triangle()->nnTdA();
    
    Matrix3d deltaNNT = ctri1.triangle()->nnTdA()*ctri1.area2d(Z) -
        ctri0.triangle()->nnTdA()*ctri0.area2d(Z);
    Matrix3d deltaNNT_expected = sDelta(jacobian, delta);
    
//    cerr << "Delta:\n" << deltaNNT << "\n";
//    cerr << "Delta expected:\n" << deltaNNT_expected << "\n";
}

static ClippedTriangle sBugsy(vector<ControlVertex> cvs)
{
    Triangle3d ctri0(cvs[0].point(), cvs[1].point(), cvs[2].point());
    SimpleMesh::Triangle* meshTri0 = new SimpleMesh::Triangle(ctri0);
    meshTri0->controlVertices(cvs[0], cvs[1], cvs[2]);
    meshTri0->edgeControlVertices(cvs[0], cvs[1], cvs[2]);
    meshTri0->cacheSensitivity();
    
    ClippedTriangle clipt0(*meshTri0);
//    cerr << "Pre-clipped:\n" << clipt0 << "\n";
    
    // Clipping:
    // y in [-0.5 0.5]
    // z in [-0.5 0.5]
    // x in [0 1]
    
    
//    Rect3d pb(-0.5, 0, -0.5, 0.5, 1, 0.5);
    Rect3d pb(1, -2.0/3, -1, 3, 2.0/3, -1.0/3);
    
    Plane3d bottomPlane(Vector3d::unit(0), -pb.p1[0]);
    Plane3d topPlane(Vector3d::unit(0), -pb.p2[0]);
    
    Plane3d topY(Vector3d::unit(Y), -pb.p2[1]),
        botY(Vector3d::unit(Y), -pb.p1[1]),
        topZ(Vector3d::unit(Z), -pb.p2[2]),
        botZ(Vector3d::unit(Z), -pb.p1[2]);
    
//    cerr << "Planes:\n"
//        << botY << "\n" << topY << "\n" << botZ << "\n" << topZ << "\n";
    
    SensitiveEdge* edgeTopY = new SensitiveEdge(*clipt0.triangle(), topY);
    SensitiveEdge* edgeBotY = new SensitiveEdge(*clipt0.triangle(), botY);
    
    clipt0 = clipt0.partAbove(Y, pb.p1[1], edgeBotY);
    clipt0 = clipt0.partBelow(Y, pb.p2[1], edgeTopY);
    
    
    SensitiveEdge* edgeAbove = new SensitiveEdge(*clipt0.triangle(), topPlane);
    SensitiveEdge* edgeBelow = new SensitiveEdge(*clipt0.triangle(), bottomPlane);
    
//    cerr << "Above: " << topPlane << " edge " << edgeAbove->line()
//        << "\n";
//    cerr << "Below: " << bottomPlane << " edge " << edgeBelow->line() << "\n";
    
    clipt0 = clipt0.partAbove(X, pb.p1[0], edgeBelow);
    clipt0 = clipt0.partBelow(X, pb.p2[0], edgeAbove);
    
    return clipt0;
}


// Stick a triangle up so its top vertex lies exactly on a clipping boundary.
// Make sure that the Jacobians are ok.
//BOOST_AUTO_TEST_CASE( VertexOnBoundary )
//{
//    // Triangle faces down.  This is a regression test and I'm imitating
//    // the problem there.
//    
//    Vector3d cv0(0,0,0);
//    Vector3d cv1(0,1,0);
//    Vector3d cv2(1,0,0);
//    
//    vector<ControlVertex> cvs0, cvs1;
//    cvs0.push_back(ControlVertex(cv0, 0));
//    cvs0.push_back(ControlVertex(cv1, 1));
//    cvs0.push_back(ControlVertex(cv2, 4));
//    
//    SimpleMesh::Triangle tri( Triangle3d(cv0, cv1, cv2) );
//    tri.controlVertices(cvs0);
//    
//    // Now clip the triangle to the y = 1 plane.
//    
//    ClippedTriangle preClipped(tri);
//    ClippedTriangle clipped = preClipped.partBelow(1, 1.0);
//}

BOOST_AUTO_TEST_CASE( BugChaser )
{
    cerr << "SKIPPING\n"; return;
    cerr << "\n\n===== BugChaser\n\n";
    
    // First triangle
//    Vector3d cv0(-0.1830127018922193, 0.81698729810778059, 0);
//    Vector3d cv1(0.68301270189221941, 0.3169872981077807, 0);
//    Vector3d cv4(-0.1830127018922193, 0.81698729810778059, 1);

    // Second triangle
//    Vector3d cv0(0.683013, 0.316987, 0);
//    Vector3d cv1(0.683013, 0.316987, 1);
//    Vector3d cv4(-0.183013, 0.816987, 1);
    
    // Third triangle
//    Vector3d cv0(1.18301, 1.18301, 0);
//    Vector3d cv1(1.18301, 1.18301, 1);
//    Vector3d cv4(0.683013, 0.316987, 1);
    
    // Fourth triangle
//    Vector3d cv0(0.683013, 0.316987, 0);
//    Vector3d cv1(1.18301, 1.18301, 0);
//    Vector3d cv4(0.683013, 0.316987, 1);

    // From standard block, Triangle 1
    Vector3d cv0(-2, -1.33333, -0.666667);
    Vector3d cv1(-2, 1.33333, -0.666667);
    Vector3d cv4(2, -1.33333, -0.666667);
    
//    Vector3d delta(1e-5, 0, 0);
    Vector3d delta(0, 1e-5, 0);
    
    vector<ControlVertex> cvs0, cvs1;
    cvs0.push_back(ControlVertex(cv0, 0));
    cvs0.push_back(ControlVertex(cv1, 1));
    cvs0.push_back(ControlVertex(cv4, 4));
    
    ClippedTriangle clipt0 = sBugsy(cvs0);
    cerr << "Clipped triangle:\n"
        << "quickPatch(" << clipt0 << ", 'g');\n";
    
    cvs1 = cvs0;
    for (int nn = 0; nn < 3; nn++)
    {
        cvs1[nn].point(cvs1[nn].point() + delta);
    }
    ClippedTriangle clipt1 = sBugsy(cvs1);
    cerr << "New one:\n"
        << "quickPatch(" << clipt1 << ", 'g');\n";
    
    map<unsigned int, Matrix3<Vector3d> > DNNT, NNTDh, ADnnT;
    // = clipt0.DNNT(X);
    DNNT = clipt0.DNNT(clipt0.triangle()->upAxis());
    clipt0.intNNTDh(NNTDh);
    clipt0.intDnnT(ADnnT);
    
    cerr << "DNNT = \n";
    printJacobian(DNNT);
    cerr << "\nNNTDh =\n";
    printJacobian(NNTDh);
    cerr << "\nADnnT =\n";
    printJacobian(ADnnT);
    
    Matrix3d nnTA0 = clipt0.nnTA();
    Matrix3d nnTA1 = clipt1.nnTA();
    
    Matrix3d deltaNNT = clipt1.nnTA() - clipt0.nnTA();
    Matrix3d deltaNNT_expected = sDelta(DNNT, delta);
    
    cerr << "deltaNNT =\n" << deltaNNT << "\n";
    cerr << "expected\n" << deltaNNT_expected << "\n";
    cerr << "\tareal\n" << sDelta(ADnnT, delta) << "\n";
    cerr << "\tbdy\n" << sDelta(NNTDh, delta) << "\n";
    
    Matrix3d DNNT_meas = deltaNNT/norm(delta);
    Matrix3d DNNT_exp = deltaNNT_expected/norm(delta);
    
    cerr << "DNNT_meas =\n" << DNNT_meas << "\n";
    cerr << "DNNT_exp =\n" << DNNT_exp << "\n";
    
    if (vectorNorm(deltaNNT_expected) > 1e-5)
    {
        CHECK_MATRIX_CLOSE(DNNT_meas, DNNT_exp);
    }
    
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE( SimpleBugChaser )
{
    cerr << "SKIPPING\n"; return;
//    cerr << "\n\n===== SimpleBugChaser\n\n";
    // Normal along Y.  Fails.
    Vector3d cv0(0,0.8,0);
    Vector3d cv1(1.0,0,0);
    Vector3d cv4(0,0.8,1);

    // Normal along X. Succeeds.
//    Vector3d cv0(0,1.0,0);
//    Vector3d cv1(0.8,0,0);
//    Vector3d cv4(0,1.0,1);
    
    Vector3d delta(1e-2, 0, 0);
    
    vector<ControlVertex> cvs0, cvs1;
    cvs0.push_back(ControlVertex(cv0, 0));
    cvs0.push_back(ControlVertex(cv1, 1));
    cvs0.push_back(ControlVertex(cv4, 4));
    
    ClippedTriangle clipt0 = sBugsy(cvs0);    
//    cerr << "Clipped triangle:\n"
//        << "quickPatch(" << clipt0 << ", 'g');\n";
//    cerr << "NNTA =\n" << clipt0.nnTA() << "\n";
    
    cvs1 = cvs0;
    for (int nn = 0; nn < 3; nn++)
    {
        cvs1[nn].point(cvs1[nn].point() + delta);
    }
    ClippedTriangle clipt1 = sBugsy(cvs1);
//    cerr << "New one:\n"
//        << "quickPatch(" << clipt1 << ", 'g');\n";
//    cerr << "NNTA =\n" << clipt1.nnTA() << "\n";
    
    map<unsigned int, Matrix3<Vector3d> > DNNT, NNTDh, ADnnT;
    // = clipt0.DNNT(X);
    DNNT = clipt0.DNNT(clipt0.triangle()->upAxis());
    clipt0.intNNTDh(NNTDh);
    clipt0.intDnnT(ADnnT);
    
//    cerr << "DNNT = \n";
//    printJacobian(DNNT);
//    cerr << "\nNNTDh =\n";
//    printJacobian(NNTDh);
//    cerr << "\nADnnT =\n";
//    printJacobian(ADnnT);
    
    Matrix3d nnTA0 = clipt0.nnTA();
    Matrix3d nnTA1 = clipt1.nnTA();
    
    Matrix3d deltaNNT = clipt1.nnTA() - clipt0.nnTA();
    Matrix3d deltaNNT_expected = sDelta(DNNT, delta);
    
//    cerr << "deltaNNT =\n" << deltaNNT << "\n";
//    cerr << "expected\n" << deltaNNT_expected << "\n";
//    cerr << "\tareal\n" << sDelta(ADnnT, delta) << "\n";
//    cerr << "\tbdy\n" << sDelta(NNTDh, delta) << "\n";
    
    BOOST_CHECK(true);
}


/*
BOOST_AUTO_TEST_CASE( IntNNTDh_ClipRight_Flipped )
{
    cerr << "\n==== IntNNTDh_ClipRight_Flipped\n\n";
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    Triangle3d tri0(o+z/2, x+z/2, y+z/2);
    Triangle xyTriangle0(sOneFlippedTri(tri0[0], tri0[1], tri0[2]));
    
    ClippedTriangle ctri0 = ClippedTriangle(xyTriangle0).partBelow(X, 0.5);
    cerr << "Tri: \n" << ctri0 << "\n";
    
    Matrix3d nnT = outerProduct(ctri0.arealNormal(), unit(ctri0.arealNormal()));
    
    map<unsigned int, Matrix3<Vector3d> > jacobian;
    ctri0.intNNTDh(jacobian);
    
    cerr << "Orientation:\n" << nnT << "\n";
    cerr << "Jacobian:\n";
    printJacobian(jacobian);

    // The integral contributions are the same as in the unflipped case, but
    // should have p1 and p2 reversed.
    
    // Edge 0:
    // dA/dp0 = [0 -3/8 0]
    // dA/dp1 = 0
    // dA/dp2 = [0 -1/8 0]
    
    // Edge 1:
    // dA/dp0 = 0
    // dA/dp1 = [3/8 3/8 0]
    // dA/dp2 = [1/8 1/8 0]
    
    // Edge 2:
    // dA/dp0 = [-1/2 0 0]
    // dA/dp1 = [-1/2 0 0]
    // dA/dp2 = 0
    
    // Total:
    // dA/dp0 = [-1/2 -3/8 0]
    // dA/dp1 = [-1/8 3/8 0]
    // dA/dp2 = [1/8 0 0]
    
    // The sensitivity in question is just the orientation factor times the
    // areal sensitivity.  When the normal vector is unit(z), the answer is:
    CHECK_CLOSE(jacobian[0](Z,Z), Vector3d(-0.5, -0.375, 0));
    CHECK_CLOSE(jacobian[1](Z,Z), Vector3d(-0.125, 0.375, 0));
    CHECK_CLOSE(jacobian[2](Z,Z), Vector3d(0.125, 0, 0));
    
    // make sure all the parts that should be zero ARE zero.
    for (int vert = 0; vert < 3; vert++)
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    if (ii != Z || jj != Z)
    {
        BOOST_CHECK_SMALL(norm(jacobian[vert](ii,jj)), 1e-6);
    }
}
*/

static Matrix3d sDelta(const map<unsigned int, Matrix3<Vector3d> > & jacobian,
    Vector3d delta)
{
    Matrix3d out;
    
    map<unsigned int, Matrix3<Vector3d> >::const_iterator itr;
    for (itr = jacobian.begin(); itr != jacobian.end(); itr++)
    {
        for (int mm = 0; mm < 9; mm++)
            out[mm] += dot(itr->second[mm], delta);
    }
    return out;
}




