/*
 *  testOrientedVoxels.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/29/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestVoxelizer

#include <boost/test/unit_test.hpp>

#include "OrientedVoxels.h"
#include "VoxelizePolyhedra.h"
#include "SimpleMesh/Search.h"
#include "UserPreferences.h"

#include <iomanip>
#include <cstdlib>

using namespace std;
using namespace RLE;
using namespace YeeUtilities;


SupportRegion3 support(int x0, int y0, int z0, int x1, int y1, int z1)
{
    SupportRegion3 supp;
    supp.mark(x0, y0, z0, x1, y1, z1);
    return supp;
}

SupportRegion3 support(const Rect3i & r)
{
    SupportRegion3 supp;
    supp.mark(r.p1[0], r.p1[1], r.p1[2], r.p2[0], r.p2[1], r.p2[2]);
    return supp;
}

//vector<TrackedPolyhedron*> getPolyhedra(
//    const map<string, TrackedPolyhedron*> & polyhedra)
//{
//    vector<TrackedPolyhedron*> out;
//    map<string, TrackedPolyhedron*>::const_iterator itr;
//    for (itr = polyhedra.begin(); itr != polyhedra.end(); itr++)
//    {
//        out.push_back(itr->second);
//    }
//    return out;
//}
//
//vector<vector<SimpleMesh::Triangle> > getPolyhedra(
//    const map<string, TrackedPolyhedron*> & polyhedra)
//{
//    vector<vector<SimpleMesh::Triangle> > out(polyhedra.size());
//    int idx = 0;
//
//    map<string, TrackedPolyhedron*>::const_iterator itr;
//    for (itr = polyhedra.begin(); itr != polyhedra.end(); itr++)
//    {
//        out[idx] = itr->second->triangles();
//    }
//    SimpleMesh::initMeshSensitivity(out);
//    return out;
//}
//
//vector<vector<SimpleMesh::Triangle> > getPolyhedra(
//    const std::vector<TrackedPolyhedron*> & polyhedra)
//{
//    vector<vector<SimpleMesh::Triangle> > out(polyhedra.size());
//
//    for (int idx = 0; idx < polyhedra.size(); idx++)
//    {
//        out[idx] = polyhedra[idx]->triangles();
//    }
//    SimpleMesh::initMeshSensitivity(out);
//    return out;
//}

BOOST_AUTO_TEST_CASE( Basics )
{
//    cerr << "SKIPPING.\n"; return;
    Rect3d bounds(0,0,0,10,10,10);
    Rect3d clipBounds(1,1,1,9,9,9);
    Rect3i yeeCells(0,0,0,9,9,9);
    
    OrientedVoxels ov1(yeeCells, bounds);
    BOOST_CHECK(ov1.voxelBounds() == yeeCells);
    BOOST_CHECK(ov1.realBounds() == bounds);
    BOOST_CHECK(ov1.clipBounds() == bounds);
    
    OrientedVoxels ov2(yeeCells, bounds, clipBounds);
    BOOST_CHECK(ov2.voxelBounds() == yeeCells);
    BOOST_CHECK(ov2.realBounds() == bounds);
    BOOST_CHECK(ov2.clipBounds() == clipBounds);
    
    OrientedVoxels ov3;
    BOOST_CHECK(ov3.voxelBounds() == Rect3i(0,0,0,0,0,0));
    BOOST_CHECK(ov3.realBounds() == Rect3d(0,0,0,0,0,0));
    BOOST_CHECK(ov3.clipBounds() == Rect3d(0,0,0,0,0,0));
    ov3.realBounds(bounds);
    ov3.clipBounds(clipBounds);
    ov3.voxelBounds(yeeCells);
    BOOST_CHECK(ov3.voxelBounds() == yeeCells);
    BOOST_CHECK(ov3.realBounds() == bounds);
    BOOST_CHECK(ov3.clipBounds() == clipBounds);
}

// Test downsampling the fill factors in a single Yee cell.
// Test each octant.
BOOST_AUTO_TEST_CASE( DownsampleFillFactor )
{
//    cerr << "SKIPPING.\n"; return;
    Rect3d bounds(0, 0, 0, 1, 1, 1);
    Rect3i yeeCells(0, 0, 0, 0, 0, 0);
    
    vector<DynamicRLE3<double> > fillFactors(1);
    fillFactors[0].mark(0,0,0,1,1,0, 1.0); // z = 0: four cells with value 1
    fillFactors[0].mark(0,0,1,1,1,1, 0.0); // z = 1: four cells with value 0
    fillFactors[0].mark(1,1,0, 0.0); // make this corner empty too.
    
    // Remember that out-of-bounds voxels clip in this downsampling!
    for (int oct = 0; oct < 3; oct++)
    {
        OrientedVoxels vox(yeeCells, bounds);
        vox.downsampleFillFactors(fillFactors, oct);
        BOOST_CHECK_EQUAL(vox.fillFactors(0).at(0,0,0), 1.0);
    }
    
    {
        OrientedVoxels vox(yeeCells, bounds);
        vox.downsampleFillFactors(fillFactors, 3);
        BOOST_CHECK_EQUAL(vox.fillFactors(0).at(0,0,0), 0.75);
    }
    
    for (int oct = 4; oct < 7; oct++)
    {
        OrientedVoxels vox(yeeCells, bounds);
        vox.downsampleFillFactors(fillFactors, oct);
        BOOST_CHECK_EQUAL(vox.fillFactors(0).at(0,0,0), 0.5);
    }
    
    {
        OrientedVoxels vox(yeeCells, bounds);
        vox.downsampleFillFactors(fillFactors, 7);
        BOOST_CHECK_EQUAL(vox.fillFactors(0).at(0,0,0), 0.375);
    }
}

// Verify that faces around the outer boundary of the voxelization region do not contribute
// to the orientation of any cells in the interior.
BOOST_AUTO_TEST_CASE( NoOrientationAtBoundary )
{
//    cerr << "SKIPPING.\n"; return;
    Rect3d bounds(0, 0, 0, 3, 3, 3);
    Rect3i yeeCells(0, 0, 0, 2, 2, 2);
    
    std::vector<SimpleMesh::ControlVertex> controlVertices;
    std::vector<SimpleMesh::Triangle> triangles;
    SimpleMesh::Test::oneBoxMesh(controlVertices, triangles, bounds);
    vector<vector<SimpleMesh::Triangle>> meshes = {triangles};
    
    OrientedVoxels vox(yeeCells, bounds);
    
    vox.initOrientations(meshes);
    
    BOOST_CHECK_EQUAL(0, vox.orientations().numRuns());
}

BOOST_AUTO_TEST_CASE( Orientation )
{
//    cerr << "SKIPPING.\n"; return;
    BOOST_CHECK(true);
    
    Rect3d bounds(0,0,0,1,1,1);
    Rect3i yeeCells(0,0,0,0,0,0);
    
    Rect3d bigBlock(0,0,0,1,1,1);
    Rect3d smallBlock(0.49,0,0,1,0.51,1);
    
    std::vector<SimpleMesh::ControlVertex> controlVertices;
    std::vector<SimpleMesh::Triangle> triangles;
    SimpleMesh::Test::oneBoxMesh(controlVertices, triangles, smallBlock);
    vector<vector<SimpleMesh::Triangle>> meshes = {triangles};
    
//
//    Assembly assembly;
////    TestingGoodies::twoBoxAssembly(assembly, bigBlock, smallBlock);
//    TestingGoodies::oneBoxAssembly(assembly, smallBlock);
//
//    vector<TrackedPolyhedron*> polyhedra = TestingGoodies::listPolyhedra(
//        assembly.polyhedra());
//    vector<vector<SimpleMesh::Triangle> > meshes = getPolyhedra(polyhedra);
    
//    int material = 0;
    
    vector<Matrix3d> expectedOrientations(8); // per octant
    
    expectedOrientations[0] = Matrix3d(
        1, 0, 0,
        0, 0, 0,
        0, 0, 0);
    expectedOrientations[1] = Matrix3d(
        1.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0);
    expectedOrientations[2] = Matrix3d(
        .51/.52, 0.0, 0.0,
        0.0, .01/.52, 0.0,
        0.0, 0.0, 0.0);
    expectedOrientations[3] = Matrix3d(
        0.5, 0.0, 0.0,
        0.0, 0.5, 0.0,
        0.0, 0.0, 0.0);
    for (int oo = 4; oo < 8; oo++)
        expectedOrientations[oo] = expectedOrientations[oo-4];
    
    for (int octant = 0; octant < 8; octant++)
    {
//        cerr << "Orientation(" << octant << "):\n";
        OrientedVoxels voxels(yeeCells, realBounds(yeeCells, octant), bounds);
        voxels.initOrientations(meshes);
        
//        cerr << "\torientation(" << octant << ") = "
//            << voxels.orientations().at(0,0,0) << "\n";
//        cerr << "\texpected " << expectedOrientations[octant] << "\n";
        for (int nn = 0; nn < 9; nn++)
        {
            BOOST_CHECK_CLOSE(voxels.orientations().at(0,0,0)[nn],
                expectedOrientations[octant][nn], 1e-9);
        }
    }
}

BOOST_AUTO_TEST_CASE( OrientationSensitivity_OneVoxel )
{
//    cerr << "SKIPPING.\n"; return;
    BOOST_CHECK(true);
    
    Rect3d bounds(0,0,0,1,1,1);
    Rect3i yeeCells(0,0,0,0,0,0);
    
    Rect3d bigBlock(0,0,0,1,1,1);
    Rect3d smallBlock(0.25,0,0,1,0.75,1);
    
    // The moving vertices of the small block have ID numbers
    //  8, 10, 11, 12, 14, 15
//    int movingVertices[] = { 8, 10, 11, 12, 14, 15 };
    std::array<int, 8> movingVertices = {0, 1, 2, 3, 4, 5, 6, 7};
    
    std::vector<SimpleMesh::ControlVertex> controlVertices, controlVertices2;
    std::vector<SimpleMesh::Triangle> triangles, triangles2;
    SimpleMesh::Test::oneBoxMesh(controlVertices, triangles, smallBlock);
    vector<vector<SimpleMesh::Triangle>> meshes = {triangles};
    
//    Assembly assembly;
//    TestingGoodies::oneBoxAssembly(assembly, smallBlock);
//    vector<TrackedPolyhedron*> polyhedra = TestingGoodies::listPolyhedra(
//        assembly.polyhedra());
//    vector<vector<SimpleMesh::Triangle> > meshes = getPolyhedra(polyhedra);
    
//    int material = 0;
    Vector2d delta(1.1e-8, 0.9e-8);
    
    for (int octant = 0; octant < 8; octant++)
    {
//        cout << "Octant " << octant << ":\n";
        
        OrientedVoxels voxels(yeeCells, realBounds(yeeCells, octant), bounds);
        voxels.initOrientations(meshes);
        
        Rect3d smallBlock2(smallBlock);
        smallBlock2.p1[0] += delta[0];
        smallBlock2.p2[1] += delta[1];
        
//        Assembly assembly2;
////        TestingGoodies::twoBoxAssembly(assembly2, bigBlock, smallBlock2);
//        TestingGoodies::oneBoxAssembly(assembly2, smallBlock2);
//        vector<TrackedPolyhedron*> polyhedra2 = TestingGoodies::listPolyhedra(
//            assembly2.polyhedra());
//        vector<vector<SimpleMesh::Triangle> > meshes2 = getPolyhedra(polyhedra2);
        
        SimpleMesh::Test::oneBoxMesh(controlVertices2, triangles2, smallBlock2);
        vector<vector<SimpleMesh::Triangle>> meshes2 = {triangles2};
        
        OrientedVoxels vox2(yeeCells, realBounds(yeeCells, octant), bounds);
        vox2.initOrientations(meshes2);

        Matrix3d deltaOrientation = vox2.orientations().at(0,0,0)
            - voxels.orientations().at(0,0,0);
        Matrix3d deltaOrientationExpected;
        
        vector<Vector3d> deltaVertices = SimpleMesh::Test::subtractPositions(controlVertices2, controlVertices);
        
        for (int vv = 0; vv < movingVertices.size(); vv++)
        for (int xyz = 0; xyz < 3; xyz++)
        {
            Matrix3<double> dOrientation;
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                dOrientation(i,j) = voxels.orientationJacobians(i, j)
                    .find(movingVertices[vv])->second.at(0,0,0)[xyz];
            }
            
            deltaOrientationExpected = deltaOrientationExpected +
                dOrientation * deltaVertices.at(movingVertices[vv])[xyz];
        }
        
        Matrix3d derivMeas = deltaOrientation / norm(delta);
        Matrix3d derivExp = deltaOrientationExpected / norm(delta);
        
        double derivErrNorm = vectorNorm(derivExp - derivMeas);
        double derivRelErr = derivErrNorm / vectorNorm(derivMeas);
        
        if (fabs(vectorNorm(derivMeas)) > 1e-9)
        {
            BOOST_CHECK_SMALL(derivRelErr, 1e-6);
        }
        else
        {
            BOOST_CHECK_SMALL(vectorNorm(derivExp), 1e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE( OrientationSensitivity_OneVoxel_15Degrees )
{
//    cerr << "SKIPPING.\n"; return;
    BOOST_CHECK(true);
    
    UserPreferences::set("noNormalizedOrientation");
    
    Rect3d bounds(0,0,0,1,1,1);
    Rect3i yeeCells(0,0,0,0,0,0);
    
    Rect3d bigBlock(0,0,0,1,1,1);
    Rect3d smallBlock(-0.5, -0.5, 0.0, 0.5, 0.5, 1.0);
    Matrix3d rot = Matrix3d::rotation(-M_PI/6, Vector3d(0,0,1));
    Vector3d translate(0.5, 1.0, 0.0);
//
//    Assembly assembly;
////    TestingGoodies::twoBoxAssembly(assembly, bigBlock, smallBlock, rot, translate);
//    TestingGoodies::oneBoxAssembly(assembly, smallBlock, rot, translate);
//
//    vector<TrackedPolyhedron*> polyhedra = TestingGoodies::listPolyhedra(assembly.polyhedra());
//    vector<vector<SimpleMesh::Triangle> > meshes = getPolyhedra(polyhedra);


    std::vector<SimpleMesh::ControlVertex> controlVertices, controlVertices2;
    std::vector<SimpleMesh::Triangle> triangles, triangles2;
    SimpleMesh::Test::oneBoxMesh(controlVertices, triangles, smallBlock, rot, translate);
    vector<vector<SimpleMesh::Triangle>> meshes = {triangles};


    Vector3d delta(1e-5, 0, 0);
    
    BOOST_CHECK(true);
    for (int octant = 0; octant < 8; octant++)
    {
//        cout << "Octant " << octant << ":\n";
        
//        cerr << "\n==== Original\n";
        OrientedVoxels vox(yeeCells, realBounds(yeeCells, octant), bounds);
        vox.initOrientations(meshes);
        
//        Assembly assembly2;
////        TestingGoodies::twoBoxAssembly(assembly2, bigBlock, smallBlock, rot, translate+delta);
//        TestingGoodies::oneBoxAssembly(assembly2, smallBlock, rot, translate+delta);
//
//        vector<TrackedPolyhedron*> polyhedra2 = TestingGoodies::listPolyhedra(
//            assembly2.polyhedra());
//        vector<vector<SimpleMesh::Triangle> > meshes2 = getPolyhedra(polyhedra2);
    
        SimpleMesh::Test::oneBoxMesh(controlVertices2, triangles2, smallBlock, rot, translate+delta);
        vector<vector<SimpleMesh::Triangle>> meshes2 = {triangles2};
    
//        cerr << "\n==== Perturbed\n";
        OrientedVoxels vox2(yeeCells, realBounds(yeeCells, octant), bounds);
        vox2.initOrientations(meshes2);

        Matrix3d deltaOrientation = vox2.orientations().at(0,0,0)
            - vox.orientations().at(0,0,0);
        Matrix3d deltaOrientationExpected;
        
        vector<Vector3d> deltaVertices = SimpleMesh::Test::subtractPositions(controlVertices2, controlVertices);
        BOOST_CHECK(true);
//        cerr << "Orientation:\n" << vox.orientations().at(0,0,0) << "\n"
//            << vox2.orientations().at(0,0,0) << "\n";
        
//        for (int vv = 0; vv < 16; vv++)
//            cerr << "cv" << vv << " = " << assembly.controlVertices()[vv] << "\n";
//        
//        for (int vv = 0; vv < deltaVertices.size(); vv++)
//        {
//            cerr << "cv" << vv << ": moved " << deltaVertices[vv] << "\n";
//        }
        
        for (int vv = 0; vv < deltaVertices.size(); vv++)
        for (int xyz = 0; xyz < 3; xyz++)
        {
            Matrix3<double> dOrientation;
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
//                std::cerr << i << ", " << j << "\n";
                BOOST_CHECK(true);
                map<unsigned int, DynamicRLE3<Vector3d > > oriJac = vox.orientationJacobians(i,j);
                map<unsigned int, DynamicRLE3<Vector3d> >::const_iterator itr;
                itr = oriJac.find(vv);
                if (itr == oriJac.end())
                {
//                    std::cerr << "No entry for " << vv << ".\n";
                }
                else
                {
//                    std::cerr << "Has " << vv << ".\n";
                    dOrientation(i,j) = itr->second.at(0,0,0)[xyz];
                }
                BOOST_CHECK(true);
            }
            
            BOOST_CHECK(true);
//            if (vectorNorm(dOrientation) > 1e-5)
//            {
//                cerr << "d/dcv" << vv << "." << char('x'+xyz) << " =\n"
//                    << dOrientation << "\n";
//            }
            
            deltaOrientationExpected = deltaOrientationExpected +
                dOrientation * deltaVertices.at(vv)[xyz];
        }
        
        BOOST_CHECK(true);
        Matrix3d derivMeas = deltaOrientation / norm(delta);
        Matrix3d derivExp = deltaOrientationExpected / norm(delta);
        
//        cerr << "Delta " << deltaOrientation << "\n";
//        cerr << "Measured " << derivMeas << "\n";
//        cerr << "Expected " << derivExp << "\n";
        
        double derivErrNorm = vectorNorm(derivExp - derivMeas);
        double derivRelErr = derivErrNorm / vectorNorm(derivMeas);
        
        if (fabs(vectorNorm(derivMeas)) > 1e-5)
        {
            BOOST_CHECK_SMALL(derivRelErr, 1e-6);
        }
        else
            BOOST_CHECK_SMALL(vectorNorm(derivExp), 1e-6);
    }
}

static void sHorribleTest(int chosenTriangle);

BOOST_AUTO_TEST_CASE( HorribleTest )
{
//    for (int chosenTriangle = 0; chosenTriangle < 12; chosenTriangle++)
//    {
//        sHorribleTest(chosenTriangle);
//    }
    sHorribleTest(-1);
}

static void sHorribleTest(int chosenTriangle)
{
//    cerr << "Chose triangle " << chosenTriangle << "\n";
//    int chosenDirection = 0;
//    system("mkfifo toMatlab");
//    UserPreferences::set("noNormalizedOrientation");
    Rect3d bounds(0,0,0,1,1,1);
    Rect3i yeeCells(0,0,0,0,0,0);
    
    Rect3d bigBlock(0,0,0,1,1,1);
    Rect3d smallBlock(-0.25, -0.25, 0.25, 0.25, 0.25, 0.75);
    Matrix3d rot = Matrix3d::rotation(-M_PI/6, Vector3d(0,0,1));
    Vector3d translate(0.75, 0.5, 0.0);
    
//    Assembly assembly;
////    TestingGoodies::twoBoxAssembly(assembly, bigBlock, smallBlock, rot, translate);
//    TestingGoodies::oneBoxAssembly(assembly, smallBlock, rot, translate);
//
//    vector<TrackedPolyhedron*> polyhedra = TestingGoodies::listPolyhedra(
//        assembly.polyhedra());
//    polyhedra.pop_back();
//    vector<vector<SimpleMesh::Triangle> > meshes = getPolyhedra(polyhedra);
//
    
    std::vector<SimpleMesh::ControlVertex> controlVertices, controlVertices2;
    std::vector<SimpleMesh::Triangle> triangles, triangles2;
    SimpleMesh::Test::oneBoxMesh(controlVertices, triangles, smallBlock, rot, translate);
    vector<vector<SimpleMesh::Triangle>> meshes = {triangles};
    
//    polyhedra[0]->dbgSaveOrientedTriangles(chosenDirection);
//    polyhedra[0]->dbgSaveUpwardTriangles(-Vector3d::unit(chosenDirection));
    
    
//    ofstream matlab("pipe", ios_base::app);
//    matlab << "figure(1)\n";
//    matlab << "plot(1:10, 'ro');" << endl;
    Vector3d delta(1e-8, 0, 0);
    
//    for (int octant = 1; octant < 2; octant++)
    for (int octant = 0; octant < 8; octant++)
//    for (int octant = 7; octant < 8; octant++)
    {
//        cout << "Octant " << octant << " tri " << chosenTriangle << ":\n";
        
//        cerr << "\n==== Original\n";
        OrientedVoxels vox(yeeCells, realBounds(yeeCells, octant), bounds);
        vox.initOrientations(meshes);
        
//        Assembly assembly2;
////        TestingGoodies::twoBoxAssembly(assembly2, bigBlock, smallBlock, rot, translate+delta);
//        TestingGoodies::oneBoxAssembly(assembly2, smallBlock, rot, translate+delta);
//
//        vector<TrackedPolyhedron*> polyhedra2 = TestingGoodies::listPolyhedra(
//            assembly2.polyhedra());
//
////        polyhedra2[0] = polyhedra2[1];
//        polyhedra2.pop_back();
//        vector<vector<SimpleMesh::Triangle> > meshes2 = getPolyhedra(polyhedra2);
        
        SimpleMesh::Test::oneBoxMesh(controlVertices2, triangles2, smallBlock, rot, translate+delta);
        vector<vector<SimpleMesh::Triangle>> meshes2 = {triangles2};
        
//        cerr << "\n==== Perturbed\n";
        OrientedVoxels vox2(yeeCells, realBounds(yeeCells, octant), bounds);
        vox2.initOrientations(meshes2);

        Matrix3d deltaOrientation = vox2.orientations().at(0,0,0)
            - vox.orientations().at(0,0,0);
        Matrix3d deltaOrientationExpected;
        
        vector<Vector3d> deltaVertices = SimpleMesh::Test::subtractPositions(controlVertices2, controlVertices);
        
//        cerr << "Orientation:\n" << vox.orientations().at(0,0,0) << "\n"
//            << vox2.orientations().at(0,0,0) << "\n";
        
//        for (int vv = 0; vv < 16; vv++)
//            cerr << "cv" << vv << " = " << assembly.controlVertices()[vv] << "\n";
//        
//        for (int vv = 0; vv < 16; vv++)
//            cerr << "cv" << vv << ": moved " << deltaVertices[vv] << "\n";
        
        for (int vv = 0; vv < deltaVertices.size(); vv++)
        for (int xyz = 0; xyz < 3; xyz++)
        {
            Matrix3<double> dOrientation;
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                dOrientation(i,j) = vox.orientationJacobians(i, j)
                    .find(vv)->second.at(0,0,0)[xyz];
            }
            
//            if (vectorNorm(dOrientation) > 1e-5)
//            {
//                cerr << "d/dcv" << vv << "." << char('x'+xyz) << " =\n"
//                    << dOrientation << "\n";
//            }
            
            deltaOrientationExpected = deltaOrientationExpected +
                dOrientation * deltaVertices.at(vv)[xyz];
        }
        
        Matrix3d derivMeas = deltaOrientation / norm(delta);
        Matrix3d derivExp = deltaOrientationExpected / norm(delta);
        
//        cerr << "Delta " << deltaOrientation << "\n";
//        cerr << "Delta exp " << deltaOrientationExpected << "\n";
//        cerr << "Measured " << derivMeas << "\n";
//        cerr << "Expected " << derivExp << "\n";
        
        double derivErrNorm = vectorNorm(derivExp - derivMeas);
        double derivRelErr = derivErrNorm / vectorNorm(derivMeas);
        
        if (fabs(vectorNorm(derivMeas)) > 1e-5)
        {
            BOOST_CHECK_SMALL(derivRelErr, 1e-6);
        }
        else
            BOOST_CHECK_SMALL(vectorNorm(derivExp), 1e-6);
    }
//    cerr << "Done with triangle " << chosenTriangle << "\n\n";
}


BOOST_AUTO_TEST_CASE( OrientationSensitivity3x3 )
{
//    cerr << "SKIPPING.\n"; return;
    cerr << "OrientationSensitivity 3x3\n";
    BOOST_CHECK(true);
    
    Rect3d bounds(0,0,0,3,3,1);
    Rect3i yeeCells(0,0,0,2,2,0);
    
    Rect3d bigBlock(0,0,0,3,3,1);
    Rect3d smallBlock(1.4, 0, 0, 3, 1.4, 1);
    
    // For the two-block test:
    // The moving vertices of the small block have ID numbers
    //  8, 10, 11, 12, 14, 15
//    std::vector<int> movingVertices = { 8, 10, 11, 12, 14, 15 };

    // For the one block test let's just move all the vertices.
    std::vector<int> movingVertices = {0, 1, 2, 3, 4, 5, 6, 7};
    
    std::vector<int> measureCell = {1, 1};
    
//    Assembly assembly;
////    TestingGoodies::twoBoxAssembly(assembly, bigBlock, smallBlock);
//    TestingGoodies::oneBoxAssembly(assembly, smallBlock);
//
//    vector<TrackedPolyhedron*> polyhedra = TestingGoodies::listPolyhedra(
//        assembly.polyhedra());
//
//    vector<vector<SimpleMesh::Triangle> > meshes = getPolyhedra(polyhedra);
    
    
    std::vector<SimpleMesh::ControlVertex> controlVertices, controlVertices2;
    std::vector<SimpleMesh::Triangle> triangles, triangles2;
    SimpleMesh::Test::oneBoxMesh(controlVertices, triangles, smallBlock);
    vector<vector<SimpleMesh::Triangle>> meshes = {triangles};
    
//    int material = 0;
    Vector2d delta(1.1e-8, 0.9e-8);
    
    for (int octant = 0; octant < 8; octant++)
    {
    //    int octant = octantE(xyz);
//        cout << "Octant " << octant << ":\n";
        
        OrientedVoxels voxels(yeeCells, realBounds(yeeCells, octant), bounds);
        voxels.initOrientations(meshes);
        
        Rect3d smallBlock2(smallBlock);
        smallBlock2.p1[0] += delta[0];
        smallBlock2.p2[1] += delta[1];
        
//        Assembly assembly2;
////        TestingGoodies::twoBoxAssembly(assembly2, bigBlock, smallBlock2);
//        TestingGoodies::oneBoxAssembly(assembly2, smallBlock2);
//        vector<TrackedPolyhedron*> polyhedra2 = TestingGoodies::listPolyhedra(
//            assembly2.polyhedra());
//
//        vector<vector<SimpleMesh::Triangle> > meshes2(polyhedra2.size());
//        for (int pp = 0; pp < polyhedra2.size(); pp++)
//        {
//            meshes2[pp] = polyhedra2[pp]->triangles();
//        }
//        SimpleMesh::initMeshSensitivity(meshes2);

        SimpleMesh::Test::oneBoxMesh(controlVertices2, triangles2, smallBlock2);
        vector<vector<SimpleMesh::Triangle>> meshes2 = {triangles2};
        
        OrientedVoxels vox2(yeeCells, realBounds(yeeCells, octant), bounds);
        vox2.initOrientations(meshes2);

        Matrix3d deltaOrientation = vox2.orientations().at(measureCell[0], measureCell[1],0)
            - voxels.orientations().at(measureCell[0], measureCell[1],0);
        Matrix3d deltaOrientationExpected;
        
        vector<Vector3d> deltaVertices = SimpleMesh::Test::subtractPositions(controlVertices2, controlVertices);
//        std::cerr << "num cvs: " << assembly.controlVertices().size() << " "
//            << assembly2.controlVertices().size() << " "
//            << deltaVertices.size() << ".\n";
        
        for (int vv = 0; vv < movingVertices.size(); vv++)
        for (int xyz = 0; xyz < 3; xyz++)
        {
            Matrix3<double> dOrientation;
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                dOrientation(i,j) = voxels.orientationJacobians(i, j)
                    .find(movingVertices[vv])->second.at(measureCell[0],measureCell[1],0)[xyz];
            }
            
            deltaOrientationExpected = deltaOrientationExpected + dOrientation * deltaVertices.at(movingVertices.at(vv))[xyz];
        }
        
        Matrix3d derivMeas = deltaOrientation / norm(delta);
        Matrix3d derivExp = deltaOrientationExpected / norm(delta);
        
        double derivErrNorm = vectorNorm(derivExp - derivMeas);
        double derivRelErr = derivErrNorm / vectorNorm(derivMeas);
        
//        cerr << "Measured:\n" << derivMeas << "\n";
//        cerr << "Expected:\n" << derivExp << "\n";
        
        if (fabs(vectorNorm(derivMeas) > 1e-9))
        {
            BOOST_CHECK_SMALL(derivRelErr, 1e-6);
        }
        else
            BOOST_CHECK_SMALL(vectorNorm(derivExp), 1e-6);
    }
}

BOOST_AUTO_TEST_CASE( FillFactorSensitivity_OneVoxel )
{
//    cerr << "SKIPPING.\n"; return;
    BOOST_CHECK(true);
    
    Rect3d bounds(0,0,0,1,1,1);
    Rect3i yeeCells(0,0,0,0,0,0);
    Rect3d bigBlock(0,0,0,1,1,1);
//    Rect3d bigBlock(100,0,0,101,1,1);
//    Rect3d smallBlock(0.25,0,0,1,0.75,1);
    Rect3d smallBlock(-1, -1, -1, 2, 2, 0.6);
    
    // The moving vertices of the small block have ID numbers
    //  8, 10, 11, 12, 14, 15
//    int movingVertices[] = { 8, 10, 11, 12, 14, 15 };
//    int numMovingVertices = 6;
//    int movingVertices[] = { 12, 13, 14, 15 };
//    int numMovingVertices = 4;
    
    std::array<int, 8> movingVertices = {0, 1, 2, 3, 4, 5, 6, 7};
    
    Vector2d delta(0.1, 0.0);
    
    int material = 0; // we'll just check this one.
    
    for (int octant = 0; octant < 8; octant++)
    {
//        cerr << "Octant " << octant << "\n";
        
//        Assembly assembly;
////        TestingGoodies::twoBoxAssembly(assembly, bigBlock, smallBlock);
//        TestingGoodies::oneBoxAssembly(assembly, smallBlock);
//
//        vector<TrackedPolyhedron*> polyhedra = TestingGoodies::listPolyhedra(
//            assembly.polyhedra());

    
        std::vector<SimpleMesh::ControlVertex> controlVertices, controlVertices2;
        std::vector<SimpleMesh::Triangle> triangles, triangles2;
        SimpleMesh::Test::oneBoxMesh(controlVertices, triangles, smallBlock);
        vector<vector<SimpleMesh::Triangle>> meshes = {triangles};
        
        vector<DynamicRLE3<double> > halfCellFillFactors, halfCellFillFactors2;
        vector<map<unsigned int, DynamicRLE3<Vector3d> > > halfCellFillJacobians;
        
        VoxelizePolyhedra vp(yeeToHalf(yeeCells), bounds);
        vp.fillFactors(meshes, halfCellFillFactors);
        vp.fillFactorJacobians(meshes, halfCellFillJacobians);
            
//        for (int mat = 0; mat < 1; mat++)
//        for (int vv = 0; vv < numMovingVertices; vv++)
//        {
//            int vId = movingVertices[vv];
//
////            int xx = 0, yy = 0, zz = 1;
//
//            if (halfCellFillJacobians[0].count(vId))
//            {
////                cerr << "Vert " << vId << " half cell sensitivity:\n";
////                cerr << halfCellFillJacobians[mat][vId] << "\n";
//            }
//        }
        
        OrientedVoxels vox(yeeCells, realBounds(yeeCells, octant), bounds);
        vox.downsampleFillFactors(halfCellFillFactors, octant);
        vox.downsampleFillJacobians(halfCellFillJacobians, octant);
//        cout << "Bounds " << vox.realBounds() << " " << vox.clipBounds() << " "
//            << vox.voxelBounds() << "\n";
//        cout << "Half fills: " << halfCellFillFactors.at(material) << "\n";
//        cerr << "Half cell Jacobian in cell (0 0 1):\n";
        
//        for (int vv = 0; vv < numMovingVertices; vv++)
//        if (halfCellFillJacobians[material].count(movingVertices[vv]))
//        {
//            double dFill = halfCellFillJacobians[material][movingVertices[vv]]
//                .markAt(0,0,1, Vector3d(0,0,0))[2];
////            cout << "component " << movingVertices[vv] << ".2 = " << dFill
////                << "\n";
//        }
        
//        cout << "Fills:\n" << vox.fillFactors(material) << "\n";
        
        Rect3d smallBlock2(smallBlock);
        smallBlock2.p2[2] += delta[0];
//        smallBlock2.p1[0] += delta[0];
//        smallBlock2.p2[1] += delta[1];
        
//        Assembly assembly2;
////        TestingGoodies::twoBoxAssembly(assembly2, bigBlock, smallBlock2);
//        TestingGoodies::oneBoxAssembly(assembly2, smallBlock2);
//        vector<TrackedPolyhedron*> polyhedra2 = TestingGoodies::listPolyhedra(
//            assembly2.polyhedra());
        
        SimpleMesh::Test::oneBoxMesh(controlVertices2, triangles2, smallBlock2);
        vector<vector<SimpleMesh::Triangle>> meshes2 = {triangles2};
        
        VoxelizePolyhedra vp2(yeeToHalf(yeeCells), bounds);
        vp2.fillFactors(meshes2, halfCellFillFactors2);
        
        OrientedVoxels vox2(yeeCells, realBounds(yeeCells, octant), bounds);
        vox2.downsampleFillFactors(halfCellFillFactors2, octant);
        //cerr << "Half fills 2:\n" << halfCellFillFactors2.at(material) << "\n";
        //cerr << "Fills 2:\n" << vox2.fillFactors(material) << "\n";
        
        double deltaFill = vox2.fillFactors().at(material).markAt(0,0,0, 0.0) -
            vox.fillFactors().at(material).markAt(0,0,0, 0.0);
        
        vector<Vector3d> deltaVertices = SimpleMesh::Test::subtractPositions(controlVertices2, controlVertices);
        
        double deltaFillExpected = 0;
        for (int vv = 0; vv < movingVertices.size(); vv++)
        for (int xyz = 0; xyz < 3; xyz++)
        {
            if (vox.fillFactorJacobians().at(material).count(movingVertices[vv]))
            {
                double dFill = vox.fillFactorJacobians().at(material)
                    .find(movingVertices[vv])->second.markAt(0,0,0, Vector3d(0,0,0))[xyz];
                
                //cerr << "component " << movingVertices[vv] << "." << xyz << " = " << dFill << "\n";
                deltaFillExpected += dFill*deltaVertices.at(movingVertices.at(vv))[xyz];
            }
        }
        
        double derivExp = deltaFillExpected/norm(delta);
        double derivMeas = deltaFill/norm(delta);
        
        if (fabs(derivMeas) < 1e-9)
        {
            BOOST_CHECK_SMALL(derivExp, 1e-9);
        }
        else
        {
            BOOST_CHECK_CLOSE(derivExp, derivMeas, 1e-5);
        }
        
        
//        cout << "Expected " << deltaFillExpected << "\n";
//        cout << "Measured " << deltaFill << "\n";
//        
//        double errVal = deltaFill - deltaFillExpected;
//        
//        cout << "Delta = " << delta << endl;
//        cout << "Error = " << errVal << "\n";
//        cout << "Relative error = " << errVal / fabs(deltaFill) << "\n";
//        cout << "Error/delta = " << errVal/delta << "\n";
//        cout << "Error/delta/delta = " << errVal/delta/delta << "\n";
    }
}

static double sDot(Vector3d lhs, Vector3d rhs)
{
    return dot(lhs, rhs);
}

BOOST_AUTO_TEST_CASE( TestCorner_FF )
{
    Rect3d bounds(0,0,0,3,3,1);
    Rect3i voxels(yeeToHalf(Rect3i(0,0,0,2,2,0)));
//    cerr << "quickRect(" << bounds << ", 'r');\n";
//    cerr << "voxels " << voxels << "\n";
    
    // These induce the weird sensitivity problem.
//    Rect3d bigBlock(0.01,0.01,0,2.99,2.99,1);
//    Rect3d smallBlock(1.5, 0.1, 0, 2.9, 1.0, 1);
    
    Rect3d bigBlock(0,0,0,3,3,1);
    Rect3d smallBlock(1.501, 0, 0, 3, 1.01, 1);
    
    double theta = 0; //M_PI/6;
    const double DELTA = 1e-8;
    double delta[] = {DELTA, DELTA};
//    cout << "DELTA " << delta[0] << " " << delta[1] << endl;
    
//    int xyz = 0;
//    int measureCell[] = {1,1};
    
    // ---------------------- INITIAL MODEL
//    Assembly assembly;
////    TestingGoodies::twoBoxAssembly(assembly, bigBlock, smallBlock, Matrix3d::rotation(theta, Vector3d::unit(2)));
//    TestingGoodies::oneBoxAssembly(assembly, smallBlock, Matrix3d::rotation(theta, Vector3d::unit(2)));
    
    std::vector<SimpleMesh::ControlVertex> controlVertices, controlVertices2;
    std::vector<SimpleMesh::Triangle> triangles, triangles2;
    SimpleMesh::Test::oneBoxMesh(controlVertices, triangles, smallBlock, Matrix3d::rotation(theta, Vector3d::unit(2)));
    vector<vector<SimpleMesh::Triangle>> meshes = {triangles};
    
//    vector<Precision::RationalFunction> materials =
//        TestingGoodies::listMaterials(assembly.polyhedra(),
//            TestingGoodies::standardMaterials());
//    vector<TrackedPolyhedron*> polyhedra = TestingGoodies::listPolyhedra(
//        assembly.polyhedra());
//    cerr << "Model:\n" << *polyhedra[0] << "\n";
    
    vector<DynamicRLE3<double> > fillFactors1, fillFactors2;
    vector<map<unsigned int, DynamicRLE3<Vector3d> > > halfCellFillJacobians;
    
    VoxelizePolyhedra vp(voxels, bounds);
    vp.fillFactors(meshes, fillFactors1);
    vp.fillFactorJacobians(meshes, halfCellFillJacobians);
    
    // ---------------------- PERTURBED MODEL
    Rect3d smallBlock2(smallBlock);
    smallBlock2.p1[0] += delta[0];
    smallBlock2.p2[1] += delta[1];
    
//    Assembly assembly2;
////    TestingGoodies::twoBoxAssembly(assembly2, bigBlock, smallBlock2,
////        Matrix3d::rotation(theta, Vector3d::unit(2)));
//    TestingGoodies::oneBoxAssembly(assembly2, smallBlock2, Matrix3d::rotation(theta, Vector3d::unit(2)));
//    vector<TrackedPolyhedron*> polyhedra2 = TestingGoodies::listPolyhedra(
//        assembly2.polyhedra());
    
    SimpleMesh::Test::oneBoxMesh(controlVertices2, triangles2, smallBlock2, Matrix3d::rotation(theta, Vector3d::unit(2)));
    vector<vector<SimpleMesh::Triangle>> meshes2 = {triangles2};
    
    VoxelizePolyhedra vp2(voxels, bounds);
    vp2.fillFactors(meshes2, fillFactors2);
    
    vector<Vector3d> deltaCV = SimpleMesh::Test::subtractPositions(controlVertices2, controlVertices);
    
    if (0)
    {
//        cout << "Delta vertices:\n";
        for (int vv = 0; vv < deltaCV.size(); vv++)
        if (norm(deltaCV[vv]) != 0)
            cout << vv << ": " << deltaCV[vv] << "\n";
    }
        
    vector<DynamicRLE3<Precision::Float> > deltaFillFactorPredicted(2);
    
    int numMaterials = meshes.size();
    for (int matl = 0; matl < numMaterials; matl++)
    {
        DynamicRLE3<double> ffPredicted(fillFactors1[matl]);
        DynamicRLE3<double> deltaFF = fillFactors2[matl] - fillFactors1[matl];
        
        map<unsigned int, DynamicRLE3<Vector3d> > & DF =
            halfCellFillJacobians[matl];
        
        for (int vv = 0; vv < deltaCV.size(); vv++)
        {
            DynamicRLE3<double> dotProd;
            RLE::scalarTransform(DF[vv], deltaCV[vv], dotProd, sDot);
            ffPredicted = ffPredicted + dotProd;
//            cerr << "v" << vv << ": dot =\n" << dotProd << "\n";
        }
            
        DynamicRLE3<double> dfdv = (ffPredicted - fillFactors1[matl]) / DELTA;
        DynamicRLE3<double> dfdv_meas = deltaFF/DELTA;
        
//        cerr << "Derivative measured:\n" << deltaFF/DELTA << "\n";
//        cerr << "Derivative predicted:\n" << dfdv << "\n";
        for (int xx = voxels.p1[0]; xx <= voxels.p2[0]; xx++)
        for (int yy = voxels.p1[1]; yy <= voxels.p2[1]; yy++)
        for (int zz = voxels.p1[2]; zz <= voxels.p2[2]; zz++)
        {
            if (dfdv.markAt(xx,yy,zz,0.0) != 0.0)
            {
                BOOST_CHECK_CLOSE(dfdv.markAt(xx,yy,zz,0.0),
                    dfdv_meas.markAt(xx,yy,zz,0.0), 1e-4);
            }
            else
            {
                BOOST_CHECK_SMALL(deltaFF.markAt(xx,yy,zz,0.0), 1e-6);
            }
        }
    }
}




