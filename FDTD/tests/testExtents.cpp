/*
 *  testExtents.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/30/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test3D

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <vector>
#include "Extents.h"
#include "YeeUtilities.h"
#include "geometry.h"

using namespace std;
using namespace YeeUtilities;

BOOST_AUTO_TEST_CASE( TestPMLExtents )
{
    Rect3i allCells(-30, -20, 10, 40, -10, 20); // num 71 11 11
    int d[] = { 5, 5, 0, 10, 0, 5 }; // PML depth
    Rect3i nonPML = inset(allCells, d[0], d[1], d[2], d[3], d[4], d[5]);
    
    PhysicalExtents grid(allCells, yeeToHalf(nonPML));
    
    BOOST_CHECK(grid.hasPML(0));
    BOOST_CHECK(grid.hasPML(1));
    BOOST_CHECK(grid.hasPML(2));
    BOOST_CHECK(!grid.hasPML(3));
    BOOST_CHECK(!grid.hasPML(4));
    BOOST_CHECK(grid.hasPML(5));
    
    BOOST_CHECK_EQUAL(grid.pmlHalfCells(0).num(), Vector3i(2*d[0], 22, 22));
    BOOST_CHECK_EQUAL(grid.pmlHalfCells(1).num(), Vector3i(2*d[3], 22, 22));
    BOOST_CHECK_EQUAL(grid.pmlHalfCells(2).num(), Vector3i(142, 2*d[1], 22));
    BOOST_CHECK_EQUAL(grid.pmlHalfCells(5).num(), Vector3i(142, 22, 2*d[5]));
}

BOOST_AUTO_TEST_CASE( TestYeeUtilities )
{
    // These tests got dropped here because I discovered that YeeUtilities
    // incorrectly did the halfToYee conversion for cells with negative
    // coordinates.  I had never used it for negative coordinates really... (-:
    BOOST_CHECK_EQUAL(halfToYee(Vector3i(0,0,0)), Vector3i(0,0,0));
    BOOST_CHECK_EQUAL(halfToYee(Vector3i(-1,0,0)), Vector3i(-1,0,0));
    BOOST_CHECK_EQUAL(halfToYee(Vector3i(-2,0,0)), Vector3i(-1,0,0));
    BOOST_CHECK_EQUAL(halfToYee(Vector3i(-5,1,4)), Vector3i(-3, 0, 2));
    
    Rect3i halfRect(-200, -40, 20, 81, -3, 1);
    Rect3i yeeRect(-100, -20, 10, 40, -2, 0);
    BOOST_CHECK_EQUAL(halfToYee(halfRect), yeeRect);
    BOOST_CHECK_EQUAL(yeeToHalf(halfToYee(halfRect)), halfRect);
    
    BOOST_CHECK_EQUAL(yeeToHalf(Rect3i(-10,-10,0,10,10,0)),
        Rect3i(-20, -20, 0, 21, 21, 1));
    BOOST_CHECK_EQUAL(yeeToHalf(Rect3i(-5,-5,-5,-5,-5,-5)),
        Rect3i(-10,-10,-10,-9,-9,-9));
}

BOOST_AUTO_TEST_CASE( TestGhostUtilities )
{
    Rect3i physicalYee(0,0,0,10,10,10);
    Rect3i calcHalf(0,0,0,19,19,19);
    Rect3i ghostHalf(-1,-2,-3,20,21,22);
    Rect3i nonPMLHalf(calcHalf);
    
    NodeExtents ext(physicalYee, calcHalf, ghostHalf, nonPMLHalf);
    
    // Test depth of ghost slab
    BOOST_CHECK_EQUAL(1, ext.numGhostHalfCells(0)); // X
    BOOST_CHECK_EQUAL(1, ext.numGhostHalfCells(1));
    BOOST_CHECK_EQUAL(2, ext.numGhostHalfCells(2)); // Y
    BOOST_CHECK_EQUAL(2, ext.numGhostHalfCells(3));
    BOOST_CHECK_EQUAL(3, ext.numGhostHalfCells(4)); // Z
    BOOST_CHECK_EQUAL(3, ext.numGhostHalfCells(5));
    
    // Test the read slabs
    BOOST_CHECK(ext.ghostReadHalfCells(0) == Rect3i(0,-2,-3,0,21,22));
    BOOST_CHECK(ext.ghostReadHalfCells(1) == Rect3i(19,-2,-3,19,21,22));
    BOOST_CHECK(ext.ghostReadHalfCells(2) == Rect3i(-1,0,-3,20,1,22));
    BOOST_CHECK(ext.ghostReadHalfCells(3) == Rect3i(-1,18,-3,20,19,22));
    BOOST_CHECK(ext.ghostReadHalfCells(4) == Rect3i(-1,-2,0,20,21,2));
    BOOST_CHECK(ext.ghostReadHalfCells(5) == Rect3i(-1,-2,17,20,21,19));
    
    // Test the write slabs
    BOOST_CHECK(ext.ghostWriteHalfCells(0) == Rect3i(-1,-2,-3,-1,21,22));
    BOOST_CHECK(ext.ghostWriteHalfCells(1) == Rect3i(20,-2,-3,20,21,22));
    BOOST_CHECK(ext.ghostWriteHalfCells(2) == Rect3i(-1,-2,-3,20,-1,22));
    BOOST_CHECK(ext.ghostWriteHalfCells(3) == Rect3i(-1,20,-3,20,21,22));
    BOOST_CHECK(ext.ghostWriteHalfCells(4) == Rect3i(-1,-2,-3,20,21,-1));
    BOOST_CHECK(ext.ghostWriteHalfCells(5) == Rect3i(-1,-2,20,20,21,22));
}


BOOST_AUTO_TEST_CASE( TestOneNodeX )
{
    vector<BoundaryType> periodic(6, kTranslationSymmetricBoundary);
    periodic[0] = kPeriodicBoundary;
    periodic[1] = kPeriodicBoundary;
    
    vector<BoundaryType> pml(6, kTranslationSymmetricBoundary);
    pml[0] = kPMLBoundary;
    pml[1] = kPMLBoundary;
    
    vector<BoundaryType> pec(6, kTranslationSymmetricBoundary);
    pec[0] = kPECBoundary;
    pec[1] = kPECBoundary;
    
    vector<BoundaryType> pmc(6, kTranslationSymmetricBoundary);
    pmc[0] = kPMCBoundary;
    pmc[1] = kPMCBoundary;
    
    Rect3i allCells(-100,0,0,100,0,0);
    BOOST_CHECK_EQUAL(yeeToHalf(allCells), Rect3i(-200, 0, 0, 201, 1, 1));
    
    PhysicalExtents gridPeriodic(allCells, yeeToHalf(allCells), periodic);
    PhysicalExtents gridPML(allCells, yeeToHalf(allCells), pml);
    PhysicalExtents gridPEC(allCells, yeeToHalf(allCells), pec);
    PhysicalExtents gridPMC(allCells, yeeToHalf(allCells), pmc);
    
    BOOST_CHECK(gridPeriodic.dimensions() == Vector3b(1,0,0));
    BOOST_CHECK(gridPML.dimensions() == Vector3b(1,0,0));
    BOOST_CHECK(gridPEC.dimensions() == Vector3b(1,0,0));
    BOOST_CHECK(gridPMC.dimensions() == Vector3b(1,0,0));
    
    // A single node in a periodic grid should add one ghost point.
    NodeExtents solePeriodicNode(gridPeriodic, allCells);
    BOOST_CHECK_EQUAL(solePeriodicNode.physicalYeeCells(), allCells);
    BOOST_CHECK_EQUAL(solePeriodicNode.nodeYeeCells(), allCells);
    BOOST_CHECK_EQUAL(solePeriodicNode.calcHalfCells(), yeeToHalf(allCells));
    BOOST_CHECK_EQUAL(solePeriodicNode.ghostHalfCells(),
        inset(yeeToHalf(allCells), -1, 0, 0, -1, 0, 0));
    BOOST_CHECK_EQUAL(solePeriodicNode.allocatedYeeCells(),
        halfToYee(solePeriodicNode.ghostHalfCells()));
    BOOST_CHECK_EQUAL(solePeriodicNode.nonPMLHalfCells(), 
        solePeriodicNode.ghostHalfCells());
    BOOST_CHECK_EQUAL(solePeriodicNode.allocatedYeeCells(),
        inset(allCells, -1, 0, 0, -1, 0, 0));
    BOOST_CHECK(solePeriodicNode.dimensions() == Vector3b(1,0,0));
    
    // A single node on a grid with PML should inset the calc region by one.
    // The non-PML region was actually set to be the entire grid so it will
    // still appear that there is no PML (a zero-depth PML is just a PEC or PMC
    // boundary condition).
    NodeExtents solePMLNode(gridPML, allCells);
    BOOST_CHECK_EQUAL(solePMLNode.physicalYeeCells(), allCells);
    BOOST_CHECK_EQUAL(solePMLNode.nodeYeeCells(), allCells);
    BOOST_CHECK_EQUAL(solePMLNode.calcHalfCells(),
        inset(yeeToHalf(allCells), 1, 0, 0, 1, 0, 0));
    BOOST_CHECK_EQUAL(solePMLNode.ghostHalfCells(), yeeToHalf(allCells));
    BOOST_CHECK_EQUAL(solePMLNode.nonPMLHalfCells(), 
        solePMLNode.ghostHalfCells());
    BOOST_CHECK_EQUAL(solePMLNode.allocatedYeeCells(), allCells);
    BOOST_CHECK(solePMLNode.dimensions() == Vector3b(1,0,0));
    
    // A single node on a grid with PEC boundaries will inset the calc region
    // by one half cell on the left and two half cells on the right.
    NodeExtents solePECNode(gridPEC, allCells);
    BOOST_CHECK_EQUAL(solePECNode.physicalYeeCells(), allCells);
    BOOST_CHECK_EQUAL(solePECNode.nodeYeeCells(), allCells);
//    BOOST_CHECK_EQUAL(clip(grid.calcHalfCells(), yeeToHalf(allCells)),
//        Rect3i(-200, 0, 0, 201, 1, 1));
    BOOST_CHECK_EQUAL(solePECNode.calcHalfCells(),
        inset(yeeToHalf(allCells), 1, 0, 0, 2, 0, 0));
    BOOST_CHECK_EQUAL(solePECNode.ghostHalfCells(), yeeToHalf(allCells));
    BOOST_CHECK_EQUAL(solePECNode.nonPMLHalfCells(), 
        solePECNode.ghostHalfCells());
    BOOST_CHECK_EQUAL(solePECNode.allocatedYeeCells(), allCells);
    BOOST_CHECK(solePECNode.dimensions() == Vector3b(1,0,0));
    
    // A single node on a grid with PMC boundaries will inset the calc region
    // by two half cells on the left and one half cell on the right.
    NodeExtents solePMCNode(gridPMC, allCells);
    BOOST_CHECK_EQUAL(solePMCNode.physicalYeeCells(), allCells);
    BOOST_CHECK_EQUAL(solePMCNode.nodeYeeCells(), allCells);
    BOOST_CHECK_EQUAL(solePMCNode.calcHalfCells(),
        inset(yeeToHalf(allCells), 2, 0, 0, 1, 0, 0));
    BOOST_CHECK_EQUAL(solePMCNode.ghostHalfCells(), yeeToHalf(allCells));
    BOOST_CHECK_EQUAL(solePMCNode.nonPMLHalfCells(), 
        solePMCNode.ghostHalfCells());
    BOOST_CHECK_EQUAL(solePMCNode.allocatedYeeCells(), allCells);
    BOOST_CHECK(solePMCNode.dimensions() == Vector3b(1,0,0));
}

BOOST_AUTO_TEST_CASE( TestManyNodesY )
{
//    int numNodes = 3;
    Rect3i allCells(0, -30, 0, 0, 30, 0);
    Rect3i partitions[3] = { Rect3i(0, -30, 0, 0, -11, 0),
        Rect3i(0, -10, 0, 0, 9, 0),
        Rect3i(0, 10, 0, 0, 30, 0) };
    
    Rect3i nonPML(0, -30, 0, 0, 20, 0); // ten cells of PML at +Y
    
    vector<BoundaryType> pecPML(6, kTranslationSymmetricBoundary);
    pecPML[2] = kPECBoundary;
    pecPML[3] = kPMLBoundary;
    
    PhysicalExtents gridPECPML(allCells, yeeToHalf(nonPML), pecPML);
    BOOST_CHECK(gridPECPML.dimensions() == Vector3b(0,1,0));
    
    NodeExtents nodes[3];
    for (int mm = 0; mm < 3; mm++)
    {
        BOOST_TEST_MESSAGE(mm);
        nodes[mm] = NodeExtents(gridPECPML, partitions[mm],
            Vector3i(0, mm, 0), Vector3i(1, 3, 1));
        BOOST_TEST_MESSAGE(mm << " again");
        BOOST_CHECK_EQUAL(nodes[mm].nodeYeeCells(), partitions[mm]);
        BOOST_CHECK(!nodes[mm].isEntirelyPML());
        BOOST_CHECK(nodes[mm].dimensions() == Vector3b(0,1,0));
    }
    
    // The bottommost node should have a ghost cell on top but none beneath.
    // The calc half cells should be inset by one from the bottom boundary to
    // implement a PEC boundary.
    BOOST_CHECK_EQUAL(nodes[0].calcHalfCells(),
        inset(yeeToHalf(partitions[0]), 0, 1, 0, 0, 0, 0));
    BOOST_CHECK_EQUAL(nodes[0].ghostHalfCells(),
        inset(yeeToHalf(partitions[0]), 0, 0, 0, 0, -1, 0));
    BOOST_CHECK_EQUAL(nodes[0].nonPMLHalfCells(), nodes[0].ghostHalfCells());
    BOOST_CHECK_EQUAL(nodes[0].allocatedYeeCells(),
        inset(partitions[0], 0, 0, 0, 0, -1, 0));
    
    // The middle node has ghost cells for MPI above and beneath.
    BOOST_CHECK_EQUAL(nodes[1].calcHalfCells(), yeeToHalf(partitions[1]));
    BOOST_CHECK_EQUAL(nodes[1].ghostHalfCells(),
        inset(yeeToHalf(partitions[1]), 0, -1, 0, 0, -1, 0));
    BOOST_CHECK_EQUAL(nodes[1].nonPMLHalfCells(), nodes[1].ghostHalfCells());
    BOOST_CHECK_EQUAL(nodes[1].allocatedYeeCells(),
        inset(partitions[1], 0, -1, 0, 0, -1, 0));
    
    // The top node has a ghost cell beneath and a PMC boundary on top (behind
    // the PML layer).
    BOOST_CHECK_EQUAL(nodes[2].calcHalfCells(),
        inset(yeeToHalf(partitions[2]), 0, 0, 0, 0, 1, 0));
    BOOST_CHECK_EQUAL(nodes[2].ghostHalfCells(),
        inset(yeeToHalf(partitions[2]), 0, -1, 0, 0, 0, 0));
    BOOST_CHECK_EQUAL(nodes[2].nonPMLHalfCells(),
        clip(nodes[2].ghostHalfCells(), yeeToHalf(nonPML)));
    BOOST_CHECK_EQUAL(nodes[2].allocatedYeeCells(),
        inset(partitions[2], 0, -1, 0, 0, 0, 0));
}

BOOST_AUTO_TEST_CASE( TestDeepPML )
{
//    int numNodes = 3;
    Rect3i allCells(0, -30, 0, 0, 30, 0);
    Rect3i partitions[3] = { Rect3i(0, -30, 0, 0, -11, 0),
        Rect3i(0, -10, 0, 0, 9, 0),
        Rect3i(0, 10, 0, 0, 30, 0) };
    
    Rect3i nonPML(0, -30, 0, 0, 0, 0); // 30 cells of PML at +Y
    
    vector<BoundaryType> pecPML(6, kTranslationSymmetricBoundary);
    pecPML[2] = kPECBoundary;
    pecPML[3] = kPMLBoundary;
    
    PhysicalExtents gridPECPML(allCells, yeeToHalf(nonPML), pecPML);
    BOOST_CHECK(!gridPECPML.hasPML(0));
    BOOST_CHECK(!gridPECPML.hasPML(1));
    BOOST_CHECK(!gridPECPML.hasPML(2));
    BOOST_CHECK(gridPECPML.hasPML(3));
    BOOST_CHECK(!gridPECPML.hasPML(4));
    BOOST_CHECK(!gridPECPML.hasPML(5));
    // the PML occupies a 2x2 half cell column, 60 half cells long along Y.
    BOOST_CHECK_EQUAL(gridPECPML.pmlHalfCells(3).num(), Vector3i(2,60,2));
    BOOST_CHECK(gridPECPML.dimensions() == Vector3b(0,1,0));
    
    NodeExtents nodes[3];
    for (int mm = 0; mm < 3; mm++)
    {
        nodes[mm] = NodeExtents(gridPECPML, partitions[mm],
            Vector3i(0, mm, 0), Vector3i(1, 3, 1));
        BOOST_CHECK_EQUAL(nodes[mm].nodeYeeCells(), partitions[mm]);
        BOOST_CHECK(nodes[mm].dimensions() == Vector3b(0,1,0));
    }
    
    // The bottommost node should have a ghost cell on top but none beneath.
    // The calc half cells should be inset by one from the bottom boundary to
    // implement a PEC boundary.
    BOOST_CHECK_EQUAL(nodes[0].calcHalfCells(),
        inset(yeeToHalf(partitions[0]), 0, 1, 0, 0, 0, 0));
    BOOST_CHECK_EQUAL(nodes[0].ghostHalfCells(),
        inset(yeeToHalf(partitions[0]), 0, 0, 0, 0, -1, 0));
    BOOST_CHECK(!nodes[0].isEntirelyPML());
    BOOST_CHECK_EQUAL(nodes[0].nonPMLHalfCells(), nodes[0].ghostHalfCells());
    BOOST_CHECK_EQUAL(nodes[0].allocatedYeeCells(),
        inset(partitions[0], 0, 0, 0, 0, -1, 0));
    
    // The middle node has ghost cells for MPI above and beneath; furthermore
    // it's got some PML in it.
    Rect3i middleNonPMLHalfCells = intersection(
        nodes[1].ghostHalfCells(), yeeToHalf(nonPML));
    BOOST_CHECK_EQUAL(nodes[1].calcHalfCells(), yeeToHalf(partitions[1]));
    BOOST_CHECK_EQUAL(nodes[1].ghostHalfCells(),
        inset(yeeToHalf(partitions[1]), 0, -1, 0, 0, -1, 0));
    BOOST_CHECK(!nodes[1].isEntirelyPML());
    BOOST_CHECK_EQUAL(nodes[1].nonPMLHalfCells(), middleNonPMLHalfCells);
    BOOST_CHECK_EQUAL(nodes[1].allocatedYeeCells(),
        inset(partitions[1], 0, -1, 0, 0, -1, 0));
    
    // The top node has a ghost cell beneath and a PMC boundary on top (behind
    // the PML layer).
    BOOST_CHECK_EQUAL(nodes[2].calcHalfCells(),
        inset(yeeToHalf(partitions[2]), 0, 0, 0, 0, 1, 0));
    BOOST_CHECK_EQUAL(nodes[2].ghostHalfCells(),
        inset(yeeToHalf(partitions[2]), 0, -1, 0, 0, 0, 0));
    BOOST_CHECK(nodes[2].isEntirelyPML());
//    BOOST_CHECK_EQUAL(nodes[2].nonPMLHalfCells(),
//        clip(nodes[2].ghostHalfCells(), yeeToHalf(nonPML))); // NOT TRUE
    BOOST_CHECK_EQUAL(nodes[2].allocatedYeeCells(),
        inset(partitions[2], 0, -1, 0, 0, 0, 0));
}



