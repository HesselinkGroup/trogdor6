/*
 *  Extents.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/30/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "Extents.h"

#include <cassert>
#include <stdexcept>

using namespace YeeUtilities;
using namespace std;

PhysicalExtents::
PhysicalExtents() :
    mPhysicalYeeCells(0,0,0,0,0,0),
    mNonPMLHalfCells(0,0,0,0,0,0),
    mBoundaries(6, kPECBoundary)
{
//    for (int nn = 0; nn < 6; nn++)
//        mBoundaries[nn] = kPECBoundary;
}

//PhysicalExtents::
//PhysicalExtents(const PhysicalExtents & copyMe) :
//    mPhysicalYeeCells(copyMe.mPhysicalYeeCells),
//    mNonPMLHalfCells(copyMe.mNonPMLHalfCells)
//{
//    for (int nn = 0; nn < 6; nn++)
//        mBoundaries[nn] = copyMe.mBoundaries[nn];
//}

PhysicalExtents::
PhysicalExtents(Rect3i physicalYeeCells,
    Rect3i nonPMLHalfCells) :
    mPhysicalYeeCells(physicalYeeCells),
    mNonPMLHalfCells(nonPMLHalfCells),
    mBoundaries(6),
    mTransverseAxis(-1)
{
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (mPhysicalYeeCells.num(xyz) > 1)
        {
            mBoundaries.at(2*xyz) = kPECBoundary;
            mBoundaries.at(2*xyz+1) = kPECBoundary;
        }
        else
        {
            mBoundaries.at(2*xyz) = kTranslationSymmetricBoundary;
            mBoundaries.at(2*xyz+1) = kTranslationSymmetricBoundary;
            mTransverseAxis = xyz;
        }
    }
}

PhysicalExtents::
PhysicalExtents(Rect3i physicalYeeCells,
    Rect3i nonPMLHalfCells,
    const vector<BoundaryType> & physicalBounds) :
    mPhysicalYeeCells(physicalYeeCells),
    mNonPMLHalfCells(nonPMLHalfCells),
    mBoundaries(6),
    mTransverseAxis(-1)
{
    if (physicalBounds.size() != 6)
        throw(std::logic_error("Must specify six boundary types"));
    mBoundaries = physicalBounds;
    
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (mPhysicalYeeCells.num(xyz) == 1)
        {
            mTransverseAxis = xyz;
        }
    }
}

//PhysicalExtents & PhysicalExtents::
//operator=(const PhysicalExtents & rhs)
//{
//    if (this == &rhs)
//        return *this;
//    mPhysicalYeeCells = rhs.mPhysicalYeeCells;
//    mNonPMLHalfCells = rhs.mNonPMLHalfCells;
//    for (int nn = 0; nn < 6; nn++)
//        mBoundaries[nn] = rhs.mBoundaries[nn];
//    return *this;
//}

Vector3b PhysicalExtents::
dimensions() const
{
    Vector3b dims;
    for (int xyz = 0; xyz < 3; xyz++)
    {
        dims[xyz] = (physicalBoundaries().at(2*xyz) != kTranslationSymmetricBoundary);
    }
    return dims;
}

bool PhysicalExtents::
hasPML(int face) const
{
    if (face%2 == 0) // low X, low Y or low Z
    {
        return (yeeToHalf(mPhysicalYeeCells).p1[face/2] !=
            mNonPMLHalfCells.p1[face/2]);
    }
    else
    {
        return (yeeToHalf(mPhysicalYeeCells).p2[face/2] !=
            mNonPMLHalfCells.p2[face/2]);
    }
}

Rect3i PhysicalExtents::
pmlHalfCells(int face) const
{
    assert(hasPML(face));
    Rect3i pmlHalfCells = yeeToHalf(mPhysicalYeeCells);
    if (face%2 == 0) // left, bottom, front side
        pmlHalfCells.p2[face/2] = mNonPMLHalfCells.p1[face/2]-1;
    else
        pmlHalfCells.p1[face/2] = mNonPMLHalfCells.p2[face/2]+1;

    return pmlHalfCells;
}

void PhysicalExtents::
physicalYeeCells(const Rect3i & r)
{
    mPhysicalYeeCells = r;
}

void PhysicalExtents::
nonPMLHalfCells(const Rect3i & r)
{
    mNonPMLHalfCells = r;
}

void PhysicalExtents::
physicalBoundaries(const vector<BoundaryType> & bounds)
{
    if (bounds.size() != 6)
        throw(std::logic_error("Must specify six boundaries"));
    mBoundaries = bounds;
}

bool operator==(const PhysicalExtents & lhs, const PhysicalExtents & rhs)
{
    if (lhs.physicalYeeCells() != rhs.physicalYeeCells() ||
        lhs.nonPMLHalfCells() != rhs.nonPMLHalfCells())
        return 0;
    for (int nn = 0; nn < 6; nn++)
        if (lhs.physicalBoundaries()[nn] != rhs.physicalBoundaries()[nn])
            return 0;
    return 1;
}

NodeExtents::
NodeExtents() :
    mPhysicalYeeCells(0,0,0,0,0,0),
    mNodeYeeCells(0,0,0,0,0,0),
    mCalcHalfCells(0,0,0,0,0,0),
    mGhostHalfCells(0,0,0,0,0,0),
    mNonPMLHalfCells(0,0,0,0,0,0),
    mNodeBoundaries(6, kPECBoundary)
{
}

NodeExtents::
NodeExtents(Rect3i physicalYeeCells, Rect3i calcHalfCells,
    Rect3i ghostHalfCells, Rect3i nonPMLHalfCells) :
    mPhysicalYeeCells(physicalYeeCells),
    mCalcHalfCells(calcHalfCells),
    mGhostHalfCells(ghostHalfCells),
    mNonPMLHalfCells(nonPMLHalfCells),
    mNodeBoundaries(6, kPECBoundary)
{
}

NodeExtents::
NodeExtents(const PhysicalExtents & gridExtents, Rect3i nodeYeeCells,
    Vector3i nodeCoordinates, Vector3i numNodes) :
    mPhysicalYeeCells(gridExtents.physicalYeeCells()),
    mNodeYeeCells(nodeYeeCells),
    mNodeBoundaries(6, kPECBoundary)
{
    mCalcHalfCells = yeeToHalf(nodeYeeCells);
    mGhostHalfCells = yeeToHalf(nodeYeeCells);
    unsigned int numGhostHalfCells = 1;
    
    getNodeBoundaries(gridExtents.physicalBoundaries(), mNodeBoundaries,
        nodeCoordinates, numNodes);
    
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (nodeBoundaries().at(2*xyz) == kMPIBoundary)
        {
            mGhostHalfCells.p1[xyz] -= numGhostHalfCells;
            mCalcHalfCells.p1[xyz] -= (numGhostHalfCells-1);
        }
        else if (nodeBoundaries().at(2*xyz) == kPeriodicBoundary)
        {
            mGhostHalfCells.p1[xyz] -= 1;
        }
        else if (nodeBoundaries().at(2*xyz) == kPECBoundary)
            mCalcHalfCells.p1[xyz] += 1;
        else if (nodeBoundaries().at(2*xyz) == kPMCBoundary)
            mCalcHalfCells.p1[xyz] += 2;
        else if (nodeBoundaries().at(2*xyz) == kPMLBoundary)
            mCalcHalfCells.p1[xyz] += 1;
        else if (nodeBoundaries().at(2*xyz) == kTranslationSymmetricBoundary)
            {}
        
        if (nodeBoundaries().at(2*xyz+1) == kMPIBoundary)
        {
            mGhostHalfCells.p2[xyz] += numGhostHalfCells;
            mCalcHalfCells.p2[xyz] += (numGhostHalfCells-1);
        }
        else if (nodeBoundaries().at(2*xyz+1) == kPeriodicBoundary)
        {
            mGhostHalfCells.p2[xyz] += 1;
        }
        else if (nodeBoundaries().at(2*xyz+1) == kPECBoundary)
            mCalcHalfCells.p2[xyz] -= 2;
        else if (nodeBoundaries().at(2*xyz+1) == kPMCBoundary)
            mCalcHalfCells.p2[xyz] -= 1;
        else if (nodeBoundaries().at(2*xyz+1) == kPMLBoundary)
            mCalcHalfCells.p2[xyz] -= 1;
        else if (nodeBoundaries().at(2*xyz+1) == kTranslationSymmetricBoundary)
            {}
    }
    
    mNonPMLHalfCells = ghostHalfCells();
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (gridExtents.physicalBoundaries().at(2*xyz) == kPMLBoundary &&
            gridExtents.nonPMLHalfCells().p1[xyz] > mNonPMLHalfCells.p1[xyz])
            mNonPMLHalfCells.p1[xyz] = gridExtents.nonPMLHalfCells().p1[xyz];
        if (gridExtents.physicalBoundaries().at(2*xyz+1) == kPMLBoundary &&
            gridExtents.nonPMLHalfCells().p2[xyz] < mNonPMLHalfCells.p2[xyz])
            mNonPMLHalfCells.p2[xyz] = gridExtents.nonPMLHalfCells().p2[xyz];
    }
}

bool NodeExtents::
isEntirelyPML() const
{
    return (!vec_ge(mNonPMLHalfCells.p2, mNonPMLHalfCells.p1));
}

Vector3b NodeExtents::
dimensions() const
{
    Vector3b dims;
    for (int xyz = 0; xyz < 3; xyz++)
        dims[xyz] = (nodeBoundaries().at(2*xyz) != kTranslationSymmetricBoundary);
    return dims;
}

int NodeExtents::
numGhostHalfCells(int side) const
{
    if (side%2 == 0) // low side
        return calcHalfCells().p1[side/2] - ghostHalfCells().p1[side/2];
    else
        return ghostHalfCells().p2[side/2] - calcHalfCells().p2[side/2];
}

Rect3i NodeExtents::
ghostReadHalfCells(int side) const
{
    Rect3i insideSlab = ghostHalfCells();
    insideSlab.p1[side/2] = calcHalfCells().p1[side/2];
    insideSlab.p2[side/2] = calcHalfCells().p2[side/2];
    
    int depth = numGhostHalfCells(side);
    
    if (side%2 == 0) // low side
    {
        insideSlab.p2[side/2] = insideSlab.p1[side/2] + depth - 1;
    }
    else
    {
        insideSlab.p1[side/2] = insideSlab.p2[side/2] - depth + 1;
    }
    
    return insideSlab;
}

Rect3i NodeExtents::
ghostWriteHalfCells(int side) const
{
    Rect3i outsideSlab = ghostHalfCells();
    
    int depth = numGhostHalfCells(side);
    
    if (side%2 == 0) // low side
    {
        outsideSlab.p2[side/2] = outsideSlab.p1[side/2] + depth - 1;
    }
    else
    {
        outsideSlab.p1[side/2] = outsideSlab.p2[side/2] - depth + 1;
    }
    
    return outsideSlab;
}

void NodeExtents::
physicalYeeCells(const Rect3i & r)
{
    mPhysicalYeeCells = r;
}

void NodeExtents::
calcHalfCells(const Rect3i & r)
{
    mCalcHalfCells = r;
}

void NodeExtents::
ghostHalfCells(const Rect3i & r)
{
    mGhostHalfCells = r;
}

void NodeExtents::
nonPMLHalfCells(const Rect3i & r)
{
    mNonPMLHalfCells = r;
}

void NodeExtents::
nodeBoundaries(const vector<BoundaryType> & bounds)
{
    if (bounds.size() != 6)
        throw(std::logic_error("Must specify six boundary types"));
    mNodeBoundaries = bounds;
}

bool operator==(const NodeExtents & lhs, const NodeExtents & rhs)
{
    if (lhs.physicalYeeCells() != rhs.physicalYeeCells() ||
        lhs.nodeYeeCells() != rhs.nodeYeeCells() ||
        lhs.calcHalfCells() != rhs.calcHalfCells() ||
        lhs.ghostHalfCells() != rhs.ghostHalfCells() ||
        lhs.nonPMLHalfCells() != rhs.nonPMLHalfCells())
        return 0;
    for (int nn = 0; nn < 6; nn++)
        if (lhs.nodeBoundaries()[nn] != rhs.nodeBoundaries()[nn])
            return 0;
    return 1;
}

void
getNodeBoundaries(const vector<BoundaryType> & physicalBoundaries,
    vector<BoundaryType> & outNodeBoundaries, Vector3i nodeCoordinates,
    Vector3i numNodes)
{
    assert(vec_ge(numNodes, Vector3i(1,1,1)));
    outNodeBoundaries.resize(6);
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (numNodes[xyz] == 1)
        {
            outNodeBoundaries.at(2*xyz) = physicalBoundaries.at(2*xyz);
            outNodeBoundaries.at(2*xyz+1) = physicalBoundaries.at(2*xyz+1);
        }
        else
        {
            if (physicalBoundaries.at(2*xyz) == kPeriodicBoundary ||
                physicalBoundaries.at(2*xyz+1) == kPeriodicBoundary)
            {
                assert(physicalBoundaries.at(2*xyz) ==
                    physicalBoundaries.at(2*xyz+1));
                outNodeBoundaries.at(2*xyz) = kMPIBoundary;
                outNodeBoundaries.at(2*xyz+1) = kMPIBoundary;
            }
            else
            {
                if (nodeCoordinates[xyz] == 0)
                    outNodeBoundaries.at(2*xyz) = physicalBoundaries.at(2*xyz);
                else
                    outNodeBoundaries.at(2*xyz) = kMPIBoundary;
                
                if (nodeCoordinates[xyz] == numNodes[xyz]-1)
                    outNodeBoundaries.at(2*xyz+1) = physicalBoundaries.at(2*xyz+1);
                else
                    outNodeBoundaries.at(2*xyz+1) = kMPIBoundary;
            }
        }
    }
}




