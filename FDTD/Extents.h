/*
 *  Extents.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/30/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _EXTENTS_
#define _EXTENTS_

#include "geometry.h"
#include "YeeUtilities.h"
#include <vector>

enum BoundaryType
{
    kPECBoundary,
    kPMCBoundary,
    kPeriodicBoundary,
    kMPIBoundary,
    kPMLBoundary,
    kTranslationSymmetricBoundary
};

class PhysicalExtents
{
public:
    PhysicalExtents();
//    PhysicalExtents(const PhysicalExtents & copyMe);
    PhysicalExtents(Rect3i physicalYeeCells, Rect3i nonPMLHalfCells);
    PhysicalExtents(Rect3i physicalYeeCells, Rect3i nonPMLHalfCells,
        const std::vector<BoundaryType> & physicalBoundaries);
    
//    PhysicalExtents & operator=(const PhysicalExtents & rhs);
    
    const Rect3i & physicalYeeCells() const { return mPhysicalYeeCells; }
    const Rect3i & nonPMLHalfCells() const { return mNonPMLHalfCells; }
    const std::vector<BoundaryType> physicalBoundaries() const
        { return mBoundaries; }
    // for 2D simulations, which way is out of plane?
    int transverseAxis() const { return mTransverseAxis; }
    /**
     * Returns true for all directions that don't have translation-symmetric
     * boundaries.
     */
    Vector3b dimensions() const;
    bool hasPML(int face) const;
    
    /**
     * Return the size of the PML zone on the given face.  It extends from just
     * outside mNonPMLHalfCells to the last half cell of mPhysicalYeeCells.
     * This extends outside the calc region--the idea is that any cell in the
     * universe that is not in mNonPMLHalfCells is defined as PML, and these
     * cells are intersected with mPhysicalYeeCells.
    **/
    Rect3i pmlHalfCells(int face) const;
    
    void physicalYeeCells(const Rect3i & r);
    void nonPMLHalfCells(const Rect3i & r);
    void physicalBoundaries(const std::vector<BoundaryType> & bounds);

private:
    Rect3i mPhysicalYeeCells;
    Rect3i mNonPMLHalfCells;
    std::vector<BoundaryType> mBoundaries;
    int mTransverseAxis;
};
bool operator==(const PhysicalExtents & lhs, const PhysicalExtents & rhs);

class NodeExtents
{
public:
    NodeExtents();
    NodeExtents(Rect3i physicalYeeCells, Rect3i calcHalfCells,
        Rect3i ghostHalfCells, Rect3i nonPMLHalfCells);
    
//    NodeExtents & operator=(const NodeExtents & rhs);
    
    /**
     * gridExtents  extents of the entire simulation grid
     * physicalBoundaries   boundary types at edge of calculation
     * nodeYeeCells intersection of grid yee cells with MPI partition edges
    **/
    NodeExtents(const PhysicalExtents & gridExtents, Rect3i nodeYeeCells,
        Vector3i nodeCoordinates = Vector3i(0,0,0),
        Vector3i numNodes = Vector3i(1,1,1));
    
    const Rect3i & physicalYeeCells() const { return mPhysicalYeeCells; }
    const Rect3i & nodeYeeCells() const { return mNodeYeeCells; }
    const Rect3i & calcHalfCells() const { return mCalcHalfCells; }
    const Rect3i & ghostHalfCells() const { return mGhostHalfCells; }
    bool isEntirelyPML() const;
    const Rect3i & nonPMLHalfCells() const { return mNonPMLHalfCells; }
//    const BoundaryType* nodeBoundaries() const { return mNodeBoundaries; }
    const std::vector<BoundaryType> nodeBoundaries() const
        { return mNodeBoundaries; }
    
    /**
     * Returns true for all directions that don't have translation-symmetric
     * boundaries.
     */
    Vector3b dimensions() const;
    
    Rect3i allocatedYeeCells() const
    {
        return YeeUtilities::halfToYee(mGhostHalfCells);
    }
    // TODO: Write unit tests for constantsHalfCells().
    Rect3i constantsHalfCells() const // region of user-defined constants
    {
        return intersection(YeeUtilities::yeeToHalf(physicalYeeCells()),
            ghostHalfCells());
    }
    
    // TODO: test me
    int numGhostHalfCells(int side) const;
    // TODO: test me
    Rect3i ghostReadHalfCells(int readSide) const;
    // TODO: test me
    Rect3i ghostWriteHalfCells(int writeSide) const;
    
    void physicalYeeCells(const Rect3i & r);
    void nodeYeeCells(const Rect3i & r);
    void calcHalfCells(const Rect3i & r);
    void ghostHalfCells(const Rect3i & r);
    void nonPMLHalfCells(const Rect3i & r);
    void nodeBoundaries(const std::vector<BoundaryType> & bounds);
    
private:
    Rect3i mPhysicalYeeCells; // all the cells of the whole simulation
    Rect3i mNodeYeeCells; // all the cells nominally belonging to this node
    Rect3i mCalcHalfCells; // all the cells where update equations are run
    Rect3i mGhostHalfCells; // all the cells where fields need to be defined
    Rect3i mNonPMLHalfCells; // hurr-durr
    
    //BoundaryType mNodeBoundaries[6];
    std::vector<BoundaryType> mNodeBoundaries;
};
//bool operator==(const NodeExtents & lhs, const NodeExtents & rhs);

void getNodeBoundaries(const std::vector<BoundaryType> & physicalBoundaries,
    std::vector<BoundaryType> & outNodeBoundaries, Vector3i nodeCoordinates,
    Vector3i numNodes);



#endif
