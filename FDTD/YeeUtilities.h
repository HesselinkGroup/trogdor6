/*
 *  YeeUtilities.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 2/13/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _YEEUTILITIES_
#define _YEEUTILITIES_

#include "PrecisionGeometry.h"

// TO-DO: call things "octants" as applicable

namespace YeeUtilities
{

// returns a number from [0,7] telling which octant of the Yee cell v is in
int octant(const Vector3i & v);

// returns a vector with origin at 0 pointing to the given Yee octant
const Vector3i & halfCellOffset(int octant);

const Vector3d & halfCellPosition(int octant);

const Vector3i & cardinal(int directionIndex);

int xyz(int octant);
bool isE(int octant);
bool isH(int octant);

int octantE(int directionIndex); // returns octant 1 for dir 0 (Ex), etc.
int octantH(int directionIndex); // returns octant 6 for dir 0 (Hx), etc.

int octantE(int ii, int jj); // returns 0 when ii != jj, octantE(ii) otherwise
int octantH(int ii, int jj); // returns 7 when ii != jj, octantH(ii) otherwise

// Functor version of octantE/octantH, for tensorial permittivity/permeability
class Octant
{
public:
    Octant() : mIsElectric(true) {}
    static Octant electric();
    static Octant magnetic();
    
    int operator()(int i, int j) const;
    
    bool isElectric() const { return mIsElectric; }
    bool isMagnetic() const { return !mIsElectric; }
    
    bool operator==(const Octant & rhs) const
    {
        return mIsElectric == rhs.mIsElectric;
    }
    
private:
    bool mIsElectric; // else, magnetic.
};

Precision::Vec3 eFieldPosition(int fieldIndex);
Precision::Vec3 hFieldPosition(int fieldIndex);
Vector3i eFieldOffset(int directionIndex);
Vector3i hFieldOffset(int directionIndex);

// returns halfCell/2 for strictly positive cell indices.
Vector3i halfToYee(const Vector3i & halfCell);

// returns  2*yeeCell + halfCellOffset
Vector3i yeeToHalf(const Vector3i & yeeCell,
	const Vector3i & halfCellOffset);

// returns 2*yeeCell + halfCellOffset(halfCellIndex)
Vector3i yeeToHalf(const Vector3i & yeeCell, int halfCellIndex);

// returns smallest Yee rect containing all points in halfRect
Rect3i halfToYee(const Rect3i & halfRect);

// returns smallest Yee rect containing all points at given halfCellIndex
// in halfRect
Rect3i halfToYee(const Rect3i & halfRect, int octant);

// returns smallest Yee rect containing all points at given halfCellOffset
// in halfRect
Rect3i halfToYee(const Rect3i & halfRect, const Vector3i & halfCellOffset);

// returns smallest half cell rect containing all points in given Yee rect
// at given halfCellIndex
Rect3i yeeToHalf(const Rect3i & yeeRect, int octant);

// returns smallest half cell rect containing all points in given Yee rect
// at given half cell offset
Rect3i yeeToHalf(const Rect3i & yeeRect, const Vector3i & halfCellOffset);

// returns smallest half cell rect containing all points in given Yee rect
Rect3i yeeToHalf(const Rect3i & yeeRect);

// equivalent to yeeToHalf(rectHalftoYee(halfRect))
Rect3i expandToYeeRect(const Rect3i & halfRect);

// returns the last layer of the rect on the side with normal along sideIndex
Rect3i edgeOfRect(const Rect3i & rect, int sideIndex);

// returns a rectangle that spans the given number of Yee cells, with each
// cell centered at the correct octant position.  This is the right way to ask
// what region of physical space is included in a block of Yee cells.
Rect3d realBounds(const Rect3i & yeeRect, int octant);

}; // namespace YeeUtilities

#endif
