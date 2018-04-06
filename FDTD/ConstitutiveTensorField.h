/*
 *  ConstitutiveTensorField.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/23/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _CONSTITUTIVETENSORFIELD_
#define _CONSTITUTIVETENSORFIELD_

#include "rle/DynamicRLE3.h"
#include "rle/SupportRegion3.h"
#include "PrecisionRationalFunction.h"
#include "Extents.h"
#include "YeeUtilities.h"

#include <vector>

class ConstitutiveTensorField
{
public:
    ConstitutiveTensorField();
    ConstitutiveTensorField(const YeeUtilities::Octant & o, Rect3i bounds);
    
    void setDimensions(Rect3i yeeCells);
    
    void extrude(Rect3i innerHalfCells, Rect3i outerHalfCells);
    void copy(const ConstitutiveTensorField & source, Rect3i copyCells,
        Rect3i pasteCells);
    
    /**
     * Select all cells where the permittivity has the given order.  If
     * numeratorOrder or denominatorOrder is -1, it will match any order.
     */
    RLE::SupportRegion3 support(int i, int j, int numeratorOrder,
        int denominatorOrder) const;
    
    RLE::DynamicRLE3<Precision::RationalFunction> filter(int i, int j,
        int numeratorOrder, int denominatorOrder) const;
    
    RLE::DynamicRLE3<Precision::RationalFunction> & operator()(int i, int j);
    const RLE::DynamicRLE3<Precision::RationalFunction> & operator()(int i, int j) const;
    
    int tensorOctant(int i, int j) const; // IMPORTANT FUNCTION
    
    const YeeUtilities::Octant & fieldOctant() const { return mOctant; }
    double bytes() const;
    
    void writeBinary(std::string fileName, Rect3i cells) const;
private:
    YeeUtilities::Octant mOctant;
    Rect3i mBounds;
    
    RLE::DynamicRLE3<Precision::RationalFunction> mComponents[3][3];
};


#endif
