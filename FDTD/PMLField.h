/*
 *  PMLField.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 11/11/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _PMLFIELD_
#define _PMLFIELD_

#include "rle/DynamicRLE3.h"
#include "Precision.h"
#include "PMLParameters.h"
#include "Extents.h"
#include "YeeUtilities.h"
#include "Map.h"

class PMLField
{
public:
    PMLField(YeeUtilities::Octant octant);
    
    void
    init(const PhysicalExtents & extents, Precision::Vec3 dxyz,
        const Map<std::string, std::string> & pmlParams);
    
//    void dimensions(Rect3i yeeCells);
    
    /**
     * Returns a three-element array, e.g. constants for Ex, Ey and Ez.
    **/
    RLE::DynamicRLE3<PMLParameters>* along(int xyz)
    {
        return parameters[xyz];
    }
    
    const RLE::DynamicRLE3<PMLParameters>* along(int xyz) const
    {
        return parameters[xyz];
    }
    
    const YeeUtilities::Octant & octant() const { return mOctant; }
    
    double bytes() const
    {
        double total = 0.0;
        for (int mm = 0; mm < 3; mm++)
        for (int nn = 0; nn < 3; nn++)
            total += parameters[mm][nn].bytes();
        return total;
    }
    
private:
    YeeUtilities::Octant mOctant;
    RLE::DynamicRLE3<PMLParameters> parameters[3][3];
};

/**
 * Return an array containing 0.0 in non-PML regions and values from 0.0 to 1.0
 * representing depth into the PML in other places.  The arrays will be 1D,
 * along the X, Y or Z directions.
**/
RLE::DynamicRLE3<double> pmlDepth(const PhysicalExtents & extents, int xyz);



#endif
