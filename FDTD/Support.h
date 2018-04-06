/*
 *  Support.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/30/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _SUPPORT_
#define _SUPPORT_

#include "rle/SupportRegion3.h"
#include "SimulationDescription.h"
#include <vector>

//RLE::SupportRegion3 support(const OutputDescription & desc,
//    Vector3i bigDimensions, int octant);
RLE::SupportRegion3 support(const SourceDescription & desc,
    Vector3i bigDimensions, int octant);
RLE::SupportRegion3 support(const CurrentSourceDescription & desc,
    Vector3i bigDimensions, int octant);
RLE::SupportRegion3 support(const HuygensSurfaceDescription & desc,
    Vector3i bigDimensions, int octant);
RLE::SupportRegion3 support(const HuygensSurfaceDescription & desc,
    int face, Vector3i bigDimensions, int octant);

RLE::SupportRegion3 support(const std::vector<Region> & regions,
    Vector3i bigDimensions = Vector3i(10,10,10));

RLE::SupportRegion3 support(const Rect3i & r,
    Vector3b dimensions = Vector3b(1,1,1));


#endif
