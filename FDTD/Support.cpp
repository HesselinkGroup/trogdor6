/*
 *  Support.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/30/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "Support.h"
#include "YeeUtilities.h"

using namespace std;
using namespace YeeUtilities;
using namespace RLE;

//RLE::SupportRegion3 
//support(const OutputDescription & desc, Vector3i dimensions, int octant)
//{
//    if (isE(octant))
//    {
//        if (desc.whichE()[xyz(octant)] ||
//            desc.whichJ()[xyz(octant)] ||
//            desc.whichD()[xyz(octant)])
//        {
//            return support(desc.regions(), dimensions);
//        }
//    }
//    else if (isH(octant))
//    {
//        if (desc.whichH()[xyz(octant)] ||
//            desc.whichB()[xyz(octant)] ||
//            desc.whichM()[xyz(octant)])
//        {
//            return support(desc.regions(), dimensions);
//        }
//    }
//    return RLE::SupportRegion3(dimensions[0] > 1,
//        dimensions[1] > 1,
//        dimensions[2] > 1);
//}

RLE::SupportRegion3 
support(const SourceDescription & desc, Vector3i dimensions, int octant)
{
    if (isE(octant))
    {
        if (desc.sourceFields().whichE()[xyz(octant)])
            return support(desc.regions(), dimensions);
    }
    else if (isH(octant))
    {
        if (desc.sourceFields().whichH()[xyz(octant)])
            return support(desc.regions(), dimensions);
    }
    return RLE::SupportRegion3(dimensions[0] > 1,
        dimensions[1] > 1,
        dimensions[2] > 1);
}

RLE::SupportRegion3 
support(const CurrentSourceDescription & desc, Vector3i dimensions,
    int octant)
{
    if (isE(octant))
    {
        if (desc.sourceCurrents().whichJ()[xyz(octant)] ||
            desc.sourceCurrents().whichJE()[xyz(octant)])
            return support(desc.regions(), dimensions);
    }
    else if (isH(octant))
    {
        if (desc.sourceCurrents().whichM()[xyz(octant)] ||
            desc.sourceCurrents().whichMH()[xyz(octant)])
            return support(desc.regions(), dimensions);
    }
    return RLE::SupportRegion3(dimensions[0] > 1,
        dimensions[1] > 1,
        dimensions[2] > 1);
}

RLE::SupportRegion3 
support(const HuygensSurfaceDescription & desc, Vector3i dimensions,
    int octant)
{
    RLE::SupportRegion3 sr(dimensions[0] > 1,
        dimensions[1] > 1,
        dimensions[2] > 1);
    
    for (unsigned int nSide = 0; nSide < 6; nSide++)
    if (false == desc.omittedSides().count(cardinal(nSide)))
    {
        Rect3i yeeCells = desc.yeeCells(nSide, octant);
        
        sr.mark(yeeCells.p1[0], yeeCells.p1[1], yeeCells.p1[2],
            yeeCells.p2[0], yeeCells.p2[1], yeeCells.p2[2]);
    }
    
    return sr;
}

RLE::SupportRegion3 
support(const HuygensSurfaceDescription & desc, int face, Vector3i dimensions,
    int octant)
{
    RLE::SupportRegion3 sr(dimensions[0] > 1,
        dimensions[1] > 1,
        dimensions[2] > 1);
    
    if (!desc.omitsSide(face))
//    if (false == desc.omittedSides().count(cardinal(face)))
    {
        Rect3i yeeCells = desc.yeeCells(face, octant);
        
        sr.mark(yeeCells.p1[0], yeeCells.p1[1], yeeCells.p1[2],
            yeeCells.p2[0], yeeCells.p2[1], yeeCells.p2[2]);
    }
    
    return sr;
}

RLE::SupportRegion3 
support(const vector<Region> & regions, Vector3i dimensions)
{
    RLE::SupportRegion3 sr(dimensions[0] > 1,   
        dimensions[1] > 1,
        dimensions[2] > 1);
    
    for (unsigned int rr = 0; rr < regions.size(); rr++)
    {
        Rect3i rect = regions[rr].yeeCells();
//        if (regions[rr].stride() != Vector3i(0,0,0))
//            LOG << "Warning: not handling output stride properly yet.\n";
        sr.mark(rect.p1[0], rect.p1[1], rect.p1[2],
            rect.p2[0], rect.p2[1], rect.p2[2]);
    }
    
    return sr;
}

RLE::SupportRegion3
support(const Rect3i & r, Vector3b dimensions)
{
    RLE::SupportRegion3 sr(dimensions[0], dimensions[1], dimensions[2]);
    sr.mark(r.p1[0], r.p1[1], r.p1[2], r.p2[0], r.p2[1], r.p2[2]);
    
    return sr;
}




