/*
 *  PMLField.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 11/11/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "PMLField.h"
#include "PhysicalConstants.h"
#include "RLEOperations.h"
#include "Log.h"

using namespace YeeUtilities;
using namespace RLE;
using namespace std;

PMLField::
PMLField(Octant octant) :
    mOctant(octant)
{
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    {
        parameters[ii][jj] = RLE::DynamicRLE3<PMLParameters>();
    }
}

//void PMLField::
//dimensions(Rect3i yeeCells)
//{
//    bool hasX = (yeeCells.num(0) != 1);
//    bool hasY = (yeeCells.num(1) != 1);
//    bool hasZ = (yeeCells.num(2) != 1);
//    
//    for (int xyz = 0; xyz < 3; xyz++)
//    for (int attenDirection = 0; attenDirection < 3; attenDirection++)
//        along(attenDirection)[xyz].dimensions(hasX, hasY, hasZ);
//}

void PMLField::
init(const PhysicalExtents & extents, Precision::Vec3 dxyz,
    const Map<std::string, std::string> & pmlParams)
{
    // Arrays to store the distance into the PML, indexed by half-cell.
    DynamicRLE3<double> depth[3];
    
    double diagonal = norm(extents.physicalYeeCells().size() * dxyz);
    
    depth[0] = pmlDepth(extents, 0);
    depth[1] = pmlDepth(extents, 1);
    depth[2] = pmlDepth(extents, 2);
    
    for (int attenDir = 0; attenDir < 3; attenDir++)
    for (int xyz = 0; xyz < 3; xyz++)
    {   
        DynamicRLE3<double> depths = downsample2(depth[attenDir],
            yeeToHalf(extents.physicalYeeCells()) +
                halfCellOffset(octant()(xyz, xyz)),
            extents.physicalYeeCells());
        
        CalculatePMLParameters calcPMLParams(dxyz[attenDir], diagonal,
            Constants::eps0, Constants::mu0);
        
        if (pmlParams.count("sigma"))
            calcPMLParams.sigmaEquation(pmlParams["sigma"]);
        if (pmlParams.count("alpha"))
            calcPMLParams.alphaEquation(pmlParams["alpha"]);
        if (pmlParams.count("kappa"))
            calcPMLParams.kappaEquation(pmlParams["kappa"]);
        
        LOGF << "PML " << char('x'+attenDir) << " field " << char('x'+xyz) << ":\n";
        LOGF << "\tsigma = " << calcPMLParams.sigmaEquation() << "\n";
        LOGF << "\tkappa = " << calcPMLParams.kappaEquation() << "\n";
        LOGF << "\talpha = " << calcPMLParams.alphaEquation() << "\n";
        
        transform(depths, along(attenDir)[xyz], calcPMLParams);
        
        if (mOctant.isElectric())
            LOGF << "E";
        else
            LOGF << "H";
        
        LOGFMORE << char('x'+xyz) << " along " << char('x'+attenDir) << ":\n";
        LOGFMORE << along(attenDir)[xyz] << "\n";
//        cerr << "Depth along " << char('x'+attenDir) << " for " << char('x'+xyz)
//            << " is:\n";
//        cerr << depths << "\n";
    }
}


DynamicRLE3<double>
pmlDepth(const PhysicalExtents & extents, int xyz)
{
    // Arrays to store the distance into the PML, indexed by half-cell.
    DynamicRLE3<double> depth(xyz == 0, xyz == 1, xyz == 2);
    
    for (int face = 2*xyz; face <= 2*xyz+1; face++)
    if (extents.hasPML(face))
    {
        Rect3i pmlHalfCells = extents.pmlHalfCells(face);        
//        LOG << "Face " << face << " PML half cells " << pmlHalfCells << "\n";
        
        if (face%2 == 0) // left, bottom, front side
        {
            int pmlDepth = pmlHalfCells.p2[xyz] - pmlHalfCells.p1[xyz] + 1;
            for (int nn = 0; nn < pmlDepth; nn++)
            {
                double d = double(nn+1)/pmlDepth; // d is in (0.0, 1.0]
                Vector3i p(0,0,0);
                p[xyz] = pmlHalfCells.p2[xyz] - nn;
                depth.mark(p[0], p[1], p[2], p[0], p[1], p[2], d);
//                LOG << "Depth " << d << " point " << p << " size "
//                    << depth.numRuns() << "\n";
            }
        }
        else
        {
            int pmlDepth = pmlHalfCells.p2[xyz] - pmlHalfCells.p1[xyz] + 1;
            for (int nn = 0; nn < pmlDepth; nn++)
            {
                double d = double(nn+1)/pmlDepth; // d is in (0.0, 1.0]
                Vector3i p(0,0,0);
                p[xyz] = pmlHalfCells.p1[xyz] + nn;
                depth.mark(p[0], p[1], p[2], p[0], p[1], p[2], d);
//                LOG << "Depth " << d << " point " << p << " size "
//                    << depth.numRuns() << "\n";
            }
        }
    }
//    LOG << "Depth has " << depth.numRuns() << " runs.\n";
    return depth;
}


