/*
 *  VoxelizePolyhedra.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/16/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "VoxelizePolyhedra.h"

#include "Voxelize/Voxelize.h"
#include "Extents.h"
#include "YeeUtilities.h"
#include "Log.h"

#include <iostream>

using namespace std;
using namespace RLE;
using namespace YeeUtilities;


// ======= Static methods

Vector3d sDivideVec(const Vector3d & v, double d)
{
    return v / d;
}

map<unsigned int, RLE::DynamicRLE3<Vector3d> > sDivide(
    const map<unsigned int, RLE::DynamicRLE3<Vector3d> > & rles,
    double val)
{
    map<unsigned int, RLE::DynamicRLE3<Vector3d> > outMap;
    
    map<unsigned int, RLE::DynamicRLE3<Vector3d> >::const_iterator itr;
    
    for (itr = rles.begin(); itr != rles.end(); itr++)
    {
        RLE::scalarTransform(itr->second, val, outMap[itr->first], sDivideVec);
    }
    
    return outMap;
}


// ======= Class methods

VoxelizePolyhedra::
VoxelizePolyhedra(const Rect3i & voxelBounds, const Rect3d & realBounds,
    Vector3b nonSymmetricDimensions, Rect3d clipBounds) :
    mVoxelBounds(voxelBounds),
    mRealBounds(realBounds),
    mClipBounds(clipBounds),
    mNonSymmetricDimensions(nonSymmetricDimensions)
{
}


void VoxelizePolyhedra::
fillFactors(
    const vector<vector<SimpleMesh::Triangle> > & polyhedra,
    vector<RLE::DynamicRLE3<double> > & outFillFactors) const
{
    outFillFactors.resize(polyhedra.size());
    for (int nn = 0; nn < polyhedra.size(); nn++)
    {
        fillFactors(polyhedra.at(nn), outFillFactors[nn]);
    }
}

void VoxelizePolyhedra::
fillFactors(
    const std::vector<SimpleMesh::Triangle> & polyhedron,
    RLE::DynamicRLE3<double> & outFillFactors) const
{
    Voxelize vox(realBounds(), voxelBounds(), nonSymmetricDimensions());
    
    Vector3d vSize = vox.voxelSize();
    double voxelVolume = vSize[0]*vSize[1]*vSize[2];
    
    outFillFactors = vox.fillFactors(polyhedron) * (1.0/voxelVolume);
    
//    checkFillFactorsBounded(outFillFactors);
}


void VoxelizePolyhedra::
orientations(
    const vector<std::vector<SimpleMesh::Triangle> > & polyhedra,
    RLE::DynamicRLE3<Matrix3d> & outProjectionMatrices) const
{
    Voxelize vox(realBounds(), voxelBounds(), nonSymmetricDimensions());
    
    for (int nn = 0; nn < polyhedra.size(); nn++)
    {
        if (outProjectionMatrices.numRuns() == 0)
        {
            outProjectionMatrices = vox.orientations(
                polyhedra.at(nn), clipBounds());
        }
        else
        {
            outProjectionMatrices = outProjectionMatrices +
                vox.orientations(polyhedra.at(nn), clipBounds());
        }
    }
}



void VoxelizePolyhedra::
fillFactorJacobians(
    const vector<std::vector<SimpleMesh::Triangle> > & polyhedra,
    vector<map<unsigned int, RLE::DynamicRLE3<Vector3d> > > &
        outJacobians) const
{
//    cerr << "Fill factor Jacobians unavailable\n";
    Voxelize vox(realBounds(), voxelBounds(), nonSymmetricDimensions());
    
    Vector3d vSize = vox.voxelSize();
    double voxelVolume = vSize[0]*vSize[1]*vSize[2];
    
    outJacobians.resize(polyhedra.size());
    for (int nn = 0; nn < polyhedra.size(); nn++)
    {
        outJacobians[nn] = sDivide(vox.fillFactorJacobians(polyhedra.at(nn)),
            voxelVolume);
    }
}

void VoxelizePolyhedra::
orientationJacobians(
    const vector<std::vector<SimpleMesh::Triangle> > & polyhedra,
    map<unsigned int, RLE::DynamicRLE3<Matrix3<Vector3d> > > &
        outJacobians) const
{
    Voxelize vox(realBounds(), voxelBounds(), nonSymmetricDimensions());
    
    for (int nn = 0; nn < polyhedra.size(); nn++)
    {
        map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > > onePolyJacobian =
            vox.orientationJacobians(polyhedra.at(nn), clipBounds());
        
        map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > >::const_iterator itr;
        for (itr = onePolyJacobian.begin(); itr != onePolyJacobian.end(); itr++)
        {
            assert(hasCorrectDimensions(itr->second));
            if (outJacobians[itr->first].numRuns() == 0)
            {
                outJacobians[itr->first] = itr->second;
            }
            else
            {
                outJacobians[itr->first] = outJacobians[itr->first] + itr->second;
            }
        }
    }
}


Rect3d VoxelizePolyhedra::
clipBounds() const
{
    if (hasClipBounds())
        return mClipBounds;
    else
        return Rect3d(-1e100, -1e100, -1e100, 1e100, 1e100, 1e100);
}
