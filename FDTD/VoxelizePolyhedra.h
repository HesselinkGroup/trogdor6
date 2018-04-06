/*
 *  VoxelizePolyhedra.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/16/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _VOXELIZEPOLYHEDRA_
#define _VOXELIZEPOLYHEDRA_

#include "rle/DynamicRLE3.h"
#include "rle/SupportRegion3.h"
#include "SimpleMesh/SimpleMesh.h"
#include "PrecisionRationalFunction.h"
#include "Pointer.h"
#include "geometry.h"

#include <map>
#include <iostream>

class NodeExtents;

class VoxelizePolyhedra
{
public:
    VoxelizePolyhedra(const Rect3i & voxelBounds, const Rect3d & realBounds,
        Vector3b nonSymmetricDimensions = Vector3b(true, true, true),
        Rect3d clipBounds = Rect3d(0,0,0,0,0,0));
    
    void setLog(std::ostream & str);
    void clearLog();
    
    void fillFactors(
        const std::vector<std::vector<SimpleMesh::Triangle> > & polyhedra,
        std::vector<RLE::DynamicRLE3<double> > & outFillFactors) const;
    
    void fillFactors(
        const std::vector<SimpleMesh::Triangle> & polyhedron,
        RLE::DynamicRLE3<double> & outFillFactors) const;
    
    void orientations(
        const std::vector<std::vector<SimpleMesh::Triangle> > & polyhedra,
        RLE::DynamicRLE3<Matrix3d> & outProjectionMatrices) const;
    
    void fillFactorJacobians(
        const std::vector<std::vector<SimpleMesh::Triangle> > & polyhedra,
        std::vector<std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > > &
            outJacobians) const;
    
    void orientationJacobians(
        const std::vector<std::vector<SimpleMesh::Triangle> > & polyhedra,
        std::map<unsigned int, RLE::DynamicRLE3<Matrix3<Vector3d> > > &
            outJacobians) const;
    
    const Rect3i & voxelBounds() const { return mVoxelBounds; }
    const Rect3d & realBounds() const { return mRealBounds; }
    Rect3d clipBounds() const;
    bool hasClipBounds() const { return mClipBounds != Rect3d(0,0,0,0,0,0); }
    const Vector3b & nonSymmetricDimensions() const { return mNonSymmetricDimensions; }
    
private:
    template<class T>
    bool hasCorrectDimensions(const RLE::RLEBase3<T> & rle) const
    {
        for (int ii = 0; ii < 3; ii++)
        {
            if (rle.hasDimension(ii) != mNonSymmetricDimensions[ii])
            {
                return false;
            }
        }
        return true;
    }
    Rect3i mVoxelBounds;
    Rect3d mRealBounds;
    Rect3d mClipBounds;
    Vector3b mNonSymmetricDimensions;
};










#endif
