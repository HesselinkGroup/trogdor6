/*
 *  Voxelize.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 4/17/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 *
 */

#ifndef VOXELIZE_H
#define VOXELIZE_H

#include "geometry.h"
#include "rle/DynamicRLE3.h"
#include "SimpleMesh/SimpleMesh.h"
#include "Raycaster.h"
#include "RowVoxelizer.h"

#include <vector>
#include <map>
#include <iostream>

class Voxelize
{
public:
    Voxelize(const Rect3d & bounds, const Rect3i & voxels, const Vector3b & nonSymmetricDimensions = Vector3b(true,true,true));
    ~Voxelize();
    
    const Rect3d & bounds() const { return mBounds; }
    const Rect3i & voxelBounds() const { return mVoxels; }
    const Vector3d voxelSize() const { return mBounds.size()/mVoxels.num(); }
    Vector2d voxelHalfWidth(int axis) const;
    Vector3b nonSymmetricDimensions() const { return mNonSymmetricDimensions; }
    
    Rect3d voxel(int ii, int jj, int kk) const;
    
    // Values are CLAMPED
    // * smooth here
    RLE::DynamicRLE3<double> fillFactors(
        const std::vector<SimpleMesh::Triangle> & triangles) const;
    
    std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > fillFactorJacobians(
        const std::vector<SimpleMesh::Triangle> & triangles) const;
    
    RLE::DynamicRLE3<Matrix3d> orientations(
        const std::vector<SimpleMesh::Triangle> & triangles,
        Rect3d clipBounds = Rect3d(-1e100, -1e100, -1e100, 1e100, 1e100, 1e100))
        const;
    std::map<unsigned int, RLE::DynamicRLE3<Matrix3<Vector3d> > >
    orientationJacobians(const std::vector<SimpleMesh::Triangle> & triangles,
        Rect3d clipBounds = Rect3d(-1e100, -1e100, -1e100, 1e100, 1e100, 1e100))
        const;
    
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
    
    void raycast(Raycaster & raycaster, RowVoxelizer & vox) const;
    
    Rect3d mBounds;
    Rect3i mVoxels;
    Vector3b mNonSymmetricDimensions;
};




#endif
