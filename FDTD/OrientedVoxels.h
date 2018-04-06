
/*
 *  OrientedVoxels.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/26/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _ORIENTEDVOXELS_
#define _ORIENTEDVOXELS_

#include "Pointer.h"
#include "rle/DynamicRLE3.h"
#include "rle/SupportRegion3.h"
#include "geometry.h"
#include "SimpleMesh/SimpleMesh.h"
#include "YeeUtilities.h"
#include <vector>

// TODO: Split into a creator class and a data containing class?

class OrientedVoxels
{
public:
    OrientedVoxels();
    OrientedVoxels(Rect3i voxelBounds, Rect3d realBounds);
    OrientedVoxels(Rect3i voxelBounds, Rect3d realBounds, Rect3d clipBounds);
    void clipBounds(const Rect3d & r) { mClipBounds = r; }
    void voxelBounds(const Rect3i & r) { mVoxelBounds = r; }
    void realBounds(const Rect3d & r) { mRealBounds = r; }
    const Rect3d & clipBounds() const { return mClipBounds; }
    const Rect3i & voxelBounds() const { return mVoxelBounds; }
    const Rect3d & realBounds() const { return mRealBounds; }
    const Vector3b & nonSymmetricDimensions() const { return mNonSymmetricDimensions; }
    
    void downsampleFillFactors(
        const std::vector<RLE::DynamicRLE3<double> > & halfCellFillFactors,
        int octant);
    void downsampleFillJacobians(
        const std::vector<std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > > &
            halfCellJacobians,
        int octant); // UNTESTED ***
    void initOrientations( // INCOMPLETELY TESTED ***
        const std::vector<std::vector<SimpleMesh::Triangle> > & meshTriangles);
    
    // Remember to filter the Jacobians with a simple predicate to exclude
    // directions that don't move.  (I should do this upstream in Voxelizer.)
    RLE::SupportRegion3 sensitiveCells(int i, int j) const;
    
    const std::vector<RLE::DynamicRLE3<double> > & fillFactors() const; // fast
    const RLE::DynamicRLE3<double> & fillFactors(int material) const; // fast
    const RLE::DynamicRLE3<Matrix3d> & orientations() const; // fast
    RLE::DynamicRLE3<double> orientations(int i, int j) const; // slow
    const std::vector<std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > > &
        fillFactorJacobians() const; // fast
    const std::map<unsigned int, RLE::DynamicRLE3<Matrix3<Vector3d> > > &
        orientationJacobians() const; // fast
    std::map<unsigned int, RLE::DynamicRLE3<Vector3d> >
        orientationJacobians(int i, int j) const; // slow

    void setOrientationSmoothing(int numCells) { mOrientationSmoothing = numCells; }
    void setFillFactorSmoothing(int numCells) { mFillFactorSmoothing = numCells; }
    int orientationSmoothing() const { return mOrientationSmoothing; }
    int fillFactorSmoothing() const { return mFillFactorSmoothing; }
    
    long bytes() const;
private:
    template<class T>
    bool hasCorrectDimensions(const RLE::RLEBase3<T> & rle) const
    {
        for (int ii = 0; ii < 3; ii++)
        {
            if (rle.hasDimension(ii) != mNonSymmetricDimensions[ii])
            {
                std::cout << rle.hasDimension(ii) << " " << mNonSymmetricDimensions[ii] << ".\n";
                return false;
            }
        }
        return true;
    }
    
    Rect3d mRealBounds;
    Rect3i mVoxelBounds;
    Rect3d mClipBounds;
    Vector3b mNonSymmetricDimensions;
    
    std::vector<RLE::DynamicRLE3<double> > mFillFactors;
    RLE::DynamicRLE3<Matrix3d> mOrientation;
    std::vector<std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > >
        mFillFactorJacobians;
    std::map<unsigned int, RLE::DynamicRLE3<Matrix3<Vector3d> > >
        mOrientationJacobians;
    
    int mOrientationSmoothing;
    int mFillFactorSmoothing;
};
















#endif
