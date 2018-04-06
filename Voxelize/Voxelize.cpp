/*
 *  Voxelize.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 4/17/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "Voxelize.h"

#include "TimeWrapper.h"

#include <limits>

using namespace std;
using namespace RLE;
using namespace TimeWrapper;

Voxelize::
Voxelize(const Rect3d & bounds, const Rect3i & voxels, const Vector3b & nonSymmetricDimensions) :
    mBounds(bounds),
    mVoxels(voxels),
    mNonSymmetricDimensions(nonSymmetricDimensions)
{
}

Voxelize::
~Voxelize()
{
}

Vector2d Voxelize::
voxelHalfWidth(int axis) const
{
    Vector3d voxelRadius = voxelSize()/2;
    
    return Vector2d(voxelRadius[(axis+1)%3], voxelRadius[(axis+2)%3]);
}

Rect3d Voxelize::
voxel(int ii, int jj, int kk) const
{
    Vector3d p1 = bounds().p1 +
        Vector3d( (ii-voxelBounds().p1[0])*voxelSize()[0],
            (jj-voxelBounds().p1[1])*voxelSize()[1],
            (kk-voxelBounds().p1[2])*voxelSize()[2]);
    
    Vector3d p2 = p1 + voxelSize();
    
    return Rect3d(p1, p2);
}

DynamicRLE3<double> Voxelize::
fillFactors(const vector<SimpleMesh::Triangle> & triangles) const
{
    int xyz = 0;
    
    DynamicRLE3<double> outFillFactors(mNonSymmetricDimensions.asArray());
    
    vector<const SimpleMesh::Triangle*> triPtrs(triangles.size());
    for (int nn = 0; nn < triangles.size(); nn++)
        triPtrs[nn] = &triangles[nn];
    
    Raycaster rc(xyz, voxelHalfWidth(xyz), triPtrs);
    FillFactorRowVoxelizer vox(xyz, bounds(), voxelBounds(), nonSymmetricDimensions(), outFillFactors);
    raycast(rc, vox);
    
    return outFillFactors;
}



std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > Voxelize::
fillFactorJacobians(const vector<SimpleMesh::Triangle> & triangles) const
{
    int xyz = 0;
    
    map<unsigned int, RLE::DynamicRLE3<Vector3d> > dFdv; // control vert to dFdv
    
    vector<const SimpleMesh::Triangle*> triPtrs;
    for (int nn = 0; nn < triangles.size(); nn++)
    if (triangles[nn].isSensitive())
        triPtrs.push_back(&triangles[nn]);
    
    Raycaster rc(xyz, voxelHalfWidth(xyz), triPtrs);
    DFillFactorRowVoxelizer vox(xyz, bounds(), voxelBounds(), nonSymmetricDimensions(), dFdv);
    raycast(rc, vox);
    
    return dFdv;
}


DynamicRLE3<Matrix3d> Voxelize::
orientations(const vector<SimpleMesh::Triangle> & triangles,
    Rect3d clipBounds) const
{
    int xyz = 0;
    
    DynamicRLE3<Matrix3d> outOrientation(mNonSymmetricDimensions.asArray());
    
//    vector<const SimpleMesh::Triangle*> triPtrs(triangles.size());
//    for (int nn = 0; nn < triangles.size(); nn++)
//        triPtrs[nn] = &triangles[nn];
    
    Rect3d clipInset(inset(clipBounds, 1e-10));
    
    vector<const SimpleMesh::Triangle*> triPtrs;
    for (int nn = 0; nn < triangles.size(); nn++)
    if (clipInset.intersects(triangles[nn].triangle()))
//    if (dominantDirection(triangles[nn].triangle().normal()) != 1)
    {
        triPtrs.push_back(&triangles[nn]);
    }
//    else
//        cerr << "Domcom of " << triangles[nn].triangle().normal() << " is "
//            << dominantDirection(triangles[nn].triangle().normal())
//            << "\n";
    
    Raycaster rc(xyz, voxelHalfWidth(xyz), triPtrs);
    OrientationRowVoxelizer vox(xyz, bounds(), voxelBounds(), nonSymmetricDimensions(), outOrientation);
    raycast(rc, vox);
    
    return outOrientation;
}

std::map<unsigned int, RLE::DynamicRLE3<Matrix3<Vector3d> > > Voxelize::
orientationJacobians(const vector<SimpleMesh::Triangle> & triangles,
    Rect3d clipBounds) const
{
//    cerr << "Calling Jacobian function.\n";
    int xyz = 0;
    
    map<unsigned int, RLE::DynamicRLE3<Matrix3<Vector3d> > > DNNT;
    
    Rect3d clipInset(inset(clipBounds, 1e-10));
    
//    cerr << "Using these triangles:\n";
    vector<const SimpleMesh::Triangle*> triPtrs;
    for (int nn = 0; nn < triangles.size(); nn++)
    if (clipInset.intersects(triangles[nn].triangle()))
    if (triangles[nn].isSensitive())
//    if (dominantDirection(triangles[nn].triangle().normal()) != 1)
    {
//        cerr << "quickPatch(" << triangles[nn].triangle() << ", 'y'); %" << nn << "\n";
        triPtrs.push_back(&triangles[nn]);
    }
    
    Raycaster rc(xyz, voxelHalfWidth(xyz), triPtrs);
    DOrientationRowVoxelizer vox(xyz, bounds(), voxelBounds(), nonSymmetricDimensions(), DNNT);
    raycast(rc, vox);
    
    return DNNT;
    
}


void Voxelize::
raycast(Raycaster & raycaster, RowVoxelizer & rowVoxelizer) const
{
//    cerr << "raycast() in " << bounds() << "\n";
    const int MAX_NUM_INTERSECTIONS = 1000;
    vector<RayIntersection> intersections(MAX_NUM_INTERSECTIONS);
    
//    Vector3i nxyz(voxelBounds().num());
    
    int u0(raycaster.u0()), u1(raycaster.u1()); //, uRay(raycaster.uRay());
    
    for (int kk = voxelBounds().p1[u1]; kk <= voxelBounds().p2[u1]; kk++)
    for (int jj = voxelBounds().p1[u0]; jj <= voxelBounds().p2[u0]; jj++)
    {
        double y0 = bounds().p1[u0] +
            voxelSize()[u0]*
            (jj - voxelBounds().p1[u0] + 0.5);
        double z0 = bounds().p1[u1] +
            voxelSize()[u1]*
            (kk - voxelBounds().p1[u1] + 0.5);
        
//        cout << "Row " << jj << ", " << kk << " cast "
//            << y0 << ", " << z0 << "\n";
        
        int numHits = raycaster.boxHits(Vector2d(y0,z0),
            intersections);
        
//        cout << "Num hits " << numHits << "\n";
        
        if (numHits > 0)
        {
            // DO YER WORK
            rowVoxelizer.voxelize(jj, kk, intersections, numHits);
        }
    }
}





