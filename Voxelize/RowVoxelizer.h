/*
 *  RowVoxelizer.h
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 7/5/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef ROWVOXELIZER_H
#define ROWVOXELIZER_H

#include "Raycaster.h"
#include "rle/DynamicRLE3.h"
#include "ClipTriangle.h"
#include <map>


struct FillRecord
{
    FillRecord()
    {}
    
    FillRecord(int inVoxel, double inArea, double inVolume) :
        voxel(inVoxel),
        area(inArea),
        volume(inVolume)
    {}
        
    int voxel;
    double area;
    double volume;
};
inline bool operator<(const FillRecord & lhs, const FillRecord & rhs)
{
    return lhs.voxel < rhs.voxel;
}
inline std::ostream & operator<<(std::ostream & str, const FillRecord & fr)
{
    str << "vox " << fr.voxel << " A " << fr.area << " V " << fr.volume;
    return str;
}

class RowVoxelizer
{
public:
    RowVoxelizer(int rowAxis, Rect3d bounds, Rect3i voxelBounds, Vector3b nonSymmetricDimensions) :
        mRowAxis(rowAxis),
        mBounds(bounds),
        mVoxelBounds(voxelBounds),
        mVoxelSize(bounds.size() / Vector3d(voxelBounds.num())),
        mNonSymmetricDimensions(nonSymmetricDimensions)
    {
    }
    
    int rowAxis() const { return mRowAxis; }
    const Rect3d & bounds() const { return mBounds; }
    const Rect3i & voxelBounds() const { return mVoxelBounds; }
    const Vector3d & voxelSize() const { return mVoxelSize; }
    const Vector3b & nonSymmetricDimensions() const { return mNonSymmetricDimensions; }
    
    int voxel(double height) const
    {
        return voxelBounds().p1[rowAxis()] +
            floor( (height - bounds().p1[rowAxis()]) / voxelSize()[rowAxis()] );
    }
    
    double voxelMin(int voxel) const
    {
        return bounds().p1[rowAxis()] +
            (voxel - voxelBounds().p1[rowAxis()])*voxelSize()[rowAxis()];
    }
    
    Rect2d pixelBounds(int jj, int kk) const
    {
        Vector3d vs = voxelSize();
        int u1 = (mRowAxis+1)%3;
        int u2 = (u1+1)%3;
        
        int jj_relative = jj - voxelBounds().p1[u1],
            kk_relative = kk - voxelBounds().p1[u2];
        
        return Rect2d( bounds().p1[u1]+jj_relative*vs[u1],
            bounds().p1[u2]+kk_relative*vs[u2],
            bounds().p1[u1]+(jj_relative+1)*vs[u1],
            bounds().p1[u2]+(kk_relative+1)*vs[u2] );
    }
    
    // clip without calculating all sensitivities of clipped edges
    ClippedTriangle clipToRow(const SimpleMesh::Triangle & tri, int jj, int kk)
        const;
    
    // clip and include sensitivities of clipped edges
    ClippedTriangle clipToRow(const SimpleMesh::Triangle & tri, int jj, int kk,
        std::vector<SimpleMesh::SensitiveEdge> & edgeBuffer)
        const;
    
    virtual void voxelize(int jj, int kk,
        const std::vector<RayIntersection> & intersections, int numHits) {}
    
private:
    int mRowAxis;
    Rect3d mBounds;
    Rect3i mVoxelBounds;
    Vector3d mVoxelSize;
    Vector3b mNonSymmetricDimensions;
};


class FillFactorRowVoxelizer : public RowVoxelizer
{
public:
    FillFactorRowVoxelizer(int rowAxis, Rect3d bounds, Rect3i voxelBounds, Vector3b nonSymmetricDimensions,
        RLE::DynamicRLE3<double> & outFillFactors);
    
    virtual void voxelize(int jj, int kk,
        const std::vector<RayIntersection> & intersections, int numHits);
    
    void voxelizeOneTri(const SimpleMesh::Triangle* tri, int jj, int kk,
        std::vector<FillRecord> & fillRecords);
    
    void handleSortedFillRecords(int jj, int kk,
        const std::vector<FillRecord> & fillRecords);
    
    static double fillThreshold() { return 1e-10; }
    
    void fill(int cell, int jj, int kk, double val)
        { fill(cell, cell, jj, kk, val); }
    
    void fill(int cell0, int cell1, int jj, int kk, double val)
    {
        assert(rowAxis() == 0);
        
        if (cell0 < voxelBounds().p1[rowAxis()])
            cell0 = voxelBounds().p1[rowAxis()];
        if (cell1 > voxelBounds().p2[rowAxis()])
            cell1 = voxelBounds().p2[rowAxis()];
        
        if (cell0 > cell1) return;
        if (fabs(val) > fillThreshold())
            mFillFactors.mark(cell0, jj, kk, cell1, jj, kk, val);
        
//        if (val > 1.00001*voxelSize()[0]*voxelSize()[1]*voxelSize()[2])
//            std::cerr << "Got one here.\n";
    }
    
private:
    RLE::DynamicRLE3<double> & mFillFactors;
};

struct OrientationRecord
{
    OrientationRecord()
    {}
    
    OrientationRecord(int inVoxel, const ClippedTriangle* clipped) :
        voxel(inVoxel),
        polygon(clipped)
    {}
        
    int voxel;
    const ClippedTriangle* polygon;
//    Matrix3d orientationTimesArea;
};
inline bool operator<(const OrientationRecord & lhs, const OrientationRecord & rhs)
{
    return lhs.voxel < rhs.voxel;
}
inline std::ostream & operator<<(std::ostream & str, const OrientationRecord & ori)
{
    str << "vox " << ori.voxel;
    return str;
}

class OrientationRowVoxelizer : public RowVoxelizer
{
public:
    OrientationRowVoxelizer(int rowAxis, Rect3d bounds, Rect3i voxelBounds, Vector3b nonSymmetricDimensions,
        RLE::DynamicRLE3<Matrix3d> & outOrientations);
    
    virtual void voxelize(int jj, int kk,
        const std::vector<RayIntersection> & intersections, int numHits);
    
    void voxelizeOneTri(const ClippedTriangle* tri, int jj, int kk,
        std::vector<OrientationRecord> & orientationRecords);
    
    void handleSortedOrientationRecords(int jj, int kk,
        const std::vector<OrientationRecord> & records);
    
    void fill(int cell, int jj, int kk, Matrix3d val)
        { fill(cell, cell, jj, kk, val); }
    
    void fill(int cell0, int cell1, int jj, int kk, const Matrix3d & val)
    {
        assert(rowAxis() == 0);
        if (cell0 < voxelBounds().p1[rowAxis()])
            cell0 = voxelBounds().p1[rowAxis()];
        if (cell1 > voxelBounds().p2[rowAxis()])
            cell1 = voxelBounds().p2[rowAxis()];
        if (cell0 > cell1) return;
        mOrientations.mark(cell0, jj, kk, cell1, jj, kk, val);
    }
    
private:
    RLE::DynamicRLE3<Matrix3d> & mOrientations;
};

struct DFillRecord
{
    DFillRecord() {}
    DFillRecord(int inVoxel, const ClippedTriangle* clipped) :
        voxel(inVoxel),
        polygon(clipped)
    {}
    
    int voxel;
    const ClippedTriangle* polygon;
};
inline bool operator<(const DFillRecord & lhs, const DFillRecord & rhs)
{
    return lhs.voxel < rhs.voxel;
}

class DFillFactorRowVoxelizer : public RowVoxelizer
{
public:
    DFillFactorRowVoxelizer(int rowAxis, Rect3d bounds, Rect3i voxelBounds, Vector3b nonSymmetricDimensions,
        std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > & outDFillFactors);
    
    virtual void voxelize(int jj, int kk,
        const std::vector<RayIntersection> & intersections, int numHits);
    
    void voxelizeOneTri(const ClippedTriangle* polygon, int jj, int kk,
        std::vector<DFillRecord> & records);
    
    void handleSortedRecords(int jj, int kk,
        const std::vector<DFillRecord> & records);
    
    void fill(int cell, int jj, int kk,
        const std::map<unsigned int, Vector3d> & val)
    {
        assert(rowAxis() == 0);
        
        if (cell < voxelBounds().p1[rowAxis()] ||
            cell > voxelBounds().p2[rowAxis()])
        {
            return;
        }
        
        std::map<unsigned int, Vector3d>::const_iterator itr;
        for (itr = val.begin(); itr != val.end(); itr++)
        {
            if (mDFillFactors.count(itr->first) == 0)
            {
                mDFillFactors[itr->first] = RLE::DynamicRLE3<Vector3d>(nonSymmetricDimensions().asArray());
            }
            RLE::DynamicRLE3<Vector3d> & rle = mDFillFactors[itr->first];
            rle.mark(cell, jj, kk, rle.markAt(cell, jj, kk, Vector3d(0,0,0)) + itr->second);
        }
    }
    
private:
    std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > & mDFillFactors;
};


struct DOrientationRecord
{
    DOrientationRecord() {}
    DOrientationRecord(int inVoxel, const ClippedTriangle* clipped, int xyz) :
        voxel(inVoxel),
        axis(xyz),
        polygon(clipped)
    {}
    
    int voxel;
    int axis;
    const ClippedTriangle* polygon;
};
inline bool operator<(const DOrientationRecord & lhs, const DOrientationRecord & rhs)
{
    return lhs.voxel < rhs.voxel;
}

class DOrientationRowVoxelizer : public RowVoxelizer
{
public:
    DOrientationRowVoxelizer(int rowAxis, Rect3d bounds, Rect3i voxelBounds, Vector3b nonSymmetricDimensions,
        std::map<unsigned int, RLE::DynamicRLE3<Matrix3<Vector3d> > > & outDOrientation);
    
    virtual void voxelize(int jj, int kk,
        const std::vector<RayIntersection> & intersections, int numHits);

    void voxelizeOneTri(const ClippedTriangle* polygon, int jj, int kk);
    
    void fill(int cell, int jj, int kk,
        const std::map<unsigned int, Matrix3<Vector3d> > & val)
    {
        assert(rowAxis() == 0);
        
        if (cell < voxelBounds().p1[rowAxis()] ||
            cell > voxelBounds().p2[rowAxis()])
        {
            return;
        }
        
        std::map<unsigned int, Matrix3<Vector3d> >::const_iterator itr;
        for (itr = val.begin(); itr != val.end(); itr++)
        {
            if (mDOrientation.count(itr->first) == 0)
            {
                mDOrientation[itr->first] = RLE::DynamicRLE3<Matrix3<Vector3d> >(nonSymmetricDimensions().asArray());
            }
            RLE::DynamicRLE3<Matrix3<Vector3d> > & rle = mDOrientation[itr->first];
            rle.mark(cell, jj, kk, rle.markAt(cell, jj, kk, mZedMatrix) + itr->second);
        }
    }
    
private:
    static Matrix3<Vector3d> mZedMatrix;
    std::map<unsigned int, RLE::DynamicRLE3<Matrix3<Vector3d> > > & mDOrientation;
};










#endif
