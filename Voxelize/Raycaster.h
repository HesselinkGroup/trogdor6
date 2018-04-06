/*
 *  Raycaster.h
 *  raycast_xml
 *
 *  Created by Paul Hansen on 2/24/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _RAYCASTER_
#define _RAYCASTER_

#include <vector>
#include <iostream>
#include "geometry.h"
#include "geometry2.h"
#include "rle/DynamicRLE3.h"
#include "KDTree.h"
#include "SimpleMesh/SimpleMesh.h"

class RayIntersection
{
public:
    RayIntersection() {}
    RayIntersection(const SimpleMesh::Triangle* facet, Vector3d location) :
        mHitTriangle(facet),
        mLocation(location)
    {
    }
    
    const SimpleMesh::Triangle* triangle() const { return mHitTriangle; }
    const Vector3d & location() const { return mLocation; }
    
    struct SortAlong
    {
    public:
        SortAlong() : mAxis(0) {}
        SortAlong(unsigned int axis) : mAxis(axis) {}
        
        bool operator()(const RayIntersection & lhs,
            const RayIntersection & rhs) const
        {
            return lhs.location()[mAxis] < rhs.location()[mAxis];
        }
        
        unsigned int mAxis;
    };
    
protected:
    const SimpleMesh::Triangle* mHitTriangle;
    Vector3d mLocation;
};

struct TriangleInBox : public KDTree::KDRect
{
    TriangleInBox() :
        mRaycastTriangle(0L)
    {
    }
    
    TriangleInBox(const SimpleMesh::Triangle * facet, int rayAxis,
        const Vector2d & radius) :
        mRaycastTriangle(facet)
    {
        int uRay = rayAxis;
        int u0 = (uRay+1)%3;
        int u1 = (u0+1)%3;
        Rect3d bounds = facet->triangle().bounds();
        x[0] = bounds.p1[u0] - radius[0];
        x[1] = bounds.p1[u1] - radius[1];
        x[2] = bounds.p2[u0] + radius[0];
        x[3] = bounds.p2[u1] + radius[1];
    }
    
    const SimpleMesh::Triangle* triangle() const
        { return mRaycastTriangle; }
private:
    const SimpleMesh::Triangle* mRaycastTriangle;
};

class Raycaster
{
public:
    Raycaster(int rayDirectionXYZ, Vector2d rayHalfWidth,
        const std::vector<const SimpleMesh::Triangle*> & triangles);
    ~Raycaster();
    
    int boxHits(const Vector2d & rayPt,
        std::vector<RayIntersection> & rayHitBuffer);
    
    Vector3d rayDirection() const { return Vector3d::unit(mxyz); }
    Vector3d vec3(const Vector2d & transverse, double normal) const
    {
        Vector3d v;
        v[u0()] = transverse[0];
        v[u1()] = transverse[1];
        v[uRay()] = normal;
        return v;
    }
    int u0() const { return (mxyz+1)%3; }
    int u1() const { return (mxyz+2)%3; }
    int uRay() const { return mxyz; }
    
protected:
    enum { kKDTREE_MAX_DEPTH = 9 };
//    std::vector<const SimpleMesh::Triangle*> mHitTriangles;
//    std::vector<RayIntersection> mRayHitBuffer;
    
    //std::vector<RaycastTriangle> mRaycastTriangles;
    std::vector<TriangleInBox> mTriangleBoxes;
    std::vector<KDTree::KDRect*> mRectPointers;
    KDTree::KDTree* mKDTree;
    //int mRaysPerVoxel;
    int mxyz;
};


#endif
