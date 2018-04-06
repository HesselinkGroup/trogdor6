/*
 *  Raycaster.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 2/24/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "Raycaster.h"

#include "TimeWrapper.h"

#include <fstream>

using namespace std;
using namespace RLE;

Raycaster::
Raycaster(int rayDirectionXYZ, Vector2d rayHalfWidth,
    const std::vector<const SimpleMesh::Triangle*> & triangles) :
    mxyz(rayDirectionXYZ)
{
    for (unsigned int tt = 0; tt < triangles.size(); tt++)
    {
        mTriangleBoxes.push_back(TriangleInBox(triangles.at(tt),
            rayDirectionXYZ, rayHalfWidth));
        
//        for (int nn = 0; nn < 4; nn++)
//            cerr << mTriangleBoxes[tt].x[nn] << " ";
//        cerr << "\n";
    }
    
    mRectPointers.resize(mTriangleBoxes.size());
    for (unsigned int nn = 0; nn < mTriangleBoxes.size(); nn++)
    {
        mRectPointers[nn] = &(mTriangleBoxes[nn]);
//        cerr << "Tri " << triangles[nn]->triangle() << " "
//            << " rect " << nn << ": " << mTriangleBoxes[nn].x[0]
//            << " " << mTriangleBoxes[nn].x[1] << " "
//            << mTriangleBoxes[nn].x[2] << " "
//            << mTriangleBoxes[nn].x[3] << "\n";
    }
    
    mKDTree = new KDTree::KDTree(mRectPointers, kKDTREE_MAX_DEPTH);
}

Raycaster::
~Raycaster()
{
    delete mKDTree;
}

int Raycaster::
boxHits(const Vector2d & rayPt, vector<RayIntersection> &
    rayHitBuffer)
{
    unsigned int numRectHits;
    const vector<KDTree::KDRect*> & intersectedRects =
        mKDTree->intersectingKDRects(rayPt[0], rayPt[1], numRectHits);
    
    if (rayHitBuffer.size() < numRectHits)
        throw(std::logic_error("Ray hit buffer is not large enough"));
    
//    cerr << "Cast pt " << rayPt << "\n";
    
    for (int nn = 0; nn < numRectHits; nn++)
    {
        TriangleInBox* box = (TriangleInBox*)intersectedRects[nn];
        Vector3d interPt = intersection(
            box->triangle()->triangle().plane(), vec3(rayPt, 0.0),
            rayDirection());
        rayHitBuffer[nn] = RayIntersection(box->triangle(), interPt);
        
//        cerr << "Tri " << box->triangle()->triangle() << "\n";
    }
    sort(rayHitBuffer.begin(), rayHitBuffer.begin() + numRectHits,
        RayIntersection::SortAlong(uRay()));
    
    return numRectHits;
}


