/*
 *  KDTree.h
 *  Trogdor6
 *
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _KDTREE_
#define _KDTREE_

#include <vector>
#include <algorithm>

namespace KDTree
{

struct KDRect
{
    KDRect() {}
    KDRect(double minX0, double maxX0, double minX1, double maxX1)
    {
        x[0] = minX0;
        x[1] = maxX0;
        x[2] = minX1;
        x[3] = maxX1;
    }
    double x[4];
};

struct Comparator
{
    Comparator();
    Comparator(int coordinate) :
        mCoordinate(coordinate)
    {
    }
    
    bool operator()(const KDRect* lhs, const KDRect* rhs)
    {
        return (lhs->x[mCoordinate] < rhs->x[mCoordinate]);
    }
    
    int mCoordinate;
};

struct Node
{
    Node() { mLeft = 0L; mRight = 0L; }
    Node* mLeft;
    Node* mRight;
    double mSplitValue;
    int mSortCoordinate; // x or y (0 or 1)
    
    std::vector<KDRect*> mKDRects;
    
    long long bytes()
    {
        long long leftBytes = 0, rightBytes = 0;
        
        if (mLeft)
            leftBytes = mLeft->bytes();
        if (mRight)
            rightBytes = mRight->bytes();
        return (leftBytes + rightBytes + sizeof(mKDRects)
            + mKDRects.size()*sizeof(KDRect*));
    }
};

class KDTree
{
public:
    KDTree(std::vector<KDRect*> rects, int maxDepth);
    ~KDTree();
    const std::vector<KDRect*> & intersectingKDRects(double x, double y,
        unsigned int & outNumHits) const;
    unsigned int size() const { return mSize; }
    const std::vector<KDRect*> & intersections() const
        { return mIntersectedRects; }
private:
    Node* newNode(std::vector<KDRect*> rects, int depth, int maxDepth);
    void deleteNode(Node* n);
    
    mutable std::vector<KDRect*> mIntersectedRects;
    unsigned int mSize;
    Node* mRoot;
};

}; // namespace KDTree

#endif
