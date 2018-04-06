/*
 *  KDTree.cpp
 *  Trogdor6
 *
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 *
 */

#include "KDTree.h"
#include <cassert>

namespace KDTree
{

KDTree::
KDTree(std::vector<KDRect*> rects, int maxDepth) :
    mSize(rects.size()),
    mRoot(0L)
{
    if (rects.size() > 0)
    {
        mRoot = newNode(rects, 0, maxDepth);
        mIntersectedRects.resize(mSize);
    }
}

KDTree::
~KDTree()
{
    if (mRoot != 0L)
        deleteNode(mRoot);
}

void KDTree::
deleteNode(Node* n)
{
    if (n->mLeft)
        deleteNode(n->mLeft);
    if (n->mRight)
        deleteNode(n->mRight);
    delete n;
}

Node* KDTree::
newNode(std::vector<KDRect*> rects, int depth, int maxDepth)
{
    assert(rects.size() > 0);
    Node* node = new Node;
    if (depth+1 >= maxDepth || rects.size() == 1)
    {
        node->mKDRects = rects;
    }
    else
    {
        int medianIndex = rects.size()/2;
        node->mSortCoordinate = depth%2;
        std::sort(rects.begin(), rects.end(), Comparator(depth%4));
        node->mSplitValue = rects[medianIndex]->x[depth%4];
        
        std::vector<KDRect*> leftKDRects;
        std::vector<KDRect*> rightKDRects;
        
        for (int nn = 0; nn < rects.size(); nn++)
        {
            if (rects[nn]->x[node->mSortCoordinate] < node->mSplitValue)
                leftKDRects.push_back(rects[nn]);
            if (rects[nn]->x[node->mSortCoordinate+2] >= node->mSplitValue)
                rightKDRects.push_back(rects[nn]);
        }
        
        if (leftKDRects.size() > 0)
            node->mLeft = newNode(leftKDRects, depth+1, maxDepth);
        if (rightKDRects.size() > 0)
            node->mRight = newNode(rightKDRects, depth+1, maxDepth);
    }
    
    return node;
}


const std::vector<KDRect*> & KDTree::
intersectingKDRects(double x, double y, unsigned int & outNumHits)
    const
{
    Node* curNode = mRoot;
    
    while (curNode != 0L && (curNode->mLeft != 0L || curNode->mRight != 0L))
    {
        if (curNode->mSortCoordinate == 0)
        {
            if (x < curNode->mSplitValue)
                curNode = curNode->mLeft;
            else
                curNode = curNode->mRight;
        }
        else
        {
            if (y < curNode->mSplitValue)
                curNode = curNode->mLeft;
            else
                curNode = curNode->mRight;
        }
    }
    
    outNumHits = 0;
    if (curNode != 0L)
    {
        for (int mm = 0; mm < curNode->mKDRects.size(); mm++)
        {
            if (x > curNode->mKDRects[mm]->x[0] &&
                x < curNode->mKDRects[mm]->x[2] &&
                y > curNode->mKDRects[mm]->x[1] &&
                y < curNode->mKDRects[mm]->x[3])
            {
                mIntersectedRects[outNumHits++] = curNode->mKDRects[mm];
            }
        }
    }
    return mIntersectedRects;
}

}; // namespace KDTree


