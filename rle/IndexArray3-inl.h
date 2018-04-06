/*
 *  IndexArray3-inl.h
 *  rle_index
 *
 *  Created by Paul Hansen on 2/2/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _INDEXARRAY3_

#include "SupportRegion3.h"

namespace RLE
{

IndexArray3::
IndexArray3() :
    RLEBase3<IASegment>()
{
}

template <class T>
IndexArray3::
IndexArray3(const DynamicRLE3<T> & rle) :
    RLEBase3<IASegment>(rle.nonSymmetricDimensions())
{
    typename DynamicRLE3<T>::ConstIterator itr;
    int64 currentIndex = 0;
    const int64 STRIDE_0 = 0;
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        mark(itr.runStart(), itr.runEnd(), currentIndex++, STRIDE_0);
    }
}

IndexArray3::
IndexArray3(const SupportRegion3 & supportRegion) :
    RLEBase3<IASegment>(supportRegion.nonSymmetricDimensions())
{   
    SupportRegion::ConstIterator itr;
    int64 currentIndex = 0;
    const int64 STRIDE_1 = 1;
    for (itr = supportRegion.begin(); itr != supportRegion.end();
        itr.nextMarkedRun())
    {
        mark(itr.runStart(), itr.runEnd(), currentIndex, STRIDE_1);
        currentIndex += itr.runEnd() - itr.runStart() + 1;
    }
}

IndexArray3::
IndexArray3(const RLEBase3<IASegment> & copyMe) :
    RLEBase3<IASegment>(copyMe)
{
}

IndexArray3 & IndexArray3::
operator=(const RLEBase3<IASegment> & rhs)
{
    ((RLEBase3<IASegment>&)*this) = rhs;
    return *this;
}

std::vector<int64> IndexArray3::
indices() const
{
    std::vector<int64> vec(numRuns());
    
    AVLConstIterator itr = mSegmentTree.begin();
    int64 nn = 0;
    while (itr.hasNext())
        vec[nn++] = itr.next().mark();
    
    return vec;
}

std::vector<uint64> IndexArray3::
lengths() const
{
    std::vector<uint64> vec(numRuns());
    
    AVLConstIterator itr = mSegmentTree.begin();
    int64 nn = 0;
    while (itr.hasNext())
        vec[nn++] = itr.next().length();
    
    return vec;
}

std::vector<uint64> IndexArray3::
strides() const
{
    std::vector<uint64> vec(numRuns());
    
    AVLConstIterator itr = mSegmentTree.begin();
    int64 nn = 0;
    while (itr.hasNext())
        vec[nn++] = itr.next().data().stride();
    
    return vec;
}

std::vector<int64> IndexArray3::
starts() const
{
    std::vector<int64> vec(numRuns());
    
    AVLConstIterator itr = mSegmentTree.begin();
    int64 nn = 0;
    while (itr.hasNext())
        vec[nn++] = itr.next().first();
    
    return vec;
}


}; // namespace RLE

#endif
