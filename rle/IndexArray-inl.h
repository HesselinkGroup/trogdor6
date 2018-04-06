/*
 *  IndexArray-inl.h
 *  refactorRLE
 *
 *  Created by Paul Hansen on 2/8/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _INDEXARRAY_

#include "DynamicRLE.h"
#include "SupportRegion.h"

namespace RLE
{

template<class T>
IndexArray::
IndexArray(const DynamicRLE<T> & rle)
{
    typename DynamicRLE<T>::ConstIterator itr;
    int64 currentIndex = 0;
    const int64 STRIDE_0 = 0;
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        mark(itr.runStart(), itr.runEnd(), currentIndex++, STRIDE_0);
    }
}

IndexArray::
IndexArray(const SupportRegion & support)
{
    SupportRegion::ConstIterator itr;
    int64 currentIndex = 0;
    const int64 STRIDE_1 = 1;
    for (itr = support.begin(); itr != support.end(); itr.nextMarkedRun())
    {
        mark(itr.runStart(), itr.runEnd(), currentIndex, STRIDE_1);
        currentIndex += itr.runEnd() - itr.runStart() + 1;
    }
}

IndexArray::
IndexArray(const RLEBase<IASegment> & copyMe) :
    RLEBase<IASegment>(copyMe)
{
}

IndexArray & IndexArray::
operator=(const RLEBase<IASegment> & rhs)
{
    ((RLEBase<IASegment>&)*this) = rhs;
    return *this;
}

std::vector<int64> IndexArray::
indices() const
{
    std::vector<int64> vec(numRuns());
    
    AVLConstIterator itr = mSegmentTree.begin();
    int64 nn = 0;
    while (itr.hasNext())
        vec[nn++] = itr.next().mark();
    
    return vec;
}

std::vector<uint64> IndexArray::
lengths() const
{
    std::vector<uint64> vec(numRuns());
    
    AVLConstIterator itr = mSegmentTree.begin();
    int64 nn = 0;
    while (itr.hasNext())
        vec[nn++] = itr.next().length();
    
    return vec;
}

std::vector<uint64> IndexArray::
strides() const
{
    std::vector<uint64> vec(numRuns());
    
    AVLConstIterator itr = mSegmentTree.begin();
    int64 nn = 0;
    while (itr.hasNext())
        vec[nn++] = itr.next().data().stride();
    
    return vec;
}

std::vector<int64> IndexArray::
starts() const
{
    std::vector<int64> vec(numRuns());
    
    AVLConstIterator itr = mSegmentTree.begin();
    int64 nn = 0;
    while (itr.hasNext())
        vec[nn++] = itr.next().first();
    
    return vec;
}



#pragma mark *** IASegment ***


IASegment::
IASegment(int64 first, int64 last) :
    BaseSegment(first, last)
{
    // Warning: mMark is undefined
    if (length() == 1)
    {
        mData.mStride = 0;
    }
}

IASegment::
IASegment(int64 first, int64 last, DataType data) :
    BaseSegment(first, last),
    mData(data)
{
    if (length() == 1)
        mData.mStride = 0;
}

IASegment::
IASegment(const IASegment & copyMe, int64 first, int64 last) :
    BaseSegment(first, last),
    mData(copyMe.mData)
{
    if (first > copyMe.first() && copyMe.data().stride() != 0)
        mData.mIndex = copyMe.mData.mIndex +
            (first - copyMe.first())*copyMe.data().stride();
    if (length() == 1)
        mData.mStride = 0;
}

bool IASegment::
operator==(const IASegment & rhs) const
{
    return (mData == rhs.mData) &&
        (mFirst == rhs.mFirst) && (mLast == rhs.mLast);
}

void IASegment::
chopTail(int64 last)
{
    BaseSegment::chopTail(last);
    if (length() == 1)
        mData.mStride = 0;
}

void IASegment::
chopHead(int64 first)
{
    BaseSegment::chopHead(first);
    if (length() == 1)
        mData.mStride = 0;
}

IASegment IASegment::
tail(int64 first) const
{
    return IASegment(*this, first, mLast);
}

IASegment IASegment::
head(int64 last) const
{
    return IASegment(*this, mFirst, last);
}

void IASegment::
absorb(const IASegment & segment)
{
    assert(mFirst <= segment.last() + 1 && mLast >= segment.first()-1);
    /*std::cout << "(" << mFirst << ", " << mLast << ") absorbing ("
        << segment.first() << ", " << segment.last() << ")\n";*/
    if (segment.first() < mFirst)
    {
        // When absorbing a segment to the left, we take its starting index.
        mFirst = segment.first();
        mData.mIndex = segment.mData.mIndex;
    }
    if (segment.last() > mLast)
    {
        mLast = segment.last();
    }
    
    if (mData.stride() != segment.mData.stride() ||
        mData.index() != segment.mData.index())
    {
        // This happens if two length-one segments meet, or if a length-one
        // segment merges with a variable length segment.
        mData.mStride = 1;
    }
    
    // This should not happen, but just in case it does, protect the invariant.
    assert(length() > 1); // double-check.
    if (length() == 1)
        mData.mStride = 0;
}

bool IASegment::
canEat(const IASegment & segment) const
{
    assert(mFirst < segment.first());
    assert(mLast > segment.last());
    return (markAt(segment.first()) == segment.markAt(segment.first())) &&
        (markAt(segment.last()) == segment.markAt(segment.last()));
}

bool IASegment::
canMerge(const IASegment & segment) const
{
    /*std::cout << "Can (" << first() << ", " << last() << ") s " << mData.stride
        << " absorb (" << segment.first() << ", " << segment.last() << ") s "
        << segment.mData.stride << "?" << std::endl;*/
    if (mData.stride() == 0 && segment.mData.stride() == 0)
    {
        // I can eat another constant segment with the same index:
        // [3  ][3  ] becomes [3     ].
        if (mData.index() == segment.mData.index())
            return 1;
        
        // If both segments are length 1, they can still merge.
        if (length() == 1 && segment.length() == 1)
        {
            if (last() == segment.first() - 1 &&
                data().index() == segment.data().index()-1)
            {
                return 1;
            }
            else if (first() == segment.last() + 1 &&
                data().index() == segment.data().index() + 1)
            {
                return 1;
            }
        }
        return 0;
    }
    else if (mData.stride() == 1 && segment.mData.stride() == 1)
    {
        // I can eat another variable segment that melds seamlessly:
        // (1 2 3)(4 5 6) becomes (1 2 3 4 5 6).
        return mData.index() ==
            segment.data().index() + first() - segment.first();
    }
    
    // One segment is stride 0 and the other is stride 1.
    // if the variable segment could add on the single cell of the constant
    // segment.  (Remember: length-1 segments are constant (stride 0) by an
    // invariant of the IASegment class.  It's guaranteed that the
    // short segment here is constant.  Handy.
    if (length() == 1 || segment.length() == 1)
    {
        if (last() == segment.first()-1 &&
            markAt(last()) == segment.markAt(segment.first())-1)
        {
            // e.g. [0](1 2 3 4) or [1][2] or (0 1 2 3)[4]
            return 1;
        }
        else if (first() == segment.last()+1 &&
            markAt(first()) == segment.markAt(segment.last())+1)
        {
            // e.g. (1 2 3)[4] or [2][3]
            return 1;
        }
    }
    
    // All options for merging have been exhausted!
    return 0;
}

}; // namespace RLE


#endif
