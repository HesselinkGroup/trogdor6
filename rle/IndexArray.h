/*
 *  IndexArray.h
 *  refactorRLE
 *
 *  Created by Paul Hansen on 2/8/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _INDEXARRAY_
#define _INDEXARRAY_

#include "RLEBase.h"
#include <vector>

namespace RLE
{
template<class T>
class DynamicRLE;

class SupportRegion;

class IASegment : public BaseSegment
{
public:
    struct IndexStride
    {
        IndexStride() {}
        IndexStride(int64 inIndex) :
            mIndex(inIndex), mStride(0) {}
        IndexStride(int64 inIndex, uint64 inStride) :
            mIndex(inIndex), mStride(inStride) {}
        
        int64 index() const { return mIndex; }
        uint64 stride() const { return mStride; }
        
        int64 mIndex;
        uint64 mStride;
        bool operator==(const IndexStride & rhs) const
        {
            return (mIndex == rhs.mIndex) && (mStride == rhs.mStride);
        }
    };
    
    // Concept requirements:
    typedef IndexStride DataType; // argument for constructor
    typedef const IndexStride & DataReturnType;
    typedef int64 MarkType;
    typedef int64 MarkReturnType; // return type for markAt(), Iterator::markAt(),...
//    enum { hasDefault = false };
    static bool isDefaultMark(const DataType & mark) { return false; }
    
    inline IASegment(int64 first, int64 last);
    inline IASegment(int64 first, int64 last, DataType data);
    inline IASegment(const IASegment & copyMe, int64 first, int64 last); // truncated seg.
    
    inline bool operator==(const IASegment & rhs) const;
    
    inline DataReturnType data() const { return mData; }
    inline MarkReturnType mark() const { return mData.index(); }
    inline MarkReturnType markAt(int64 position) const
    {
        return mData.index() + mData.stride()*(position - first());
    }
    
    // Override these two functions to protect the invariant that length-1
    // segments are always constant (stride = 0).
    inline void chopTail(int64 last);
    inline void chopHead(int64 first);
    
    inline IASegment tail(int64 first) const;
    inline IASegment head(int64 last) const;
    inline void absorb(const IASegment& segment);
    inline bool canEat(const IASegment & segment) const;
    inline bool canMerge(const IASegment & segment) const;
protected:
    DataType mData;
};

class IndexArray : public RLEBase<IASegment>
{
public:
    IndexArray() {}
    
    template<class T>
    IndexArray(const DynamicRLE<T> & rle);
    
    inline IndexArray(const SupportRegion & support);
    
    inline IndexArray(const RLEBase<IASegment> & copyMe);
    inline IndexArray & operator=(const RLEBase<IASegment> & rhs);
    
    void mark(int64 first, int64 last, int64 index, uint64 stride)
    {
        RLEBase<IASegment>::mark(first, last,
            IASegment::DataType(index, stride));
    }
    
    inline std::vector<int64> indices() const;
    inline std::vector<uint64> lengths() const;
    inline std::vector<uint64> strides() const;
    inline std::vector<int64> starts() const;
};

}; // namespace RLE

#include "IndexArray-inl.h"

#endif
