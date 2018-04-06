/*
 *  IndexArray3.h
 *  rle_index
 *
 *  Created by Paul Hansen on 2/2/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _INDEXARRAY3_
#define _INDEXARRAY3_

#include "RLEBase3.h"
#include "IndexArray.h"

namespace RLE
{

template<class T>
class DynamicRLE3;

class SupportRegion3;

class IndexArray3 : public RLEBase3<IASegment>
{
public:
    inline IndexArray3();
    
    template<class T>
    IndexArray3(const DynamicRLE3<T> & rle);
    
    inline IndexArray3(const SupportRegion3 & support);
    
    inline IndexArray3(const RLEBase3<IASegment> & copyMe);
    inline IndexArray3 & operator=(const RLEBase3<IASegment> & rhs);
    
    void mark(int64 first, int64 last, int64 index, int64 stride)
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

#include "IndexArray3-inl.h"


#endif
