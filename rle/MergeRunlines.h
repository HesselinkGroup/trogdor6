/*
 *  Runline.h
 *  rle_index
 *
 *  Created by Paul Hansen on 1/29/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _RUNLINE_
#define _RUNLINE_

#include <vector>

namespace RLE
{


/**
 * Merge runlines when all the runlines coincide (there are no gaps).
 *
 * This is the OLD implementation.  It is a little slower and a little bit
 * more wasteful of memory (requiring the separate index, stride and length
 * arrays).
 */
template<class Callback>
void merge(std::vector<int64> indexArray[],
    std::vector<uint64> strideArray[],
    std::vector<uint64> lengthArray[],
    uint64 numStreams,
    Callback callback);
    
/**
 * Merge runlines when all the runlines coincide (there are no gaps).
 *
 * This is the newer, faster implementation.
 */
template<class Callback, class IndArray>
void merge(IndArray const* arrays[], uint64 numStreams, Callback callback);


/**
 * Merge runlines when there may be gaps in places (not every runline covers
 * every position), but when all the arrays were the same size and are supposed
 * to overlap cell-by-cell.
 *
 * This is the OLD implementation; it's a little slower than the new one.
 */
template<class Callback>
void mergeOverlapping(std::vector<int64> indexArray[],
    std::vector<uint64> strideArray[],
    std::vector<uint64> lengthArray[],
    std::vector<int64> startArray[],
    uint64 numStreams,
    Callback callback);

/**
 * Merge runlines when there may be gaps in places (not every runline covers
 * every position), but when all the arrays were the same size and are supposed
 * to overlap cell-by-cell.
 *
 * This is the NEW implementation, faster and ultimately more thrifty with RAM.
 */
template<class Callback, class IndArray>
void mergeOverlapping(IndArray const* arrays[], uint64 numStreams,
    Callback callback);

// Nice default callback to use.  At least it does *something*.
struct PrintingCallback
{
    PrintingCallback() :
        mCurrentIndex(0)
    {
    }
    
    void operator()(const int64 indices[],
        const uint64 strides[],
        uint64 length,
        uint64 numStreams)
    {
        std::cout << "[" << mCurrentIndex << ", " << mCurrentIndex+length-1
            << "]\n";
        mCurrentIndex += length;
    }
    
    int64 mCurrentIndex;
};

struct ThoroughCallback
{
    ThoroughCallback()
    {
    }
    
    void operator()(const int64 indices[],
        const uint64 strides[],
        uint64 length,
        uint64 numStreams)
    {
        std::cout << "Indices: [";
        for (int nn = 0; nn < numStreams; nn++)
            std::cout << indices[nn] << " ";
        std::cout << "] Strides [";
        for (int nn = 0; nn < numStreams; nn++)
            std::cout << strides[nn] << " ";
        std::cout << "] Length " << length << "\n";
    }
};


}; // namespace RLE

#include "MergeRunlines-inl.h"

#endif
