/*
 *  MultiTransform.h
 *  rle
 *
 *  Created by Paul C Hansen on 10/12/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef MULTIAPPLY_H
#define MULTIAPPLY_H

namespace RLE
{

enum {
    IncludeGaps,
    ExcludeGaps
};
/**
 * Apply a callback to a set of RLE arrays.  Much like the mergeOverlapping()
 * function, multiApply() will pass segments of constant data to the callback,
 * along with the beginnings and ends of these piecewise constant runs.  While
 * mergeOverlapping() passes strides and indices to the callback, multiApply()
 * passes data.
 *
 * I wonder if this is really the same as mergeOverlapping().
 *
 * The callback should implement the following:
 *
 * Callback::operator()(int64 first, int64 last, Segment::MarkType const* data,
 *     int const* markedArrays, int numStreams)
 */
template<class Array, class Callback>
void multiApply(Array const* arrays[], int numArrays, Callback cb,
    int gapPolicy = RLE::IncludeGaps);

}; // namespace RLE

#include "MultiApply-inl.h"

#endif
