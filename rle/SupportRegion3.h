/*
 *  SupportRegion3.h
 *  rle_index
 *
 *  Created by Paul Hansen on 2/2/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _SUPPORTREGION3_
#define _SUPPORTREGION3_

#include "RLEBase3.h"
#include "SupportRegion.h"

namespace RLE
{

template<class T>
class DynamicRLE3;

class SupportRegion3 : public RLEBase3<SRSegment>
{
public:
    inline SupportRegion3();
    inline SupportRegion3(bool hasX, bool hasY, bool hasZ);
    
    template<class Scalar>
    SupportRegion3(const Scalar nonSymmetricDimensions[]);
        
    template<class Segment>
    SupportRegion3(const RLEBase3<Segment> & rle);
    
    //inline SupportRegion3(const RLEBase3<SRSegment> & copyMe);
    inline SupportRegion3 & operator=(const RLEBase3<SRSegment> & rhs);
        
    inline void mark(int64 x0, int64 y0, int64 z0, int64 x1, int64 y1, int64 z1);
    inline void mark(int64 x0, int64 y0, int64 z0);
    
    // presently only supports same-dimension arguments
    inline bool encloses(const SupportRegion3 & rhs) const;
    
    inline void operator*=(const SupportRegion3 & rhs);
    inline void operator+=(const SupportRegion3 & rhs);
    inline void operator-=(const SupportRegion3 & rhs);
};


#pragma mark *** Boolean operations ***

inline void join(const SupportRegion3 & lhs, const SupportRegion3 & rhs,
    SupportRegion3 & result);
inline void difference(const SupportRegion3 & lhs, const SupportRegion3 & rhs,
    SupportRegion3 & result);
inline void intersection(const SupportRegion3 & lhs, const SupportRegion3& rhs,
    SupportRegion3 & result);

inline SupportRegion3
operator+(const SupportRegion3 & lhs, const SupportRegion3& rhs);
inline SupportRegion3
operator-(const SupportRegion3 & lhs, const SupportRegion3& rhs);
inline SupportRegion3
operator*(const SupportRegion3 & lhs, const SupportRegion3& rhs);


template<class T>
SupportRegion3 operator==(const DynamicRLE3<T> & lhs, const T & rhs);
template<class T>
SupportRegion3 operator!=(const DynamicRLE3<T> & lhs, const T & rhs);

// Shortcuts for restriction
template<class RunLengthEncoded>
RunLengthEncoded restriction(const RunLengthEncoded & lhs,
    const SupportRegion3 & rhs);

// These don't work ("ambiguous overload")
//template<class RunLengthEncoded>
//RunLengthEncoded operator*(const SupportRegion3 & lhs,
//    const RunLengthEncoded & rhs);
//
//template<class RunLengthEncoded>
//RunLengthEncoded operator*(const RunLengthEncoded & lhs,
//    const SupportRegion3 & rhs);

}; // namespace RLE

#include "SupportRegion3-inl.h"

#endif
