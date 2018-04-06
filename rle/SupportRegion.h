/*
 *  SupportRegion.h
 *  refactorRLE
 *
 *  Created by Paul Hansen on 2/8/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _SUPPORTREGION_
#define _SUPPORTREGION_

#include "RLEBase.h"

namespace RLE
{

class SRSegment : public BaseSegment
{
public:
    typedef bool DataType;
    typedef bool DataReturnType;
    typedef bool MarkType;
    typedef bool MarkReturnType;
//    enum { hasDefault = true };
    
    inline SRSegment(int64 first, int64 last);
    inline SRSegment(int64 first, int64 last, const DataType & data);
    inline SRSegment(const SRSegment & copyMe, int64 first, int64 last); // truncation
    
    inline bool operator==(const SRSegment & rhs) const;
    
    inline DataReturnType data() const { return true; }
    inline MarkReturnType mark() const { return true; }
    inline MarkReturnType markAt(int64 position) const { return true; }
    static bool isDefaultMark(const DataType & mark) { return mark == false; }
    
    inline SRSegment tail(int64 first) const;
    inline SRSegment head(int64 last) const;
    inline void absorb(const SRSegment& segment);
    inline bool canEat(const SRSegment & segment) const { return true; }
    inline bool canMerge(const SRSegment & segment) const { return true; }
};

template<class T>
class DynamicRLE;

class SupportRegion : public RLEBase<SRSegment>
{
public:
    inline SupportRegion();
    
    template<class T>
    SupportRegion(const DynamicRLE<T> & rle);
    
    inline SupportRegion(const RLEBase<SRSegment> & copyMe);
    inline SupportRegion & operator=(const RLEBase<SRSegment> & rhs);
    
    inline void mark(int64 first, int64 last);
    
    inline bool encloses(const SupportRegion & rhs) const;
    
    inline void operator*=(const SupportRegion & rhs);
    inline void operator+=(const SupportRegion & rhs);
    inline void operator-=(const SupportRegion & rhs);
};

#pragma mark *** Boolean operations ***

// useful predicate
template<class T>
bool equal_to(const T & lhs, const T & rhs)
{
    return lhs == rhs;
}

// another useful predicate
template<class T>
bool not_equal_to(const T & lhs, const T & rhs)
{
    return lhs != rhs;
}

inline void join(const SupportRegion & lhs, const SupportRegion & rhs,
    SupportRegion & result);
inline void difference(const SupportRegion & lhs, const SupportRegion & rhs,
    SupportRegion & result);
inline void intersection(const SupportRegion & lhs, const SupportRegion& rhs,
    SupportRegion & result);

inline SupportRegion
operator*(const SupportRegion & lhs, const SupportRegion & rhs);
inline SupportRegion
operator+(const SupportRegion & lhs, const SupportRegion & rhs);
inline SupportRegion
operator-(const SupportRegion & lhs, const SupportRegion & rhs);

template<class T>
SupportRegion operator==(const DynamicRLE<T> & lhs, const T & rhs);
template<class T>
SupportRegion operator!=(const DynamicRLE<T> & lhs, const T & rhs);

// Shortcuts for restriction
template<class RunLengthEncoded>
RunLengthEncoded restriction(const RunLengthEncoded & lhs,
    const SupportRegion & rhs);

// these don't work ("ambiguous overload")
//template<class RunLengthEncoded>
//RunLengthEncoded operator*(const SupportRegion & lhs,
//    const RunLengthEncoded & rhs);
//
//template<class RunLengthEncoded>
//RunLengthEncoded operator*(const RunLengthEncoded & lhs,
//    const SupportRegion & rhs);

}; // namespace RLE

#include "SupportRegion-inl.h"

#endif
