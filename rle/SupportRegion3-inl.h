/*
 *  SupportRegion3-inl.h
 *  rle_index
 *
 *  Created by Paul Hansen on 2/2/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _SUPPORTREGION3_

namespace RLE
{


SupportRegion3::
SupportRegion3() :
    RLEBase3<SRSegment>()
{
}

SupportRegion3::
SupportRegion3(bool hasX, bool hasY, bool hasZ) :
    RLEBase3<SRSegment>(hasX, hasY, hasZ)
{
}

template<class Scalar>
SupportRegion3::
SupportRegion3(const Scalar nonSymmetricDimensions[]) :
    RLEBase3<SRSegment>(nonSymmetricDimensions)
{
}

template<class Segment>
SupportRegion3::
SupportRegion3(const RLEBase3<Segment> & rle) :
    RLEBase3<SRSegment>(rle.nonSymmetricDimensions())
{
    typename RLEBase3<Segment>::ConstIterator itr;
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
        RLEBase<SRSegment>::mark(itr.runStart(), itr.runEnd(), true);
}

/*
template<class T>
SupportRegion3::
SupportRegion3(const DynamicRLE3<T> & rle) :
    RLEBase3<SRSegment>(rle.capacity())
{
    typename DynamicRLE<T>::ConstIterator itr;
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
        RLEBase<SRSegment>::mark(itr.runStart(), itr.runEnd(), true);
}

SupportRegion3::
SupportRegion3(const RLEBase3<SRSegment> & copyMe) :
    RLEBase3<SRSegment>(copyMe)
{
}
*/

SupportRegion3 & SupportRegion3::
operator=(const RLEBase3<SRSegment> & rhs)
{
    ((RLEBase3<SRSegment> &) *this) = rhs;
    return *this;
}

void SupportRegion3::
mark(int64 x0, int64 y0, int64 z0, int64 x1, int64 y1, int64 z1)
{
    RLEBase3<SRSegment>::mark(x0, y0, z0, x1, y1, z1, true);
}

void SupportRegion3::
mark(int64 x0, int64 y0, int64 z0)
{
    RLEBase3<SRSegment>::mark(x0, y0, z0, x0, y0, z0, true);
}

bool SupportRegion3::
encloses(const SupportRegion3 & rhs) const
{
    assert(stride(0) == rhs.stride(0) && stride(1) == rhs.stride(1)
        && stride(2) == rhs.stride(2));
    SupportRegion3::ConstIterator itr;
    
    for (itr = rhs.begin(); itr != rhs.end(); itr.nextMarkedRun())
    {
        SupportRegion3::ConstIterator itr2 =
            RLEBase<SRSegment>::location(itr.runStart());
        if (false == itr2.isMarked())
            return false;
        if (!itr2->encloses(*itr))
            return false;
    }
    return true;
}

void SupportRegion3::
operator*=(const SupportRegion3 & rhs)
{
    // It'll be a little gross to be invalidating iterators with erasures, so
    // I'll do this the slow, stupid way.  Be warned.
    SupportRegion3 sr = *this * rhs;
    *this = sr;
    
//    SupportRegion newThis;
//    intersection(*this, rhs, newThis);
//    *this = newThis;
}

void SupportRegion3::
operator+=(const SupportRegion3 & rhs)
{
    SupportRegion3 sr = *this + rhs;
    *this = sr;
//    SupportRegion::ConstIterator itr;
//    for (itr = rhs.begin(); itr != rhs.end(); itr.nextMarkedRun())
//        mark(itr.runStart(), itr.runEnd());
}

void SupportRegion3::
operator-=(const SupportRegion3 & rhs)
{
    SupportRegion3 sr = *this - rhs;
    *this = sr;
    // I don't want to think of a fast way to do this right now.
//    SupportRegion::ConstIterator itr;
//    for (itr = rhs.begin(); itr != rhs.end(); itr.nextMarkedRun())
//        erase(itr.runStart(), itr.runEnd());
}


#pragma mark *** Boolean operations ***

void join(const SupportRegion3 & lhs, const SupportRegion3 & rhs,
    SupportRegion3 & result)
{
    if (lhs.numRuns() > 0 && rhs.numRuns() > 0)
    {
        assert(lhs.stride(0) == rhs.stride(0) && lhs.stride(1) == rhs.stride(1)
            && lhs.stride(2) == rhs.stride(2));
    }
    
    result.clear();
    
    if (lhs.numRuns() > 0)
    {
        result.dimensions(lhs.nonSymmetricDimensions());
    }
    else
    {
        result.dimensions(rhs.nonSymmetricDimensions());
    }
    
    join((const SupportRegion&) lhs, (const SupportRegion&) rhs,
        (SupportRegion&) result);
}

// Predicates for transform(RLEBase3, RLEBase3, RLEBase3, BinaryOp).
// Used to implement difference(), since the STL doesn't have these in quite
// this form and I don't want to have to cobble them together.
static bool logical_a_and_not_b(bool a, bool b) { return a && (!b); }
static bool logical_not_a_and_b(bool a, bool b) { return (!a) && b; }

void difference(const SupportRegion3 & lhs, const SupportRegion3 & rhs,
    SupportRegion3 & result)
{
    result.clear();
    
    if (lhs.stride(0) == rhs.stride(0) && lhs.stride(1) == rhs.stride(1)
        && lhs.stride(2) == rhs.stride(2))
    {
        result.dimensions(lhs.nonSymmetricDimensions());
        difference((const SupportRegion&) lhs, (const SupportRegion&) rhs,
            (SupportRegion&) result);
    }
    else if (lhs.numDimensions() < rhs.numDimensions())
    {
        throw(std::logic_error("Cannot subtract higher-dimensional region from"
            " lower-dimensional region"));
        //transformHighLow(rhs, lhs, result, logical_not_a_and_b, false, false);
        //differenceLowHigh(lhs, rhs, result);
    }
    else if (lhs.numDimensions() > rhs.numDimensions())
    {
        transformHighLow(lhs, rhs, result, logical_a_and_not_b, false, false);
        //differenceHighLow(lhs, rhs, result);
    }
}

void intersection(const SupportRegion3 & lhs, const SupportRegion3 & rhs,
    SupportRegion3 & result)
{
    result.clear();
    
    if (lhs.stride(0) == rhs.stride(0) && lhs.stride(1) == rhs.stride(1)
        && lhs.stride(2) == rhs.stride(2))
    {
        result.dimensions(lhs.nonSymmetricDimensions());
        intersection((const SupportRegion&) lhs, (const SupportRegion&) rhs,
            (SupportRegion&) result);
    }
    else if (lhs.numDimensions() < rhs.numDimensions())
    {
        transformHighLow(rhs, lhs, result, std::logical_and<bool>(),
            false, false);
    }
    else if (lhs.numDimensions() > rhs.numDimensions())
    {
        transformHighLow(lhs, rhs, result, std::logical_and<bool>(),
            false, false);
    }
}

SupportRegion3 operator*(const SupportRegion3 & lhs, const SupportRegion3 & rhs)
{
    SupportRegion3 result;
    intersection(lhs, rhs, result);
    return result;
}

SupportRegion3 operator+(const SupportRegion3 & lhs, const SupportRegion3 & rhs)
{
    SupportRegion3 result;
    join(lhs, rhs, result);
    return result;
}

SupportRegion3 operator-(const SupportRegion3 & lhs, const SupportRegion3 & rhs)
{
    SupportRegion3 result;
    difference(lhs, rhs, result);
    return result;
}

template<class T>
SupportRegion3 operator==(const DynamicRLE3<T> & lhs, const T & rhs)
{
    SupportRegion3 result;
    transform(lhs, rhs, result, equal_to<T>);
    return result;
}

template<class T>
SupportRegion3 operator!=(const DynamicRLE3<T> & lhs, const T & rhs)
{
    SupportRegion3 result;
    transform(lhs, rhs, result, not_equal_to<T>);
    return result;
}

template<class RunLengthEncoded>
RunLengthEncoded restriction(const RunLengthEncoded & lhs,
    const SupportRegion3 & rhs)
{
    RunLengthEncoded rle;
    restriction(lhs, rhs, rle);
    return rle;
}
//
//template<class RunLengthEncoded>
//RunLengthEncoded operator*(const SupportRegion3 & lhs,
//    const RunLengthEncoded & rhs)
//{
//    return restriction(rhs, lhs);
//}
//
//template<class RunLengthEncoded>
//RunLengthEncoded operator*(const RunLengthEncoded & lhs,
//    const SupportRegion3 & rhs)
//{
//    return restriction(lhs, rhs);
//}

}; // namespace RLE



#endif
