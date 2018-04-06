/*
 *  SupportRegion-inl.h
 *  refactorRLE
 *
 *  Created by Paul Hansen on 2/8/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _SUPPORTREGION_

namespace RLE
{

SupportRegion::
SupportRegion() :
    RLEBase<SRSegment>()
{
}

template<class T>
SupportRegion::
SupportRegion(const DynamicRLE<T> & rle)
{
    typename DynamicRLE<T>::ConstIterator itr;
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
        mark(itr.runStart(), itr.runEnd());
}

SupportRegion::
SupportRegion(const RLEBase<SRSegment> & copyMe) :
    RLEBase<SRSegment>(copyMe)
{
}

SupportRegion & SupportRegion::
operator=(const RLEBase<SRSegment> & rhs)
{
    ((RLEBase<SRSegment> &) *this) = rhs;
    return *this;
}

void SupportRegion::
mark(int64 first, int64 last)
{
    RLEBase<SRSegment>::mark(first, last, true);
}

bool SupportRegion::
encloses(const SupportRegion & rhs) const
{
    SupportRegion::ConstIterator itr;
    
    for (itr = rhs.begin(); itr != rhs.end(); itr.nextMarkedRun())
    {
        SupportRegion::ConstIterator itr2 = location(itr.runStart());
        if (false == itr2.isMarked())
            return false;
        if (!itr2->encloses(*itr))
            return false;
    }
    return true;
}

void SupportRegion::
operator*=(const SupportRegion & rhs)
{
    // It'll be a little gross to be invalidating iterators with erasures, so
    // I'll do this the slow, stupid way.  Be warned.
    SupportRegion newThis;
    intersection(*this, rhs, newThis);
    *this = newThis;
}

void SupportRegion::
operator+=(const SupportRegion & rhs)
{
    SupportRegion::ConstIterator itr;
    for (itr = rhs.begin(); itr != rhs.end(); itr.nextMarkedRun())
        mark(itr.runStart(), itr.runEnd());
}

void SupportRegion::
operator-=(const SupportRegion & rhs)
{
    // I don't want to think of a fast way to do this right now.
    
    SupportRegion::ConstIterator itr;
    for (itr = rhs.begin(); itr != rhs.end(); itr.nextMarkedRun())
        erase(itr.runStart(), itr.runEnd());
}


#pragma mark *** Segment ***

SRSegment::
SRSegment(int64 first, int64 last) :
    BaseSegment(first, last)
{
}

SRSegment::
SRSegment(int64 first, int64 last, const DataType & data) :
    BaseSegment(first, last)
{
}

SRSegment::
SRSegment(const SRSegment & copyMe, int64 first, int64 last) :
    BaseSegment(first, last)
{
}

bool SRSegment::
operator==(const SRSegment & rhs) const
{
    return (mFirst == rhs.mFirst) && (mLast == rhs.mLast);
}

SRSegment SRSegment::
tail(int64 first) const
{
    return SRSegment(*this, first, mLast);
}

SRSegment SRSegment::
head(int64 last) const
{
    return SRSegment(*this, mFirst, last);
}

void SRSegment::
absorb(const SRSegment & segment)
{
    assert(mFirst <= segment.last() + 1 && mLast >= segment.first()-1);
    
    if (segment.first() < mFirst)
        mFirst = segment.first();
    if (segment.last() > mLast)
        mLast = segment.last();
}

#pragma mark *** Operations ***


void join(const SupportRegion & lhs, const SupportRegion & rhs,
    SupportRegion & result)
{
    result = lhs;
    
    SupportRegion::ConstIterator itr;
    for (itr = rhs.begin(); itr != rhs.end(); itr.nextMarkedRun())
        result.mark(itr.runStart(), itr.runEnd());
}

void difference(const SupportRegion & lhs, const SupportRegion & rhs,
    SupportRegion & result)
{
    result = lhs;
    
    SupportRegion::ConstIterator itr;
    for (itr = rhs.begin(); itr != rhs.end(); itr.nextMarkedRun())
        result.erase(itr.runStart(), itr.runEnd());
}

void intersection(const SupportRegion & lhs, const SupportRegion & rhs,
    SupportRegion & result)
{
    RLE::transform(lhs, rhs, result, std::logical_and<bool>(), false, false);
}


SupportRegion operator*(const SupportRegion & lhs, const SupportRegion & rhs)
{
    SupportRegion result;
    intersection(lhs, rhs, result);
    return result;
}

SupportRegion operator+(const SupportRegion & lhs, const SupportRegion & rhs)
{
    SupportRegion result;
    join(lhs, rhs, result);
    return result;
}

SupportRegion operator-(const SupportRegion & lhs, const SupportRegion & rhs)
{
    SupportRegion result;
    difference(lhs, rhs, result);
    return result;
}

template<class T>
SupportRegion operator==(const DynamicRLE<T> & lhs, const T & rhs)
{
    SupportRegion result;
    transform(lhs, rhs, result, equal_to<T>);
    return result;
}

template<class T>
SupportRegion operator!=(const DynamicRLE<T> & lhs, const T & rhs)
{
    SupportRegion result;
    transform(lhs, rhs, result, not_equal_to<T>);
    return result;
}


template<class RunLengthEncoded>
RunLengthEncoded restriction(const RunLengthEncoded & lhs,
    const SupportRegion & rhs)
{
    RunLengthEncoded rle;
    restriction(lhs, rhs, rle);
    return rle;
}
// ambiguous overload
//template<class RunLengthEncoded>
//RunLengthEncoded operator*(const SupportRegion & lhs,
//    const RunLengthEncoded & rhs)
//{
//    return restriction(rhs, lhs);
//}
//
//template<class RunLengthEncoded>
//RunLengthEncoded operator*(const RunLengthEncoded & lhs,
//    const SupportRegion & rhs)
//{
//    return restriction(lhs, rhs);
//}


}; // namespace RLE

#endif
