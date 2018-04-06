/*
 *  RLEBase-inl.h
 *  refactorRLE
 *
 *  Created by Paul Hansen on 2/8/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _RLEBASE_

#include <stdexcept>
#include <limits>
#include <boost/type_traits/is_arithmetic.hpp>

namespace RLE
{

template<typename S, typename T, bool condition>
class AssignIf
{
public:
    static void assign(S & variable, const T & argument)
    {
        variable = S(argument);
    }
};

template<typename S, typename T>
class AssignIf<S, T, false>
{
public:
    static void assign(S & variable, const T & argument)
    {
    }
};


template<class Segment>
RLEBase<Segment>::
RLEBase()
{
//    if (Segment::hasDefault)
//    {
//        AssignIf<typename Segment::MarkType, int,
//            boost::is_arithmetic<typename Segment::MarkType>::value >::assign(
//            mDefaultMark, 0);
//    }
}

//template<class Segment>
//RLEBase<Segment>::
//RLEBase(const typename Segment::DataType & defaultMark) :
//    mDefaultMark(defaultMark)
//{
//}

//template<class Segment>
//typename Segment::MarkReturnType RLEBase<Segment>::
//defaultMark() const
//{
//    return mDefaultMark;
//}

template<bool Condition>
struct ConditionalEquality
{
public:
    template<class T, class S>
    static bool isEqual(const T & lhs, const S & rhs)
    {
        return false;
    }
};

template<>
struct ConditionalEquality<true>
{
public:
    template<class T, class S>
    static bool isEqual(const T & lhs, const S & rhs)
    {
        return lhs == rhs;
    }
};


template<class Segment>
void RLEBase<Segment>::
clear()
{
    mSegmentTree.clear();
}

//template<class Segment>
//void RLEBase<Segment>::
//clear(const typename Segment::MarkType & newDefaultMark)
//{
//    mSegmentTree.clear();
//    mDefaultMark = newDefaultMark;
//}

template<class Segment>
void RLEBase<Segment>::
mark(int64 first, int64 last, const typename Segment::DataType & mark)
{
//    if (ConditionalEquality<Segment::hasDefault>::isEqual(mark, mDefaultMark))
//    {
//        erase(first, last);
//        return;
//    }
    assert(last >= first);
    if (Segment::isDefaultMark(mark))
    {
        erase(first, last);
        return;
    }
    Segment segment(first, last, mark);
    
    AVLIterator itr = mSegmentTree.location(segment);
    
    bool didAbsorbNew = false;
    if (itr.hasPrevious() && itr.peekPrevious().surrounds(segment))
    {
        // New segment is in the middle of an existing segment, with some
        // of the old segment on both sides.
        if (itr.peekPrevious().canEat(segment))
        {
            // New segment has exactly the same value as the old data in this
            // place, so we can skip the insertion altogether.
        }
        else
        {
            // New segment will divide an existing segment int64o two halves:
            Segment tail = itr.peekPrevious().tail(segment.last()+1);
            itr.peekPrevious().chopTail(segment.first() - 1);
            mSegmentTree.insert(tail);
            mSegmentTree.insert(segment);
        }
    }
    else
    {
        // Resolve the relationship with segments that touch on the left.
        if (itr.hasPrevious() && itr.peekPrevious().rightTouches(segment))
        {
            // A segment to the left is in contact with the new segment.  First
            // resolve this contact.
            if (itr.peekPrevious().canMerge(segment))
            {
                itr.peekPrevious().absorb(segment);
                didAbsorbNew = true;
            }
            else if (itr.peekPrevious().rightOverlaps(segment))
                itr.peekPrevious().chopTail(segment.first()-1);
        }
        
        // Move to the right and eat everything we cover completely.
        while (itr.hasNext() && segment.encloses(itr.peekNext()))
            itr.eraseNext();
        
        // Resolve the relationship with any segment that touches on the right.
        if (itr.hasNext() && itr.peekNext().leftTouches(segment))
        {
            if (itr.peekNext().canMerge(segment))
            {
                if (didAbsorbNew)
                {
                    assert(itr.peekPrevious().rightTouches(itr.peekNext()));
                    itr.peekPrevious().absorb(itr.peekNext());
                    itr.eraseNext();
                }
                else
                {
                    itr.peekNext().absorb(segment);
                    didAbsorbNew = true;
                }
            }
            else if (itr.peekNext().leftOverlaps(segment))
            {
                itr.peekNext().chopHead(segment.last()+1);
            }
        }
        
        if (didAbsorbNew == false)
        {
            mSegmentTree.insert(segment);
        }
    }
}

template<class Segment>
void RLEBase<Segment>::
erase(int64 first, int64 last)
{
    assert(last >= first);
    AVLIterator itr = mSegmentTree.location(first);
    
    if (itr.hasPrevious() && itr.peekPrevious().surrounds(first, last))
    {
        // New segment will divide an existing segment int64o two halves:
        Segment tail = itr.peekPrevious().tail(last+1);
        itr.peekPrevious().chopTail(first - 1);
        mSegmentTree.insert(tail);
    }
    else
    {
        // Truncate any segment that touches on the left.
        if (itr.hasPrevious() && itr.peekPrevious().rightOverlaps(first, last))
            itr.peekPrevious().chopTail(first-1);
        
        // Erase all the completely overlapped segments
        while (itr.hasNext() && itr.peekNext().last() <= last)
            itr.eraseNext();
        
        // Truncate any segment that touches on the right.
        if (itr.hasNext() && itr.peekNext().leftOverlaps(first, last))
            itr.peekNext().chopHead(last+1);
    }
}


template<class Segment>
bool RLEBase<Segment>::
isMarked(int64 index) const
{
    AVLConstIterator itr(mSegmentTree.location(index));
    
    if (itr.hasPrevious() && itr.peekPrevious().last() >= index)
        return 1;
    return (itr.hasNext() && itr.peekNext().first() == index);
}

template<class Segment>
typename Segment::MarkReturnType RLEBase<Segment>::
markAt(int64 index) const
{
    AVLConstIterator itr(mSegmentTree.location(index));
    
    if (itr.hasPrevious() && itr.peekPrevious().last() >= index)
        return itr.peekPrevious().markAt(index);
    else if (itr.hasNext() && itr.peekNext().first() == index)
        return itr.peekNext().markAt(index);
    
    throw(std::logic_error("No mark"));
}

template<class Segment>
typename Segment::MarkType RLEBase<Segment>::
markAt(int64 index, const typename Segment::MarkType & defaultMark) const
{
    AVLConstIterator itr(mSegmentTree.location(index));
    
    if (itr.hasPrevious() && itr.peekPrevious().last() >= index)
        return itr.peekPrevious().markAt(index);
    else if (itr.hasNext() && itr.peekNext().first() == index)
        return itr.peekNext().markAt(index);
    
    return defaultMark;
}

template<class Segment>
typename Segment::MarkReturnType RLEBase<Segment>::
operator()(int64 index) const
{
    AVLConstIterator itr(mSegmentTree.location(index));
    
    if (itr.hasPrevious() && itr.peekPrevious().last() >= index)
        return itr.peekPrevious().markAt(index);
    else if (itr.hasNext() && itr.peekNext().first() == index)
        return itr.peekNext().markAt(index);
    
    throw(std::logic_error("No mark"));
}

template<class Segment>
typename Segment::MarkType RLEBase<Segment>::
at(int64 index) const
{
    AVLConstIterator itr(mSegmentTree.location(index));
    
    if (itr.hasPrevious() && itr.peekPrevious().last() >= index)
        return itr.peekPrevious().markAt(index);
    else if (itr.hasNext() && itr.peekNext().first() == index)
        return itr.peekNext().markAt(index);
    
    return typename Segment::MarkType();
}

template<class Segment>
uint64 RLEBase<Segment>::
bytes() const
{
    return sizeof(mSegmentTree) + mSegmentTree.bytes();
}

template<class Segment>
uint64 RLEBase<Segment>::
numRuns() const
{
    return mSegmentTree.size();
}

template<class Segment>
uint64 RLEBase<Segment>::
length() const
{
    uint64 totalLength = 0;
    
    AVLConstIterator itr = mSegmentTree.begin();
    
    while (itr.hasNext())
        totalLength += itr.next().length();
    
    return totalLength;
}

template<class Segment>
bool RLEBase<Segment>::
operator==(const RLEBase<Segment> & rhs) const
{
    return rhs.mSegmentTree == mSegmentTree;
}

template<class Segment>
bool RLEBase<Segment>::
operator!=(const RLEBase<Segment> & rhs) const
{
    return !(rhs.mSegmentTree == mSegmentTree);
}

template<class Segment>
typename RLEBase<Segment>::Iterator RLEBase<Segment>::
begin()
{
    return Iterator(this);
}

template<class Segment>
typename RLEBase<Segment>::Iterator RLEBase<Segment>::
end()
{
    const int64 UNUSED_ARGUMENT_MEANS_END = 0;
    return Iterator(this, UNUSED_ARGUMENT_MEANS_END);
}

template<class Segment>
typename RLEBase<Segment>::Iterator RLEBase<Segment>::
unmarkedBegin()
{
    Iterator itr(this);
    if (mSegmentTree.size() > 0)
        itr.move(-1);
    return itr;
}

template<class Segment>
typename RLEBase<Segment>::Iterator RLEBase<Segment>::
location(int64 position)
{
    Iterator itr(this);
    itr.seek(position);
    return itr;
}

template<class Segment>
typename RLEBase<Segment>::ConstIterator RLEBase<Segment>::
begin() const
{
    return ConstIterator(this);
}

template<class Segment>
typename RLEBase<Segment>::ConstIterator RLEBase<Segment>::
end() const
{
    const int64 UNUSED_ARGUMENT_MEANS_END = 0;
    return ConstIterator(this, UNUSED_ARGUMENT_MEANS_END);
}

template<class Segment>
typename RLEBase<Segment>::ConstIterator RLEBase<Segment>::
unmarkedBegin() const
{
    ConstIterator itr(this);
    if (mSegmentTree.size() > 0)
        itr.move(-1);
    return itr;
}

template<class Segment>
typename RLEBase<Segment>::ConstIterator RLEBase<Segment>::
location(int64 position) const
{
    ConstIterator itr(this);
    itr.seek(position);
    return itr;
}

template<class Segment>
typename RLEBase<Segment>::ConstIterator RLEBase<Segment>::
segment(int64 which) const
{
    ConstIterator itr(begin());
    
    for (int64 nn = 0; nn < which && itr != end(); nn++)
        itr.nextMarkedRun();
    
    return itr;
}


#pragma mark *** Iterator ***

template<class Segment>
RLEBase<Segment>::Iterator::
Iterator() :
    mData(0L)
{
}

template<class Segment>
RLEBase<Segment>::Iterator::
Iterator(RLEBase* data) :
    mData(data)
{
    mIterator = data->mSegmentTree.begin();
    if (mIterator.hasNext())
        mCurrentPosition = mIterator.peekNext().first();
    else
        mCurrentPosition = 0; // uh... default start.  Yah.
}

template<class Segment>
RLEBase<Segment>::Iterator::
Iterator(RLEBase* data, int64 unused) :
    mData(data)
{
    mIterator = data->mSegmentTree.end();
    if (mIterator.hasPrevious())
        mCurrentPosition = mIterator.peekPrevious().last()+1;
    else
        mCurrentPosition = 0; // default end position
}



template<class Segment>
typename RLEBase<Segment>::Iterator RLEBase<Segment>::Iterator::
operator++(int unused)
{
    Iterator returnValue(*this);
    move(1);
    return returnValue;
}

template<class Segment>
typename RLEBase<Segment>::Iterator & RLEBase<Segment>::Iterator::
operator++()
{
    move(1);
    return *this;
}

template<class Segment>
typename RLEBase<Segment>::Iterator RLEBase<Segment>::Iterator::
operator--(int unused)
{
    Iterator returnValue(*this);
    move(-1);
    return returnValue;
}

template<class Segment>
typename RLEBase<Segment>::Iterator & RLEBase<Segment>::Iterator::
operator--()
{
    move(-1);
    return *this;
}

template<class Segment>
void RLEBase<Segment>::Iterator::
seek(int64 position)
{
    mCurrentPosition = position;
    mIterator.seek(mCurrentPosition);
}

template<class Segment>
void RLEBase<Segment>::Iterator::
move(int64 offset)
{
    mCurrentPosition += offset;
    mIterator.seek(mCurrentPosition);
}

template<class Segment>
int64 RLEBase<Segment>::Iterator::
position() const
{
    return mCurrentPosition;
}

template<class Segment>
bool RLEBase<Segment>::Iterator::
valid() const
{
    return mIterator.valid();
}

template<class Segment>
bool RLEBase<Segment>::Iterator::
isMarked() const
{
    if (mIterator.hasPrevious() &&
        mCurrentPosition <= mIterator.peekPrevious().last())
    {
        return 1;
    }
    else if (mIterator.hasNext() &&
        mCurrentPosition >= mIterator.peekNext().first())
    {
        return 1;
    }
    return 0;
}

template<class Segment>
Segment & RLEBase<Segment>::Iterator::
operator*() const
{
    if (mIterator.hasPrevious() &&
        mCurrentPosition <= mIterator.peekPrevious().last())
    {
        return mIterator.peekPrevious();
    }
    else if (mIterator.hasNext() &&
        mCurrentPosition >= mIterator.peekNext().first())
    {
        return mIterator.peekNext();
    }
    throw(std::logic_error("No current segment"));
}

template<class Segment>
Segment * RLEBase<Segment>::Iterator::
operator->() const
{
    if (mIterator.hasPrevious() &&
        mCurrentPosition <= mIterator.peekPrevious().last())
    {
        return &mIterator.peekPrevious();
    }
    else if (mIterator.hasNext() &&
        mCurrentPosition >= mIterator.peekNext().first())
    {
        return &mIterator.peekNext();
    }
    throw(std::logic_error("No current segment"));
}





template<class Segment>
void RLEBase<Segment>::Iterator::
nextRun()
{
    if (mIterator.hasPrevious())
    {
        if (mCurrentPosition <= mIterator.peekPrevious().last())
            seek(mIterator.peekPrevious().last()+1);
        else if (mIterator.hasNext())
        {
            if (mCurrentPosition >= mIterator.peekNext().first())
                seek(mIterator.peekNext().last()+1);
            else
                seek(mIterator.peekNext().first());
        }
        else
        {
            //std::cerr << "In last run already (a).\n";
            mCurrentPosition = mIterator.peekPrevious().last()+1;
            assert(*this == mData->end());
        }
    }
    else if (mIterator.hasNext())
    {
        if (mCurrentPosition == mIterator.peekNext().first())
            seek(mIterator.peekNext().last()+1);
        else if (mCurrentPosition < mIterator.peekNext().first())
            seek(mIterator.peekNext().first());
        else
            throw(std::logic_error("Iterator has invalid position!"));
    }
    else
    {
        //std::cerr << "In last run already (c).\n";
        assert(*this == mData->end());
    }
}

template<class Segment>
void RLEBase<Segment>::Iterator::
previousRun()
{
    if (mIterator.hasPrevious())
    {
        if (mCurrentPosition <= mIterator.peekPrevious().last())
            seek(mIterator.peekPrevious().first()-1);
        else if (mIterator.hasNext())
        {
            if (mCurrentPosition >= mIterator.peekNext().first())
                seek(mIterator.peekNext().first()-1);
            else
                seek(mIterator.peekPrevious().last());
        }
        else
            seek(mIterator.peekPrevious().last());
    }
    else if (mIterator.hasNext())
    {
        if (mCurrentPosition == mIterator.peekNext().first())
            seek(mIterator.peekNext().first()-1);
        else
            throw(std::logic_error("Iterator has invalid position!"));
    }
    else
    {
        //std::cerr << "In first run already (b).\n";
        assert(*this == mData->unmarkedBegin());
    }
}

template<class Segment>
void RLEBase<Segment>::Iterator::
nextMarkedRun()
{
    if (mIterator.hasNext())
    {
        if (mCurrentPosition < mIterator.peekNext().first())
            mCurrentPosition = mIterator.peekNext().first();
        else
        {
            mIterator++; // will only move if not at end yet.
            if (mIterator.hasNext())
                mCurrentPosition = mIterator.peekNext().first();
            else
            {
                //std::cerr << "At end.\n";
                mCurrentPosition = mIterator.peekPrevious().last()+1;
                assert(*this == mData->end());
            }
        }
    }
    else
    {
        mCurrentPosition = mIterator.peekPrevious().last()+1;
        //std::cerr << "At end.\n";
        assert(*this == mData->end());
    }
}

template<class Segment>
void RLEBase<Segment>::Iterator::
previousMarkedRun()
{
    if (mIterator.hasPrevious())
    {
        if (mCurrentPosition > mIterator.peekPrevious().last())
            mCurrentPosition = mIterator.peekPrevious().last();
        else
        {
            mIterator--;
            if (mIterator.hasPrevious())
                mCurrentPosition = mIterator.peekPrevious().last();
            else
            {
                //std::cerr << "At beginning.\n";
                mCurrentPosition = mIterator.peekNext().first();
                assert(*this == mData->begin());
            }
        }
    }
    else
    {
        mCurrentPosition = mIterator.peekNext().first();
        //std::cerr << "At beginning.\n";
        assert(*this == mData->begin());
    }
}


template<class Segment>
int64 RLEBase<Segment>::Iterator::
runStart() const
{
    if (mIterator.hasPrevious())
    {
        if (mCurrentPosition <= mIterator.peekPrevious().last())
        {
            return mIterator.peekPrevious().first();
        }
        else if (mIterator.hasNext() &&
            mCurrentPosition >= mIterator.peekNext().first())
        {
            return mIterator.peekNext().first();
        }
        else
            return mIterator.peekPrevious().last()+1;
    }
    else if (mIterator.hasNext() &&
        mCurrentPosition >= mIterator.peekNext().first())
    {
        return mIterator.peekNext().first();
    }
    else
    {
        //std::cerr << "In unbounded run.\n";
        throw(std::logic_error("Length undefined in unbounded run.\n"));
    }
}

template<class Segment>
int64 RLEBase<Segment>::Iterator::
runEnd() const
{
    if (mIterator.hasPrevious())
    {
        if (mCurrentPosition <= mIterator.peekPrevious().last())
        {
            return mIterator.peekPrevious().last();
        }
        else if (mIterator.hasNext())
        {
            if (mCurrentPosition >= mIterator.peekNext().first())
            {
                return mIterator.peekNext().last();
            }
            else
                return mIterator.peekNext().first()-1;
        }
        else
        {
            //std::cerr << "In unbounded run.\n";
            throw(std::logic_error("Length undefined in unbounded run.\n"));
        }
    }
    else if (mIterator.hasNext())
    {
        if (mCurrentPosition >= mIterator.peekNext().first())
            return mIterator.peekNext().last();
        else
            return mIterator.peekNext().first()-1;
    }
    else
    {
        //std::cerr << "In unbounded run.\n";
        throw(std::logic_error("Length undefined in unbounded run.\n"));
        //return -1;
    }
}

template<class Segment>
bool RLEBase<Segment>::Iterator::
runIsBoundedAtStart() const
{
    if (isMarked())
        return true;
    return mIterator.hasPrevious();
}

template<class Segment>
bool RLEBase<Segment>::Iterator::
runIsBoundedAtEnd() const
{
    if (isMarked())
        return true;
    return mIterator.hasNext();
}

template<class Segment>
bool RLEBase<Segment>::Iterator::
operator==(const RLEBase<Segment>::Iterator & rhs) const
{
    return (mData == rhs.mData &&
        mCurrentPosition == rhs.mCurrentPosition);
}

template<class Segment>
bool RLEBase<Segment>::Iterator::
operator!=(const RLEBase<Segment>::Iterator & rhs) const
{
    return (mData != rhs.mData ||
        mCurrentPosition != rhs.mCurrentPosition);
}

template<class Segment>
bool RLEBase<Segment>::Iterator::
operator<(const RLEBase<Segment>::Iterator & rhs) const
{
    return (mData == rhs.mData &&
        mCurrentPosition < rhs.mCurrentPosition);
}
template<class Segment>
bool RLEBase<Segment>::Iterator::
operator>(const RLEBase<Segment>::Iterator & rhs) const
{
    return (mData == rhs.mData &&
        mCurrentPosition > rhs.mCurrentPosition);
}
template<class Segment>
bool RLEBase<Segment>::Iterator::
operator<=(const RLEBase<Segment>::Iterator & rhs) const
{
    return (mData == rhs.mData &&
        mCurrentPosition <= rhs.mCurrentPosition);
}
template<class Segment>
bool RLEBase<Segment>::Iterator::
operator>=(const RLEBase<Segment>::Iterator & rhs) const
{
    return (mData == rhs.mData &&
        mCurrentPosition >= rhs.mCurrentPosition);
}

template<class Segment>
typename Segment::MarkReturnType RLEBase<Segment>::Iterator::
mark() const
{
    if (mIterator.hasPrevious() &&
        mCurrentPosition <= mIterator.peekPrevious().last())
    {
        return mIterator.peekPrevious().markAt(mCurrentPosition);
    }
    else if (mIterator.hasNext() &&
        mCurrentPosition >= mIterator.peekNext().first())
    {
        return mIterator.peekNext().markAt(mCurrentPosition);
    }
    
    throw(std::logic_error("No current segment"));
    
//    if (Segment::hasDefault)
//        return mData->defaultMark();
}

#pragma mark *** ConstIterator ***

template<class Segment>
RLEBase<Segment>::ConstIterator::
ConstIterator() :
    mData(0L)
{
}

template<class Segment>
RLEBase<Segment>::ConstIterator::
ConstIterator(const RLEBase* data) :
    mData(data)
{
    mIterator = data->mSegmentTree.begin();
    if (mIterator.hasNext())
        mCurrentPosition = mIterator.peekNext().first();
    else
        mCurrentPosition = 0; // uh... default start.  Yah.
}

template<class Segment>
RLEBase<Segment>::ConstIterator::
ConstIterator(const RLEBase* data, int64 unused) :
    mData(data)
{
    mIterator = data->mSegmentTree.end();
    if (mIterator.hasPrevious())
        mCurrentPosition = mIterator.peekPrevious().last()+1;
    else
        mCurrentPosition = 0; // default end position
}

template<class Segment>
RLEBase<Segment>::ConstIterator::
ConstIterator(const Iterator & itr) :
    mData(itr.mData),
    mIterator(itr.mIterator),
    mCurrentPosition(itr.mCurrentPosition)
{
}

template<class Segment>
typename RLEBase<Segment>::ConstIterator RLEBase<Segment>::ConstIterator::
operator++(int unused)
{
    ConstIterator returnValue(*this);
    move(1);
    return returnValue;
}

template<class Segment>
typename RLEBase<Segment>::ConstIterator & RLEBase<Segment>::ConstIterator::
operator++()
{
    move(1);
    return *this;
}

template<class Segment>
typename RLEBase<Segment>::ConstIterator RLEBase<Segment>::ConstIterator::
operator--(int unused)
{
    ConstIterator returnValue(*this);
    move(-1);
    return returnValue;
}

template<class Segment>
typename RLEBase<Segment>::ConstIterator & RLEBase<Segment>::ConstIterator::
operator--()
{
    move(-1);
    return *this;
}

template<class Segment>
void RLEBase<Segment>::ConstIterator::
seek(int64 position)
{
    mCurrentPosition = position;
    mIterator.seek(mCurrentPosition);
}

template<class Segment>
void RLEBase<Segment>::ConstIterator::
move(int64 offset)
{
    mCurrentPosition += offset;
    mIterator.seek(mCurrentPosition);
}

template<class Segment>
int64 RLEBase<Segment>::ConstIterator::
position() const
{
    return mCurrentPosition;
}

template<class Segment>
bool RLEBase<Segment>::ConstIterator::
valid() const
{
    return mIterator.valid();
}

template<class Segment>
bool RLEBase<Segment>::ConstIterator::
isMarked() const
{
    if (mIterator.hasPrevious() &&
        mCurrentPosition <= mIterator.peekPrevious().last())
    {
        return 1;
    }
    else if (mIterator.hasNext() &&
        mCurrentPosition >= mIterator.peekNext().first())
    {
        return 1;
    }
    return 0;
}

template<class Segment>
const Segment & RLEBase<Segment>::ConstIterator::
operator*() const
{
    if (mIterator.hasPrevious() &&
        mCurrentPosition <= mIterator.peekPrevious().last())
    {
        return mIterator.peekPrevious();
    }
    else if (mIterator.hasNext() &&
        mCurrentPosition >= mIterator.peekNext().first())
    {
        return mIterator.peekNext();
    }
    throw(std::logic_error("No current segment"));
}

template<class Segment>
const Segment * RLEBase<Segment>::ConstIterator::
operator->() const
{
    if (mIterator.hasPrevious() &&
        mCurrentPosition <= mIterator.peekPrevious().last())
    {
        return &mIterator.peekPrevious();
    }
    else if (mIterator.hasNext() &&
        mCurrentPosition >= mIterator.peekNext().first())
    {
        return &mIterator.peekNext();
    }
    throw(std::logic_error("No current segment"));
}

template<class Segment>
void RLEBase<Segment>::ConstIterator::
nextRun()
{
    if (mIterator.hasPrevious())
    {
        if (mCurrentPosition <= mIterator.peekPrevious().last())
            seek(mIterator.peekPrevious().last()+1);
        else if (mIterator.hasNext())
        {
            if (mCurrentPosition >= mIterator.peekNext().first())
                seek(mIterator.peekNext().last()+1);
            else
                seek(mIterator.peekNext().first());
        }
        else
        {
            //std::cerr << "In last run already (a).\n";
            mCurrentPosition = mIterator.peekPrevious().last()+1;
            assert(*this == mData->end());
        }
    }
    else if (mIterator.hasNext())
    {
        if (mCurrentPosition == mIterator.peekNext().first())
            seek(mIterator.peekNext().last()+1);
        else if (mCurrentPosition < mIterator.peekNext().first())
            seek(mIterator.peekNext().first());
        else
            throw(std::logic_error("Iterator has invalid position!"));
    }
    else
    {
        //std::cerr << "In last run already (c).\n";
        assert(*this == mData->end());
    }
}

template<class Segment>
void RLEBase<Segment>::ConstIterator::
previousRun()
{
    if (mIterator.hasPrevious())
    {
        if (mCurrentPosition <= mIterator.peekPrevious().last())
            seek(mIterator.peekPrevious().first()-1);
        else if (mIterator.hasNext())
        {
            if (mCurrentPosition >= mIterator.peekNext().first())
                seek(mIterator.peekNext().first()-1);
            else
                seek(mIterator.peekPrevious().last());
        }
        else
            seek(mIterator.peekPrevious().last());
    }
    else if (mIterator.hasNext())
    {
        if (mCurrentPosition == mIterator.peekNext().first())
            seek(mIterator.peekNext().first()-1);
        else
            throw(std::logic_error("Iterator has invalid position!"));
    }
    else
    {
        //std::cerr << "In first run already (b).\n";
        assert(*this == mData->unmarkedBegin());
    }
}

template<class Segment>
void RLEBase<Segment>::ConstIterator::
nextMarkedRun()
{
    if (mIterator.hasNext())
    {
        if (mCurrentPosition < mIterator.peekNext().first())
            mCurrentPosition = mIterator.peekNext().first();
        else
        {
            mIterator++; // will only move if not at end yet.
            if (mIterator.hasNext())
                mCurrentPosition = mIterator.peekNext().first();
            else
            {
                //std::cerr << "At end.\n";
                mCurrentPosition = mIterator.peekPrevious().last()+1;
                assert(*this == mData->end());
            }
        }
    }
    else
    {
        mCurrentPosition = mIterator.peekPrevious().last()+1;
        //std::cerr << "At end.\n";
        assert(*this == mData->end());
    }
}

template<class Segment>
void RLEBase<Segment>::ConstIterator::
previousMarkedRun()
{
    if (mIterator.hasPrevious())
    {
        if (mCurrentPosition > mIterator.peekPrevious().last())
            mCurrentPosition = mIterator.peekPrevious().last();
        else
        {
            mIterator--;
            if (mIterator.hasPrevious())
                mCurrentPosition = mIterator.peekPrevious().last();
            else
            {
                //std::cerr << "At beginning.\n";
                mCurrentPosition = mIterator.peekNext().first();
                assert(*this == mData->begin());
            }
        }
    }
    else
    {
        mCurrentPosition = mIterator.peekNext().first();
        //std::cerr << "At beginning.\n";
        assert(*this == mData->begin());
    }
}


template<class Segment>
int64 RLEBase<Segment>::ConstIterator::
runStart() const
{
    if (mIterator.hasPrevious())
    {
        if (mCurrentPosition <= mIterator.peekPrevious().last())
        {
            return mIterator.peekPrevious().first();
        }
        else if (mIterator.hasNext() &&
            mCurrentPosition >= mIterator.peekNext().first())
        {
            return mIterator.peekNext().first();
        }
        else
            return mIterator.peekPrevious().last()+1;
    }
    else if (mIterator.hasNext() &&
        mCurrentPosition >= mIterator.peekNext().first())
    {
        return mIterator.peekNext().first();
    }
    else
    {
        //std::cerr << "In unbounded run.\n";
        throw(std::logic_error("Length undefined in unbounded run.\n"));
    }
}

template<class Segment>
int64 RLEBase<Segment>::ConstIterator::
runEnd() const
{
    if (mIterator.hasPrevious())
    {
        if (mCurrentPosition <= mIterator.peekPrevious().last())
        {
            return mIterator.peekPrevious().last();
        }
        else if (mIterator.hasNext())
        {
            if (mCurrentPosition >= mIterator.peekNext().first())
            {
                return mIterator.peekNext().last();
            }
            else
                return mIterator.peekNext().first()-1;
        }
        else
        {
            //std::cerr << "In unbounded run.\n";
            throw(std::logic_error("Length undefined in unbounded run.\n"));
        }
    }
    else if (mIterator.hasNext())
    {
        if (mCurrentPosition >= mIterator.peekNext().first())
            return mIterator.peekNext().last();
        else
            return mIterator.peekNext().first()-1;
    }
    else
    {
        //std::cerr << "In unbounded run.\n";
        throw(std::logic_error("Length undefined in unbounded run.\n"));
        //return -1;
    }
}

template<class Segment>
bool RLEBase<Segment>::ConstIterator::
runIsBoundedAtStart() const
{
    if (isMarked())
        return true;
    return mIterator.hasPrevious();
}

template<class Segment>
bool RLEBase<Segment>::ConstIterator::
runIsBoundedAtEnd() const
{
    if (isMarked())
        return true;
    return mIterator.hasNext();
}

template<class Segment>
bool RLEBase<Segment>::ConstIterator::
operator==(const RLEBase<Segment>::ConstIterator & rhs) const
{
    return (mData == rhs.mData &&
        mCurrentPosition == rhs.mCurrentPosition);
}

template<class Segment>
bool RLEBase<Segment>::ConstIterator::
operator!=(const RLEBase<Segment>::ConstIterator & rhs) const
{
    return (mData != rhs.mData ||
        mCurrentPosition != rhs.mCurrentPosition);
}

template<class Segment>
bool RLEBase<Segment>::ConstIterator::
operator<(const RLEBase<Segment>::ConstIterator & rhs) const
{
    return (mData == rhs.mData &&
        mCurrentPosition < rhs.mCurrentPosition);
}
template<class Segment>
bool RLEBase<Segment>::ConstIterator::
operator>(const RLEBase<Segment>::ConstIterator & rhs) const
{
    return (mData == rhs.mData &&
        mCurrentPosition > rhs.mCurrentPosition);
}
template<class Segment>
bool RLEBase<Segment>::ConstIterator::
operator<=(const RLEBase<Segment>::ConstIterator & rhs) const
{
    return (mData == rhs.mData &&
        mCurrentPosition <= rhs.mCurrentPosition);
}
template<class Segment>
bool RLEBase<Segment>::ConstIterator::
operator>=(const RLEBase<Segment>::ConstIterator & rhs) const
{
    return (mData == rhs.mData &&
        mCurrentPosition >= rhs.mCurrentPosition);
}

template<class Segment>
typename Segment::MarkReturnType RLEBase<Segment>::ConstIterator::
mark() const
{
    if (mIterator.hasPrevious() &&
        mCurrentPosition <= mIterator.peekPrevious().last())
    {
        return mIterator.peekPrevious().markAt(mCurrentPosition);
    }
    else if (mIterator.hasNext() &&
        mCurrentPosition >= mIterator.peekNext().first())
    {
        return mIterator.peekNext().markAt(mCurrentPosition);
    }
//    if (Segment::hasDefault)
//        return mData->defaultMark();
//    else
    throw(std::logic_error("No current segment"));
}


#pragma mark *** Operations ***

template<class Segment>
void translate(const RLEBase<Segment> & rle, RLEBase<Segment> & result,
    int64 distance)
{
    result.clear();
//    if (Segment::hasDefault)
//        result.clear(rle.defaultMark());
//    else
//        result.clear();
    
    typename RLEBase<Segment>::ConstIterator itr;
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        result.mark(itr.runStart()+distance, itr.runEnd()+distance,
            itr->data());
    }
}

template<class InSegment, class OutSegment, class UnaryFunction>
void transform(const RLEBase<InSegment> & rle, RLEBase<OutSegment> & result,
    UnaryFunction op)
{
    result.clear();
    
    typename RLEBase<InSegment>::ConstIterator itr;
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
        result.mark(itr.runStart(), itr.runEnd(), op(itr.mark()));
}

template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transform(const RLEBase<InSeg1> & lhs, const RLEBase<InSeg2> & rhs,
    RLEBase<OutSeg> & result, BinaryFunction op)
{
    result.clear();
    
    typename RLEBase<InSeg1>::ConstIterator itrLeft = lhs.unmarkedBegin();
    typename RLEBase<InSeg2>::ConstIterator itrRight = rhs.unmarkedBegin();
    typename RLEBase<InSeg1>::ConstIterator endLeft = lhs.end();
    typename RLEBase<InSeg2>::ConstIterator endRight = rhs.end();
    
    int64 currentPosition = std::min(itrLeft.position(), itrRight.position());
    
    itrLeft.seek(currentPosition);
    itrRight.seek(currentPosition);
    
    while (itrLeft != endLeft || itrRight != endRight)
    {
        int64 last = std::numeric_limits<int64>::max();
        if (itrLeft != endLeft)
            last = std::min(last, itrLeft.runEnd());
        if (itrRight != endRight)
            last = std::min(last, itrRight.runEnd());
        assert(last < std::numeric_limits<int64>::max());
        
        if (itrLeft.isMarked() && itrRight.isMarked())
        {
            result.mark(currentPosition, last,
                op(itrLeft.mark(), itrRight.mark()));
        }
        else
            throw(std::logic_error("Undefined value"));
        
        currentPosition = last+1;
        
        if (itrLeft != endLeft)
            itrLeft.seek(currentPosition);
        if (itrRight != endRight)
            itrRight.seek(currentPosition);
    }
}

template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transform(const RLEBase<InSeg1> & lhs, const RLEBase<InSeg2> & rhs,
    RLEBase<OutSeg> & result, BinaryFunction op,
    const typename InSeg1::MarkType & default1,
    const typename InSeg2::MarkType & default2)
{
    result.clear();
    
    typename RLEBase<InSeg1>::ConstIterator itrLeft = lhs.unmarkedBegin();
    typename RLEBase<InSeg2>::ConstIterator itrRight = rhs.unmarkedBegin();
    typename RLEBase<InSeg1>::ConstIterator endLeft = lhs.end();
    typename RLEBase<InSeg2>::ConstIterator endRight = rhs.end();
    
    int64 currentPosition = std::min(itrLeft.position(), itrRight.position());
    
    itrLeft.seek(currentPosition);
    itrRight.seek(currentPosition);
    
    while (itrLeft != endLeft || itrRight != endRight)
    {
        int64 last = std::numeric_limits<int64>::max();
        if (itrLeft != endLeft)
            last = std::min(last, itrLeft.runEnd());
        if (itrRight != endRight)
            last = std::min(last, itrRight.runEnd());
        assert(last < std::numeric_limits<int64>::max());
        
        if (itrLeft.isMarked())
        {
            if (itrRight.isMarked())
                result.mark(currentPosition, last,
                    op(itrLeft.mark(), itrRight.mark()));
            else
                result.mark(currentPosition, last,
                    op(itrLeft.mark(), default2));
        }
        else if (itrRight.isMarked())
            result.mark(currentPosition, last, op(default1, itrRight.mark()));
        
        currentPosition = last+1;
        
        if (itrLeft != endLeft)
            itrLeft.seek(currentPosition);
        if (itrRight != endRight)
            itrRight.seek(currentPosition);
    }
}

template<class InSeg, class OutSeg, class BinaryFunction>
void transform(const RLEBase<InSeg> & lhs,
    const typename InSeg::MarkType & scalar,
    RLEBase<OutSeg> & result, BinaryFunction op)
{
    result.clear();
    
    typename RLEBase<InSeg>::ConstIterator itr;
    for (itr = lhs.begin(); itr != lhs.end(); itr.nextMarkedRun())
        result.mark(itr.runStart(), itr.runEnd(), op(itr.mark(), scalar));
}

template<class InSeg, class S, class OutSeg, class BinaryFunction>
void scalarTransform(const RLEBase<InSeg> & lhs, const S & scalar,
    RLEBase<OutSeg> & result, BinaryFunction op)
{
    result.clear();
    
    typename RLEBase<InSeg>::ConstIterator itr;
    for (itr = lhs.begin(); itr != lhs.end(); itr.nextMarkedRun())
        result.mark(itr.runStart(), itr.runEnd(), op(itr.mark(), scalar));
}

/*
template<class InSegment, class OutSegment, class UnaryFunction>
void transform(const RLEBase<InSegment> & rle, RLEBase<OutSegment> & result,
    UnaryFunction op)
{
    if (OutSegment::hasDefault)
        result.clear(op(rle.defaultMark()));
    else
        result.clear();
    
    typename RLEBase<InSegment>::ConstIterator itr;
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
        result.mark(itr.runStart(), itr.runEnd(), op(itr.mark()));
}

template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transform(const RLEBase<InSeg1> & lhs, const RLEBase<InSeg2> & rhs,
    RLEBase<OutSeg> & result, BinaryFunction op)
{
    if (OutSeg::hasDefault)
        result.clear(op(lhs.defaultMark(), rhs.defaultMark()));
    else
        result.clear();
    
    typename RLEBase<InSeg1>::ConstIterator itrLeft = lhs.unmarkedBegin();
    typename RLEBase<InSeg2>::ConstIterator itrRight = rhs.unmarkedBegin();
    typename RLEBase<InSeg1>::ConstIterator endLeft = lhs.end();
    typename RLEBase<InSeg2>::ConstIterator endRight = rhs.end();
    
    int64 currentPosition = std::min(itrLeft.position(), itrRight.position());
    
    itrLeft.seek(currentPosition);
    itrRight.seek(currentPosition);
    
    while (itrLeft != endLeft || itrRight != endRight)
    {
        int64 last = std::numeric_limits<int64>::max();
        if (itrLeft != endLeft)
            last = std::min(last, itrLeft.runEnd());
        if (itrRight != endRight)
            last = std::min(last, itrRight.runEnd());
        assert(last < std::numeric_limits<int64>::max());
        
        if (itrLeft.isMarked() || itrRight.isMarked())
        {
            result.mark(currentPosition, last,
                op(itrLeft.mark(), itrRight.mark()));
        }
        
        currentPosition = last+1;
        
        if (itrLeft != endLeft)
            itrLeft.seek(currentPosition);
        if (itrRight != endRight)
            itrRight.seek(currentPosition);
    }
}

template<class InSeg, class S, class OutSeg, class BinaryFunction>
void scalarTransform(const RLEBase<InSeg> & lhs, const S & scalar,
    RLEBase<OutSeg> & result, BinaryFunction op)
{
    if (OutSeg::hasDefault)
        result.clear(op(lhs.defaultMark(), scalar));
    else
        result.clear();
    
    typename RLEBase<InSeg>::ConstIterator itr;
    for (itr = lhs.begin(); itr != lhs.end(); itr.nextMarkedRun())
        result.mark(itr.runStart(), itr.runEnd(), op(itr.mark(), scalar));
}*/

template<class InSeg1, class InSeg2>
void restriction(const RLEBase<InSeg1> & lhs, const RLEBase<InSeg2> & rhs,
    RLEBase<InSeg1> & result)
{
    result.clear();
//    if (InSeg1::hasDefault)
//        result.clear(lhs.defaultMark());
//    else
//        result.clear();
    
    typename RLEBase<InSeg1>::ConstIterator itrLeft = lhs.unmarkedBegin();
    typename RLEBase<InSeg2>::ConstIterator itrRight = rhs.unmarkedBegin();
    typename RLEBase<InSeg1>::ConstIterator endLeft = lhs.end();
    typename RLEBase<InSeg2>::ConstIterator endRight = rhs.end();
    
    int64 currentPosition = std::min(itrLeft.position(), itrRight.position());
    
    itrLeft.seek(currentPosition);
    itrRight.seek(currentPosition);
    
    while (itrLeft != endLeft || itrRight != endRight)
    {
        int64 last = std::numeric_limits<int64>::max();
        if (itrLeft != endLeft)
            last = std::min(last, itrLeft.runEnd());
        if (itrRight != endRight)
            last = std::min(last, itrRight.runEnd());
        assert(last < std::numeric_limits<int64>::max());
        
        if (itrLeft.isMarked() && itrRight.isMarked())
        {
            InSeg1 truncated(*itrLeft, currentPosition, last);
            result.mark(currentPosition, last, truncated.data());
        }
        
        currentPosition = last+1;
        
        if (itrLeft != endLeft)
            itrLeft.seek(currentPosition);
        if (itrRight != endRight)
            itrRight.seek(currentPosition);
    }
}

template<class Segment, class UnaryPredicate>
void filter(const RLEBase<Segment> & lhs, RLEBase<Segment> & result,
    UnaryPredicate op)
{
    result.clear();
    
    typename RLEBase<Segment>::ConstIterator itr;
    
    for (itr = lhs.begin(); itr != lhs.end(); itr.nextMarkedRun())
    {
        if (op(itr.mark()))
        {
            result.mark(itr.runStart(), itr.runEnd(), itr.mark());
        }
    }
}


}; // namespace RLE

#endif
