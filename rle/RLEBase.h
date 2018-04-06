/*
 *  RLEBase.h
 *  refactorRLE
 *
 *  Created by Paul Hansen on 2/8/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _RLEBASE_
#define _RLEBASE_

#include "AVL.h"
#include <stdexcept>
#include <vector>
#include <cstdint>

namespace RLE
{
typedef std::int64_t int64;
typedef std::uint64_t uint64;

template <class Segment>
class RLEBase
{
public:
    typedef Segment SegmentType; // makes it accessible from the outside
    RLEBase();
    
    // Delete all of the runs.
    void clear();
    
    void mark(int64 first, int64 last, const typename Segment::DataType & val);
    void erase(int64 first, int64 last); // should never be called...
    bool isMarked(int64 index) const;
    typename Segment::MarkReturnType markAt(int64 index) const;
    typename Segment::MarkType markAt(int64 index,
        const typename Segment::MarkType & defaultMark) const;
    
    // Fast by-reference accessor.  Will throw an exception if index is
    // unmarked.
    typename Segment::MarkReturnType operator()(int64 index) const;
    
    // Slower by-value accessor.  Will return Segment::MarkType() if the
    // index is unmarked.  (May not always make sense!)
    typename Segment::MarkType at(int64 index) const; // will return a default
    
    uint64 bytes() const;
    uint64 numRuns() const;
    uint64 length() const; // sum lengths of all runs
    
    bool operator==(const RLEBase & rhs) const;
    bool operator!=(const RLEBase & rhs) const;
    
    typedef typename AVLTree<Segment>::Iterator AVLIterator;
    typedef typename AVLTree<Segment>::ConstIterator AVLConstIterator;
    AVLConstIterator avlBegin() const { return mSegmentTree.begin(); }
    AVLConstIterator avlEnd() const { return mSegmentTree.end(); }
    
    class Iterator;
    class ConstIterator
    {
    public:
        ConstIterator();
        ConstIterator(const RLEBase::Iterator & itr);
        ConstIterator(const RLEBase* data); // puts at start of data
        ConstIterator(const RLEBase* data, int64 unused); // puts at end of data
        
        ConstIterator operator++(int unused);
        ConstIterator & operator++();
        ConstIterator operator--(int unused);
        ConstIterator & operator--();
        void seek(int64 position);
        void move(int64 offset);
        int64 position() const;
        
        bool valid() const;
        bool isMarked() const;
        
        const Segment & operator*() const;
        const Segment* operator->() const;
        
        /*
            DynamicRLE: mark is a value (const reference)
            IndexArray: mark is an index-stride pair (const reference)
            SupportRegion: mark is a bool (by value)
            
            Thus, the traits class must have a typedef for this:
            ReturnMarkType, perhaps?
        */
        typename Segment::MarkReturnType mark() const;
        /*
            DynamicRLE: markAt is a value (const reference for run segments)
            IndexArray: markAt is an index-stride pair
            SupportRegion: markAt is a bool (by value)
        */
        typename Segment::MarkReturnType markAt(int64 position) const;
        
        // Move to the first position of the next marked or unmarked segment
        void nextRun();
        // Move to the last position of the previous marked or unmarked segment
        void previousRun();
        // Move to the first position of the next marked segment
        void nextMarkedRun();
        // Move to the last position of the previous marked segment
        void previousMarkedRun();
        
        int64 runStart() const;
        int64 runEnd() const;
        bool runIsBoundedAtStart() const;
        bool runIsBoundedAtEnd() const;
        
        bool operator==(const ConstIterator & rhs) const;
        bool operator!=(const ConstIterator & rhs) const;
        bool operator<(const ConstIterator & rhs) const;
        bool operator>(const ConstIterator & rhs) const;
        bool operator<=(const ConstIterator & rhs) const;
        bool operator>=(const ConstIterator & rhs) const;
        
    protected:
        const RLEBase* mData;
        // mConstIterator is consistent with mData->mSegmentTree.location().
        AVLConstIterator mIterator;
        int64 mCurrentPosition; // compare to mIterator.valueAtPrevious().first()
    };
    friend class ConstIterator;
    
    class Iterator
    {
    public:
        Iterator();
        Iterator(RLEBase* data); // puts at start of data
        Iterator(RLEBase* data, int64 unused); // puts at end of data
        
        Iterator operator++(int unused);
        Iterator & operator++();
        Iterator operator--(int unused);
        Iterator & operator--();
        void seek(int64 position);
        void move(int64 offset);
        int64 position() const;
        
        bool valid() const;
        bool isMarked() const;
        
        Segment & operator*() const;
        Segment* operator->() const;
        
        typename Segment::MarkReturnType mark() const;
        typename Segment::MarkReturnType markAt(int64 position) const;
        
        // Move to the first position of the next marked or unmarked segment
        void nextRun();
        // Move to the last position of the previous marked or unmarked segment
        void previousRun();
        // Move to the first position of the next marked segment
        void nextMarkedRun();
        // Move to the last position of the previous marked segment
        void previousMarkedRun();
        
        int64 runStart() const;
        int64 runEnd() const;
        bool runIsBoundedAtStart() const;
        bool runIsBoundedAtEnd() const;
        
        bool operator==(const Iterator & rhs) const;
        bool operator!=(const Iterator & rhs) const;
        bool operator<(const Iterator & rhs) const;
        bool operator>(const Iterator & rhs) const;
        bool operator<=(const Iterator & rhs) const;
        bool operator>=(const Iterator & rhs) const;
        
        friend class ConstIterator;
    protected:
        RLEBase* mData;
        // mIterator is consistent with mData->mSegmentTree.location().
        AVLIterator mIterator;
        int64 mCurrentPosition; // compare to mIterator.valueAtPrevious().first()
    };
    friend class Iterator;
    
    // For iterating over all marked segments.
    Iterator begin();
    Iterator end();
    ConstIterator begin() const;
    ConstIterator end() const;
    
    // For iterating over all segments INCLUDING the unmarked ones.
    Iterator unmarkedBegin();
    ConstIterator unmarkedBegin() const;
    
    // The fact is that unmarkedEnd() is not well defined.  To iterate
    // backwards over all runs, you should start at end(), I guess... ick.
    //Iterator unmarkedLast(); // can't get to the end
    
    Iterator location(int64 position);
    ConstIterator location(int64 position) const;
    
    // Jump to a given MARKED segment.  This is for debugging convenience and
    // nothing else—it's a nasty, slow serial search.
    ConstIterator segment(int64 which) const;
    
protected:
    AVLTree<Segment> mSegmentTree;
};


template<class Segment>
std::ostream & operator<<(std::ostream & str,
    const RLEBase<Segment> & rle)
{
    typename RLEBase<Segment>::ConstIterator itr;
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        str << "[" << itr->first() << ", " << itr->last() << "]: "
            << itr->mark() << "\t";
    }
    return str;
}

class BaseSegment
{
public:
    // Concept requirements:
    // typedef DataType
    // typedef DataReturnType (probably DataType or DataType&)
    // typedef MarkType (possibly same as DataType)
    // typedef MarkReturnType (probably MarkType or MarkType&)
    // enum { hasDefault };
    
    BaseSegment(int64 first, int64 last) : mFirst(first), mLast(last) {}
    
    int64 first() const { return mFirst; }
    int64 last() const { return mLast; }
    int64 length() const { return mLast - mFirst + 1; }
    
    // chopped up pieces and such
    void chopTail(int64 last) { mLast = last; }
    void chopHead(int64 first) { mFirst = first; }
    
    // extent predicates
    bool leftTouches(const BaseSegment & segment) const
    {
        return (mFirst >= segment.first() && mFirst <= segment.last()+1);
    }
    bool rightTouches(int64 first, int64 last) const
    {
        return (mLast >= first-1 && mLast <= last);
    }
    bool rightTouches(const BaseSegment & segment) const
    {
        return (mLast >= segment.first()-1 && mLast <= segment.last());
    }
    bool leftOverlaps(int64 first, int64 last) const
    {
        return mFirst <= last;
    }
    bool leftOverlaps(const BaseSegment & segment) const
    {
        return mFirst <= segment.last();
    }
    bool rightOverlaps(int64 first, int64 last) const
    {
        return mLast >= first;
    }
    bool rightOverlaps(const BaseSegment & segment) const
    {
        return mLast >= segment.first();
    }
    bool encloses(int64 first, int64 last) const
    {
        return (mFirst <= first && mLast >= last);
    }
    bool encloses(const BaseSegment & segment) const
    {
        return mFirst <= segment.first() && mLast >= segment.last();
    }
    bool surrounds(int64 first, int64 last) const
    {
        return (mFirst < first && mLast > last);
    }
    bool surrounds(const BaseSegment & segment) const
    {
        return (mFirst < segment.first() && mLast > segment.last());
    }
    bool inside(int64 first, int64 last) const
    {
        return mFirst >= first && mLast <= last;
    }
protected:
    int64 mFirst;
    int64 mLast;
};

inline bool operator<(const BaseSegment & lhs, const BaseSegment & rhs)
{
    return lhs.first() < rhs.first();
}

inline bool operator<(const BaseSegment & lhs, int64 rhs)
{
    return lhs.first() < rhs;
}

inline bool operator<(int64 lhs, const BaseSegment & rhs)
{
    return lhs < rhs.first();
}


// Translation
template<class Segment>
void translate(const RLEBase<Segment> & rle, RLEBase<Segment> & result,
    int64 distance);


// -------------------------------------
// Unary transformations

// outRLE = op(RLE)
//template<class RLE, class UnaryFunction>
//RLEBase<InSegment> transform(const RLEBase<InSegment> & rle, UnaryFunction op);

// op(RLE, outRLE)
template<class InSegment, class OutSegment, class UnaryFunction>
void transform(const RLEBase<InSegment> & rle, RLEBase<OutSegment> & result,
    UnaryFunction op);

// -------------------------------------
// Binary transformations

// outRLE = op(RLE, RLE)
//template<class InSegment, class BinaryFunction>
//RLEBase<InSegment> transform(const RLEBase<InSegment> & in1,
//    const RLEBase<InSegment> & in2, BinaryFunction op);

// op(RLE, RLE, outRLE)
template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transform(const RLEBase<InSeg1> & in1, const RLEBase<InSeg2> & in2,
    RLEBase<OutSeg> & result, BinaryFunction op);

// op(RLE, RLE, outRLE) with default marks
template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transform(const RLEBase<InSeg1> & in1, const RLEBase<InSeg2> & in2,
    RLEBase<OutSeg> & result, BinaryFunction op,
    const typename InSeg1::MarkType & default1,
    const typename InSeg2::MarkType & default2);

// outRLE = op(RLE, scalar)
//template<class InSegment, class BinaryFunction>
//RLEBase<InSegment> transform(const RLEBase<InSegment> & in1,
//    const InSegment::MarkType & scalar, BinaryFunction op);

// op(RLE, scalar, outRLE)
template<class InSegment, class OutSegment, class BinaryFunction>
void transform(const RLEBase<InSegment> & in1,
    const typename InSegment::MarkType & scalar, RLEBase<OutSegment> & result,
    BinaryFunction op);

// op(RLE, scalar, outRLE) for ANY scalar type T
// This function needs a weird name to avoid template ambiguity.
template<class InSegment, class T, class OutSegment, class BinaryFunction>
void scalarTransform(const RLEBase<InSegment> & in1,
    const T & scalar, RLEBase<OutSegment> & result, BinaryFunction op);

// Special type of binary transformation
template<class InSeg1, class InSeg2>
void restriction(const RLEBase<InSeg1> & lhs, const RLEBase<InSeg2> & rhs,
    RLEBase<InSeg1> & result);

// Quick implementation that won't work for runs of varying data like IndexArray
template<class Segment, class UnaryPredicate>
void filter(const RLEBase<Segment> & lhs, RLEBase<Segment> & result,
    UnaryPredicate op);



}; // namespace RLE

#include "RLEBase-inl.h"

#endif
