/*
 *  RLEBase3.h
 *  refactorRLE
 *
 *  Created by Paul Hansen on 2/9/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _RLEBASE3_
#define _RLEBASE3_

#include "RLEBase.h"

namespace RLE
{

template <class Segment>
class RLEBase3 : public RLEBase<Segment>
{
public:
    // TODO: change language: symmetry
    RLEBase3(bool hasX = true, bool hasY = true, bool hasZ = true);
    
    RLEBase3(const RLEBase3<Segment> & copyMe);
    
    // TODO: symmetry
    template<class Scalar>
    RLEBase3(const Scalar nonSymmetricDimensions[]);

// This is too vague, I can't enforce that dims is an array.    
//    template<class Array>
//    RLEBase3(const Array & dims); // must have operator[]
    
    void mark(int64 x0, int64 y0, int64 z0, int64 x1, int64 y1, int64 z1,
        const typename Segment::DataType & val);
    void mark(int64 x, int64 y, int64 z, const typename Segment::DataType & val);
    void erase(int64 x0, int64 y0, int64 z0, int64 x1, int64 y1, int64 z1);
    void erase(int64 x, int64 y, int64 z);
    bool isMarked(int64 x0, int64 y0, int64 z0) const;
    typename Segment::MarkReturnType markAt(int64 ii, int64 jj, int64 kk) const;
    typename Segment::MarkType markAt(int64 ii, int64 jj, int64 kk,
        const typename Segment::MarkType & defaultMark) const;
    
    // Fast by-reference accessor.  Will throw an exception if index is
    // unmarked.
    typename Segment::MarkReturnType operator()(int64 ii, int64 jj, int64 kk)
        const;
    
    // Slower by-value accessor.  Will return Segment::MarkType() if the
    // index is unmarked.  (May not always make sense!)
    typename Segment::MarkType at(int64 ii, int64 jj, int64 kk) const;
    
    // TODO: think; symmetry
    // Set the capacities and strides.  This function does not change non-empty
    // arrays!
    void dimensions(bool hasX, bool hasY, bool hasZ);
    
    // TODO: think; symmetry
//    template<class Scalar>
//    void dimensions(const Scalar size[]);
    void dimensions(const bool nonSymmetricDimensions[]);
    
    // TODO: think; symmetry
    bool hasDimension(int nn) const { return mNonSymmetricDimensions[nn]; }
    //bool hasDimension(int nn) const { return mStride[nn] != 0; }
    
    // TODO: think; symmetry
    // Change the dimension.
    void changeDimensions(bool xx, bool yy, bool zz);
    template<class Scalar>
    void changeDimensions(const Scalar nonSymmetricDimensions[]);
    
    // TODO: keep capacity().
    const int64* capacity() const { return mCapacity; }
    const int64* minCoord() const { return mMinCoord; }
    const int64* stride() const { return mStride; }
    const bool* nonSymmetricDimensions() const { return mNonSymmetricDimensions; }
    int64 capacity(int nn) const { return mCapacity[nn]; }
    int64 minCoord(int nn) const { return mMinCoord[nn]; }
    int64 stride(int nn) const { return mStride[nn]; }
    
    int64 axis() const { return mAxis; }
    
    // TODO: symmetry
    int64 numDimensions() const;
    
    // TODO: subscript/index language?  consider it.
    int64 linearIndex(int64 ii, int64 jj, int64 kk) const;
    
    template<class Scalar>
    void cartesianCoordinates(int64 index, Scalar & ii, Scalar & jj, Scalar & kk) const;
    
    bool operator==(const RLEBase3<Segment> & rhs) const;
    bool operator!=(const RLEBase3<Segment> & rhs) const;
    
    typename RLEBase<Segment>::Iterator location(int64 ii, int64 jj, int64 kk);
    typename RLEBase<Segment>::ConstIterator location(int64 ii, int64 jj, int64 kk)
        const;
    
protected:
    int64 mCapacity[3];
    int64 mMinCoord[3];
    int64 mStride[3];
    bool mNonSymmetricDimensions[3];
    
    int64 mAxis;
};


template<class Segment>
std::ostream & operator<<(std::ostream & str,
    const RLEBase3<Segment> & rle)
{
    int64 ii, jj, kk;
    typename RLEBase<Segment>::ConstIterator itr;
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        rle.cartesianCoordinates(itr->first(), ii, jj, kk);
        str << "[(" << ii << ", " << jj << ", " << kk << " ), (";
        rle.cartesianCoordinates(itr->last(), ii, jj, kk);
        str << ii << ", " << jj << ", " << kk << ")]: "
            << itr->mark() << "\t";
    }
    return str;
}

// Translation
template<class Segment, class Scalar>
void translate(const RLEBase3<Segment> & rle, RLEBase3<Segment> & result,
    Scalar dx, Scalar dy, Scalar dz);

template<class Segment, class Scalar>
void translate(const RLEBase3<Segment> & rle, RLEBase3<Segment> & result,
    const Scalar* dxyz);

// -------------------------------------
// Unary transformations

// op(RLE, outRLE)
template<class InSegment, class OutSegment, class UnaryFunction>
void transform(const RLEBase3<InSegment> & rle, RLEBase3<OutSegment> & result,
    UnaryFunction op);
    
// -------------------------------------
// Binary transformations

// op(RLE, RLE, outRLE)
template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transform(const RLEBase3<InSeg1> & in1, const RLEBase3<InSeg2> & in2,
    RLEBase3<OutSeg> & result, BinaryFunction op);

// op(RLE, RLE, outRLE) with default marks
template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transform(const RLEBase3<InSeg1> & in1, const RLEBase3<InSeg2> & in2,
    RLEBase3<OutSeg> & result, BinaryFunction op,
    const typename InSeg1::MarkType & default1,
    const typename InSeg2::MarkType & default2);

// op(RLE, scalar, outRLE)
template<class InSegment, class OutSegment, class BinaryFunction>
void transform(const RLEBase3<InSegment> & in1,
    const typename InSegment::MarkType & scalar, RLEBase3<OutSegment> & result,
    BinaryFunction op);

// op(RLE, scalar, outRLE) for ANY scalar type T
// This function needs a weird name to avoid template ambiguity.
template<class InSegment, class T, class OutSegment, class BinaryFunction>
void scalarTransform(const RLEBase3<InSegment> & in1,
    const T & scalar, RLEBase3<OutSegment> & result, BinaryFunction op);

// TODO: figure out how to handle arbitrary symmetries e.g. X-symm x Y-symm.
// lhs is of higher dimensionality than rhs.
template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transformHighLow(const RLEBase3<InSeg1> & lhs,
    const RLEBase3<InSeg2> & rhs, RLEBase3<OutSeg> & result, BinaryFunction op);

// lhs is of higher dimensionality than rhs.
template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transformHighLow(const RLEBase3<InSeg1> & lhs,
    const RLEBase3<InSeg2> & rhs, RLEBase3<OutSeg> & result, BinaryFunction op,
    const typename InSeg1::MarkType & default1,
    const typename InSeg2::MarkType & default2);

/*
// Unary transformation, acting on the mark/data type
template<class InSegment, class OutSegment, class UnaryFunction>
void transform(const RLEBase3<InSegment> & rle, RLEBase3<OutSegment> & result,
    UnaryFunction op);

// Binary transformation (RLE, RLE), acting on the mark/data type
template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transform(const RLEBase3<InSeg1> & lhs, const RLEBase3<InSeg2> & rhs,
    RLEBase3<OutSeg> & result, BinaryFunction op);

// Binary transformation (RLE, RLE), acting on the mark/data type
// lhs is of higher dimensionality than rhs.
template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transformHighLow(const RLEBase3<InSeg1> & lhs,
    const RLEBase3<InSeg2> & rhs, RLEBase3<OutSeg> & result, BinaryFunction op);

// Binary transformation (RLE, scalar), acting on the mark/data type
template<class InSeg1, class S, class OutSeg, class BinaryFunction>
void scalarTransform(const RLEBase3<InSeg1> & lhs, const S & scalar,
    RLEBase3<OutSeg> & result, BinaryFunction op);
*/

// Special type of binary transformation
template<class InSeg1, class InSeg2>
void restriction(const RLEBase3<InSeg1> & lhs, const RLEBase3<InSeg2> & rhs,
    RLEBase3<InSeg1> & result);

template<class InSeg1, class InSeg2>
void restrictionHighLow(const RLEBase3<InSeg1> & lhs,
    const RLEBase3<InSeg2> & rhs, RLEBase3<InSeg1> & result);

template<class InSeg1, class InSeg2>
void restrictionLowHigh(const RLEBase3<InSeg1> & lhs,
    const RLEBase3<InSeg2> & rhs, RLEBase3<InSeg1> & result);

template<class Segment, class UnaryPredicate>
void filter(const RLEBase3<Segment> & lhs, RLEBase3<Segment> & result,
    UnaryPredicate op);


}; // namespace RLE


#include "RLEBase3-inl.h"


#endif

