/*
 *  RLEBase3-inl.h
 *  refactorRLE
 *
 *  Created by Paul Hansen on 2/9/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _RLEBASE3_

#include <cmath>

namespace RLE
{

template<class Segment>
RLEBase3<Segment>::
RLEBase3(bool hasX, bool hasY, bool hasZ) :
    RLEBase<Segment>()
{
    dimensions(hasX, hasY, hasZ);
}

template<class Segment>
RLEBase3<Segment>::
RLEBase3(const RLEBase3<Segment> & copyMe) :
    RLEBase<Segment>()
{
//    dimensions(copyMe.capacity());
    *this = copyMe;
}

template<class Segment>
template<class Scalar>
RLEBase3<Segment>::
RLEBase3(const Scalar nonSymmetricDimensions[]) :
    RLEBase<Segment>()
{
    dimensions(nonSymmetricDimensions[0], nonSymmetricDimensions[1], nonSymmetricDimensions[2]);
}


//template<class Segment>
//template<class Array>
//RLEBase3<Segment>::
//RLEBase3(const Array & dims) :
//    RLEBase<Segment>()
//{
//    dimensions(dims[0] > 1, dims[1] > 1, dims[2] > 1);
//}

template<class Segment>
void RLEBase3<Segment>::
dimensions(bool hasX, bool hasY, bool hasZ)
{
    if (RLEBase<Segment>::numRuns() != 0)
        throw(std::logic_error("Setting dimensions of non-empty array"));
    
    mNonSymmetricDimensions[0] = hasX;
    mNonSymmetricDimensions[1] = hasY;
    mNonSymmetricDimensions[2] = hasZ;
    
    // try to find an equitable division of memory.
    int64 numDims = int64(hasX) + int64(hasY) + int64(hasZ);
    if (numDims == 0)
    {
//        hasX = true;
//        numDims = 1;
        
        // The weird case!!
        for (int ii = 0; ii < 3; ii++)
        {
            mCapacity[ii] = 1;
            mMinCoord[ii] = 0;
            mStride[ii] = 0;
            mAxis = 0;
        }
        return;
    }
    assert(numDims > 0);
    assert(numDims <= 3);
    
    int64 int64Bits = log2(std::numeric_limits<int64>::max());
    assert(int64Bits > 0);
    int64 eachCapacity = (1 << int64Bits/numDims) - 1;
    
    mCapacity[0] = hasX ? eachCapacity : 1;
    mCapacity[1] = hasY ? eachCapacity : 1;
    mCapacity[2] = hasZ ? eachCapacity : 1;
    
    mMinCoord[0] = -mCapacity[0]/2;
    mMinCoord[1] = -mCapacity[1]/2;
    mMinCoord[2] = -mCapacity[2]/2;
    
    mStride[0] = hasX*1;
    mStride[1] = hasY*mCapacity[0];
    mStride[2] = hasZ*mCapacity[0] * mCapacity[1];
    
    if (mStride[0] == 0 && mStride[1] == 0)
        mAxis = 2;
    else if (mStride[0] == 0)
        mAxis = 1;
    else
        mAxis = 0;
}

template<class Segment>
//template<class Scalar>
void RLEBase3<Segment>::
dimensions(const bool nonSymmetricDimensions[])
{
    dimensions(nonSymmetricDimensions[0], nonSymmetricDimensions[1], nonSymmetricDimensions[2]);
}

template<class Segment>
void RLEBase3<Segment>::
changeDimensions(bool xx, bool yy, bool zz)
{
    bool xyz[] = {xx, yy, zz};
    changeDimensions(xyz);
}

template<class Segment>
template<class Scalar>
void RLEBase3<Segment>::
changeDimensions(const Scalar nonSymmetricDimensions[])
{
    RLEBase3<Segment> tmp = *this;
    RLEBase<Segment>::clear();
    dimensions(nonSymmetricDimensions[0], nonSymmetricDimensions[1], nonSymmetricDimensions[2]);
    
    typename RLEBase<Segment>::ConstIterator itr;
    for (itr = tmp.begin(); itr != tmp.end(); itr.nextMarkedRun())
    {
        int64 x0, y0, z0, x1, y1, z1;
        tmp.cartesianCoordinates(itr.runStart(), x0, y0, z0);
        tmp.cartesianCoordinates(itr.runEnd(), x1, y1, z1);
        
        mark(x0, y0, z0, x1, y1, z1, itr->data());
    }
}

template<class Segment>
int64 RLEBase3<Segment>::
numDimensions() const
{
    return (mNonSymmetricDimensions[0] != 0) + (mNonSymmetricDimensions[1] != 0) + (mNonSymmetricDimensions[2] != 0);
}

template<class Segment>
void RLEBase3<Segment>::
mark(int64 x0, int64 y0, int64 z0, int64 x1, int64 y1, int64 z1,
    const typename Segment::DataType & val)
{
    if (mAxis == 0)
    {
        for (int64 zz = z0; zz <= z1; zz++)
        for (int64 yy = y0; yy <= y1; yy++)
        {
            RLEBase<Segment>::mark(linearIndex(x0, yy, zz),
                linearIndex(x1, yy, zz), val);
        }
    }
    else if (mAxis == 1)
    {
        for (int64 zz = z0; zz <= z1; zz++)
        {
            RLEBase<Segment>::mark(linearIndex(0, y0, zz),
                linearIndex(0, y1, zz), val);
        }
    }
    else if (mAxis == 2)
    {
        RLEBase<Segment>::mark(linearIndex(0, 0, z0),
            linearIndex(0, 0, z1), val);
    }
}

template<class Segment>
void RLEBase3<Segment>::
mark(int64 x, int64 y, int64 z, const typename Segment::DataType & val)
{
    mark(x, y, z, x, y, z, val);
}

template<class Segment>
void RLEBase3<Segment>::
erase(int64 x0, int64 y0, int64 z0, int64 x1, int64 y1, int64 z1)
{
    if (mAxis == 0)
    {
        for (int64 zz = z0; zz <= z1; zz++)
        for (int64 yy = y0; yy <= y1; yy++)
        {
            RLEBase<Segment>::erase(linearIndex(x0, yy, zz),
                linearIndex(x1, yy, zz));
        }
    }
    else if (mAxis == 1)
    {
        for (int64 zz = z0; zz <= z1; zz++)
        {
            RLEBase<Segment>::erase(linearIndex(0, y0, zz),
                linearIndex(0, y1, zz));
        }
    }
    else if (mAxis == 2)
    {
        RLEBase<Segment>::erase(linearIndex(0, 0, z0), linearIndex(0, 0, z1));
    }
}

template<class Segment>
void RLEBase3<Segment>::
erase(int64 x, int64 y, int64 z)
{
    erase(x, y, z, x, y, z);
}

template<class Segment>
bool RLEBase3<Segment>::
isMarked(int64 x0, int64 y0, int64 z0) const
{
    return RLEBase<Segment>::isMarked(linearIndex(x0, y0, z0));
}

template<class Segment>
typename Segment::MarkReturnType RLEBase3<Segment>::
markAt(int64 ii, int64 jj, int64 kk) const
{
    return RLEBase<Segment>::markAt(linearIndex(ii, jj, kk));
}

template<class Segment>
typename Segment::MarkType RLEBase3<Segment>::
markAt(int64 ii, int64 jj, int64 kk,
    const typename Segment::MarkType & defaultMark) const
{
    return RLEBase<Segment>::markAt(linearIndex(ii, jj, kk), defaultMark);
}

template<class Segment>
typename Segment::MarkReturnType RLEBase3<Segment>::
operator()(int64 ii, int64 jj, int64 kk) const
{
    return RLEBase<Segment>::operator()(linearIndex(ii, jj, kk));
}

template<class Segment>
typename Segment::MarkType RLEBase3<Segment>::
at(int64 ii, int64 jj, int64 kk) const
{
    return RLEBase<Segment>::at(linearIndex(ii, jj, kk));
}


template<class Segment>
int64 RLEBase3<Segment>::
linearIndex(int64 ii, int64 jj, int64 kk) const
{
    ii -= mMinCoord[0];
    jj -= mMinCoord[1];
    kk -= mMinCoord[2];

    // Do we need to check these?  Hrrm.  Depends on what I do with out of
    // bounds indices int64ernally... this has to do with 3D indices int64o 2D
    // arrays, which I think should be ok (the coordinate gets projected down).
//    assert(ii >= 0);
//    assert(jj >= 0);
//    assert(kk >= 0);
    
    return ii*mStride[0] + jj*mStride[1] + kk*mStride[2];
}

template<class Segment>
template<class Scalar>
void RLEBase3<Segment>::
cartesianCoordinates(int64 index, Scalar & ii, Scalar & jj, Scalar & kk) const
{
    if (mStride[2])
        kk = index / mStride[2];
    else
        kk = 0;
    index -= kk*mStride[2];
    
    if (mStride[1])
        jj = index / mStride[1];
    else
        jj = 0;
    index -= jj*mStride[1];
    
    if (mStride[0])
        ii = index;
    else
        ii = 0;
    
    ii += mMinCoord[0];
    jj += mMinCoord[1];
    kk += mMinCoord[2];
}

template<class Segment>
bool RLEBase3<Segment>::
operator==(const RLEBase3<Segment> & rhs) const
{
    return mCapacity[0] == rhs.mCapacity[0] && mCapacity[1] == rhs.mCapacity[1]
        && mCapacity[2] == mCapacity[2] && RLEBase<Segment>::operator==(rhs);
}

template<class Segment>
bool RLEBase3<Segment>::
operator!=(const RLEBase3<Segment> & rhs) const
{
    return !(*this == rhs);
}

template<class Segment>
typename RLEBase<Segment>::Iterator RLEBase3<Segment>::
location(int64 ii, int64 jj, int64 kk)
{
    return RLEBase<Segment>::location(linearIndex(ii,jj,kk));
}

template<class Segment>
typename RLEBase<Segment>::ConstIterator RLEBase3<Segment>::
location(int64 ii, int64 jj, int64 kk) const
{
    return RLEBase<Segment>::location(linearIndex(ii,jj,kk));
}



#pragma mark *** Operations ***



template<class Segment, class Scalar>
void translate(const RLEBase3<Segment> & rle, RLEBase3<Segment> & result,
    Scalar dx, Scalar dy, Scalar dz)
{
    result.dimensions(rle.nonSymmetricDimensions());
        
    translate( (const RLEBase<Segment>&)rle, (RLEBase<Segment>&)result, 
        dx*rle.stride(0) + dy*rle.stride(1) + dz*rle.stride(2));
}


template<class Segment, class Scalar>
void translate(const RLEBase3<Segment> & rle, RLEBase3<Segment> & result,
    const Scalar* dxyz)
{
    result.dimensions(rle.nonSymmetricDimensions());
        
    translate( (const RLEBase<Segment>&)rle, (RLEBase<Segment>&)result, 
        dxyz[0]*rle.stride(0) + dxyz[1]*rle.stride(1) + dxyz[2]*rle.stride(2));
}




template<class InSegment, class OutSegment, class UnaryFunction>
void transform(const RLEBase3<InSegment> & rle, RLEBase3<OutSegment> & result,
    UnaryFunction op)
{
    result.dimensions(rle.nonSymmetricDimensions());
    transform((const RLEBase<InSegment>&)rle, (RLEBase<OutSegment>&)result, op);
}

template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transform(const RLEBase3<InSeg1> & lhs, const RLEBase3<InSeg2> & rhs,
    RLEBase3<OutSeg> & result, BinaryFunction op)
{
    if (lhs.capacity(0) != rhs.capacity(0))
        throw(std::logic_error("Mismatched capacity X"));
    if (lhs.capacity(1) != rhs.capacity(1))
        throw(std::logic_error("Mismatched capacity Y"));
    if (lhs.capacity(2) != rhs.capacity(2))
        throw(std::logic_error("Mismatched capacity Z"));
    
    result.dimensions(lhs.nonSymmetricDimensions());
    RLE::transform( (const RLEBase<InSeg1>&)lhs,
        (const RLEBase<InSeg2>&)rhs, (RLEBase<OutSeg>&) result, op);
}

template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transform(const RLEBase3<InSeg1> & lhs, const RLEBase3<InSeg2> & rhs,
    RLEBase3<OutSeg> & result, BinaryFunction op,
    const typename InSeg1::MarkType & default1,
    const typename InSeg2::MarkType & default2)
{
    if (lhs.capacity(0) != rhs.capacity(0))
        throw(std::logic_error("Mismatched capacity X"));
    if (lhs.capacity(1) != rhs.capacity(1))
        throw(std::logic_error("Mismatched capacity Y"));
    if (lhs.capacity(2) != rhs.capacity(2))
        throw(std::logic_error("Mismatched capacity Z"));
    
    result.dimensions(lhs.nonSymmetricDimensions());
    RLE::transform( (const RLEBase<InSeg1>&)lhs,
        (const RLEBase<InSeg2>&)rhs, (RLEBase<OutSeg>&) result, op,
        default1, default2);
}

template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transformHighLow(const RLEBase3<InSeg1> & lhs,
    const RLEBase3<InSeg2> & rhs, RLEBase3<OutSeg> & result, BinaryFunction op)
{
    assert(lhs.numDimensions() > rhs.numDimensions());
    assert( (lhs.stride(0) && rhs.stride(0)) ||
        (lhs.stride(1) && rhs.stride(1)) ||
        (lhs.stride(2) && rhs.stride(2)) );
    
    result.clear();
    
    result.dimensions(lhs.nonSymmetricDimensions());
    
    typename RLEBase<InSeg1>::ConstIterator itrLeft(lhs.unmarkedBegin());
    typename RLEBase<InSeg1>::ConstIterator endLeft(lhs.end());
    typename RLEBase<InSeg2>::ConstIterator itrRight(rhs.unmarkedBegin());
    typename RLEBase<InSeg2>::ConstIterator endRight(rhs.end());
    
    int64 indexLeft, indexRight, curX, curY, curZ;
    
    indexLeft = itrLeft.position();
    lhs.cartesianCoordinates(indexLeft, curX, curY, curZ);
    indexRight = rhs.linearIndex(curX, curY, curZ);
    
    itrLeft.seek(indexLeft);
    itrRight.seek(indexRight);
    
    while (itrLeft != endLeft)
    {
        int64 length = itrLeft.runEnd() - indexLeft + 1;
        if (rhs.axis() == lhs.axis() && itrRight < endRight)
        {
            // The length of the rhs runline is:
            //  infinity    if rhs.axis() != lhs.axis() (extrude!)
            //  infinity    if itrRight >= endRight (unbounded last runline)
            //  something   otherwise
            length = std::min(length, itrRight.runEnd() - indexRight + 1);
        }
        
        if (itrLeft.isMarked() && itrRight.isMarked())
        {
            result.RLEBase<OutSeg>::mark(indexLeft, indexLeft + length - 1,
                op(itrLeft.mark(), itrRight.mark()));
        }
        else
            throw(std::logic_error("Undefined value"));
        
        indexLeft += length;
        lhs.cartesianCoordinates(indexLeft, curX, curY, curZ);
        indexRight = rhs.linearIndex(curX, curY, curZ);
        
        if (itrLeft != endLeft)
            itrLeft.seek(indexLeft);
        itrRight.seek(indexRight);
    }    
}

template<class InSeg1, class InSeg2, class OutSeg, class BinaryFunction>
void transformHighLow(const RLEBase3<InSeg1> & lhs,
    const RLEBase3<InSeg2> & rhs, RLEBase3<OutSeg> & result, BinaryFunction op,
    const typename InSeg1::MarkType & default1,
    const typename InSeg2::MarkType & default2)
{
    assert(lhs.numDimensions() > rhs.numDimensions());
    assert( (lhs.stride(0) && rhs.stride(0)) ||
        (lhs.stride(1) && rhs.stride(1)) ||
        (lhs.stride(2) && rhs.stride(2)) );
    
    result.clear();
    result.dimensions(lhs.nonSymmetricDimensions());
    
    typename RLEBase<InSeg1>::ConstIterator itrLeft(lhs.unmarkedBegin());
    typename RLEBase<InSeg1>::ConstIterator endLeft(lhs.end());
    typename RLEBase<InSeg2>::ConstIterator itrRight(rhs.unmarkedBegin());
    typename RLEBase<InSeg2>::ConstIterator endRight(rhs.end());
    
    int64 indexLeft, indexRight, curX, curY, curZ;
    
    indexLeft = itrLeft.position();
    lhs.cartesianCoordinates(indexLeft, curX, curY, curZ);
    indexRight = rhs.linearIndex(curX, curY, curZ);
    
    itrLeft.seek(indexLeft);
    itrRight.seek(indexRight);
    
    while (itrLeft != endLeft)
    {
        uint64 length = itrLeft.runEnd() - indexLeft + 1;
        if (rhs.axis() == lhs.axis() && itrRight < endRight)
        {
            // The length of the rhs runline is:
            //  infinity    if rhs.axis() != lhs.axis() (extrude!)
            //  infinity    if itrRight >= endRight (unbounded last runline)
            //  something   otherwise
            length = std::min((int64) length, itrRight.runEnd() - indexRight + 1);
        }
        
        if (itrLeft.isMarked())
        {
            if (itrRight.isMarked())
            {
                result.RLEBase<OutSeg>::mark(indexLeft, indexLeft + length - 1,
                    op(itrLeft.mark(), itrRight.mark()));
            }
            else
            {
                result.RLEBase<OutSeg>::mark(indexLeft, indexLeft + length - 1,
                    op(itrLeft.mark(), default2));
            }
        }
        else if (itrRight.isMarked())
        {
            result.RLEBase<OutSeg>::mark(indexLeft, indexLeft + length - 1,
                op(default1, itrRight.mark()));
        }
        
        indexLeft += length;
        lhs.cartesianCoordinates(indexLeft, curX, curY, curZ);
        indexRight = rhs.linearIndex(curX, curY, curZ);
        
        if (itrLeft != endLeft)
            itrLeft.seek(indexLeft);
        itrRight.seek(indexRight);
    }
}

template<class InSeg, class OutSeg, class BinaryFunction>
void transform(const RLEBase3<InSeg> & lhs,
    const typename InSeg::MarkType & scalar,
    RLEBase3<OutSeg> & result, BinaryFunction op)
{
    result.dimensions(lhs.nonSymmetricDimensions());
    transform( (const RLEBase<InSeg>&)lhs, scalar,
        (RLEBase<OutSeg>&) result, op);
}

template<class InSeg, class S, class OutSeg, class BinaryFunction>
void scalarTransform(const RLEBase3<InSeg> & lhs, const S & scalar,
    RLEBase3<OutSeg> & result, BinaryFunction op)
{
    result.dimensions(lhs.nonSymmetricDimensions());
    scalarTransform( (const RLEBase<InSeg>&)lhs, scalar,
        (RLEBase<OutSeg>&) result, op);
}

template<class InSeg1, class InSeg2>
void restriction(const RLEBase3<InSeg1> & lhs, const RLEBase3<InSeg2> & rhs,
    RLEBase3<InSeg1> & result)
{
    if (lhs.numDimensions() > rhs.numDimensions())
        restrictionHighLow(lhs, rhs, result);
    else if (lhs.numDimensions() < rhs.numDimensions())
        restrictionLowHigh(lhs, rhs, result);
    else
    {
        assert(lhs.capacity(0) == rhs.capacity(0) &&
            lhs.capacity(1) == rhs.capacity(1) &&
            lhs.capacity(2) == rhs.capacity(2));
        
        result.dimensions(lhs.nonSymmetricDimensions());
        
        restriction((const RLEBase<InSeg1>&)lhs, (const RLEBase<InSeg2>&) rhs,
            (RLEBase<InSeg1>&)result);
    }
}

template<class InSeg1, class InSeg2>
void restrictionHighLow(const RLEBase3<InSeg1> & lhs,
    const RLEBase3<InSeg2> & rhs, RLEBase3<InSeg1> & result)
{
    assert(lhs.numDimensions() > rhs.numDimensions());
    // check that lhs and rhs share at least one non-symmetrical dimension.
//    std::cout << "LHS stride " << lhs.stride(0) << " " << lhs.stride(1)
//        << " " << lhs.stride(2) << "\n";
//    std::cout << "RHS stride " << rhs.stride(0) << " " << rhs.stride(1)
//        << " " << rhs.stride(2) << "\n";
    assert( (lhs.stride(0) && rhs.stride(0)) ||
        (lhs.stride(1) && rhs.stride(1)) ||
        (lhs.stride(2) && rhs.stride(2)) );
    result.dimensions(lhs.nonSymmetricDimensions());
    result.clear();
    
    typename RLEBase<InSeg1>::ConstIterator itrLeft(lhs.unmarkedBegin());
    typename RLEBase<InSeg1>::ConstIterator endLeft(lhs.end());
    typename RLEBase<InSeg2>::ConstIterator itrRight(rhs.unmarkedBegin());
    typename RLEBase<InSeg2>::ConstIterator endRight(rhs.end());
    
    int64 indexLeft, indexRight, curX, curY, curZ;
    
    indexLeft = itrLeft.position();
    lhs.cartesianCoordinates(indexLeft, curX, curY, curZ);
    indexRight = rhs.linearIndex(curX, curY, curZ);
    
//    std::cout << "Start:\n";
//    std::cout << "L " << indexLeft << " xyz " << curX << " " << curY
//        << " " << curZ << " R " << indexRight << "\n";
    
    itrLeft.seek(indexLeft);
    itrRight.seek(indexRight);
    
    while (itrLeft != endLeft)
    {
        int64 length = itrLeft.runEnd() - indexLeft + 1;
        if (rhs.axis() == lhs.axis() && itrRight < endRight)
        {
            // The length of the rhs runline is:
            //  infinity    if rhs.axis() != lhs.axis() (extrude!)
            //  infinity    if itrRight >= endRight (unbounded last runline)
            //  something   otherwise
            length = std::min(length, itrRight.runEnd() - indexRight + 1);
        }
        
        if (itrLeft.isMarked() && itrRight.isMarked())
        {
            InSeg1 truncated(*itrLeft, indexLeft, indexLeft + length - 1);
            result.RLEBase<InSeg1>::mark(truncated.first(),
                truncated.last(), truncated.data());
        }
        
        indexLeft += length;
        lhs.cartesianCoordinates(indexLeft, curX, curY, curZ);
        indexRight = rhs.linearIndex(curX, curY, curZ);
        
//        std::cout << "L " << indexLeft << " xyz " << curX << " " << curY
//            << " " << curZ << " R " << indexRight << "\n";
        
        if (itrLeft != endLeft)
            itrLeft.seek(indexLeft);
        itrRight.seek(indexRight);
    }
}

template<class InSeg1, class InSeg2>
void restrictionLowHigh(const RLEBase3<InSeg1> & lhs,
    const RLEBase3<InSeg2> & rhs, RLEBase3<InSeg1> & result)
{
    assert(lhs.numDimensions() < rhs.numDimensions());
    // check that lhs and rhs share at least one non-symmetrical dimension.
//    std::cout << "LHS stride " << lhs.stride(0) << " " << lhs.stride(1)
//        << " " << lhs.stride(2) << "\n";
//    std::cout << "RHS stride " << rhs.stride(0) << " " << rhs.stride(1)
//        << " " << rhs.stride(2) << "\n";
    assert( (lhs.stride(0) && rhs.stride(0)) ||
        (lhs.stride(1) && rhs.stride(1)) ||
        (lhs.stride(2) && rhs.stride(2)) );
    result.dimensions(rhs.nonSymmetricDimensions());
    result.clear();
    
    typename RLEBase<InSeg1>::ConstIterator itrLeft(lhs.unmarkedBegin());
    typename RLEBase<InSeg1>::ConstIterator endLeft(lhs.end());
    typename RLEBase<InSeg2>::ConstIterator itrRight(rhs.unmarkedBegin());
    typename RLEBase<InSeg2>::ConstIterator endRight(rhs.end());
    
    int64 indexLeft, indexRight, curX, curY, curZ;
    
    indexRight = itrRight.position();
    rhs.cartesianCoordinates(indexRight, curX, curY, curZ);
    indexLeft = lhs.linearIndex(curX, curY, curZ);
    
    itrLeft.seek(indexLeft);
    itrRight.seek(indexRight);
    
    while (itrRight != endRight)
    {
        int64 length = itrRight.runEnd() - indexRight + 1;
        if (rhs.axis() == lhs.axis() && itrLeft < endLeft)
        {
            // The length of the rhs runline is:
            //  infinity    if rhs.axis() != lhs.axis() (extrude!)
            //  infinity    if itrRight >= endRight (unbounded last runline)
            //  something   otherwise
            length = std::min(length, itrLeft.runEnd() - indexLeft + 1);
        }
        
        if (itrLeft.isMarked() && itrRight.isMarked())
        {
            // These cases are NOT the same!  Truncation is essential to getting
            // this procedure to work right.
            if (rhs.axis() == lhs.axis()) // shared runline direction
            {
                // extrusion will result in many parallel copies of the runline
                InSeg1 truncated(*itrLeft, indexLeft, indexLeft + length - 1);
                result.RLEBase<InSeg1>::mark(indexRight, indexRight+length-1,
                    truncated.data());
            }
            else
            {
                // extrusion is perpendicular to existing runline direction
                result.RLEBase<InSeg1>::mark(indexRight,
                    indexRight + length - 1, itrLeft->markAt(indexLeft));
            }
            
//            result.RLEBase<InSeg1>::mark(indexRight,
//                indexRight + length - 1, itrLeft->markAt(indexLeft));
        }
        
        indexRight += length;
        rhs.cartesianCoordinates(indexRight, curX, curY, curZ);
        indexLeft = lhs.linearIndex(curX, curY, curZ);
        
        if (itrRight != endRight)
            itrRight.seek(indexRight);
        itrLeft.seek(indexLeft);
    }
}

template<class Segment, class UnaryPredicate>
void filter(const RLEBase3<Segment> & lhs, RLEBase3<Segment> & result,
    UnaryPredicate op)
{
    result.dimensions(lhs.nonSymmetricDimensions());
    result.clear();
    
    filter((const RLEBase<Segment> &) lhs, (RLEBase<Segment>&) result, op);
}

}; // namespace RLE


#endif
