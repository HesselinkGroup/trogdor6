/*
 *  DynamicRLE-inl.h
 *  refactorRLE
 *
 *  Created by Paul Hansen on 2/8/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _DYNAMICRLE_

namespace RLE
{

template<class T>
DynamicRLE<T>::
DynamicRLE() :
    RLEBase<DRLESegment<T> >()
{
}

template<class T>
DynamicRLE<T>::
DynamicRLE(const RLEBase<DRLESegment<T> > & copyMe) :
    RLEBase<DRLESegment<T> >(copyMe)
{
}

template<class T>
DynamicRLE<T> & DynamicRLE<T>::
operator=(const RLEBase<DRLESegment<T> > & rhs)
{
    ((RLEBase<DRLESegment<T> >&)*this) = rhs;
    return *this;
}

template<class T>
std::vector<T> DynamicRLE<T>::
values() const
{
    std::vector<T> vec(RLEBase<DRLESegment<T> >::numRuns());
    
    int64 curPos = 0;
    typename RLEBase<DRLESegment<T> >::AVLConstIterator itr;
    itr = RLEBase<DRLESegment<T> >::mSegmentTree.begin();
    while (itr.hasNext())
    {
        vec[curPos++] = itr.next().mark();
    }
    return vec;
}


#pragma mark *** Segment ***


template<class T>
DRLESegment<T>::
DRLESegment(int64 first, int64 last) :
    BaseSegment(first, last)
{
    // Warning: mMark is undefined
}

template<class T>
DRLESegment<T>::
DRLESegment(int64 first, int64 last, const T & val) :
    BaseSegment(first, last),
    mData(val)
{
}

template<class T>
DRLESegment<T>::
DRLESegment(const DRLESegment & copyMe, int64 first, int64 last) :
    BaseSegment(first, last),
    mData(copyMe.mData)
{
}

template<class T>
bool DRLESegment<T>::
operator==(const DRLESegment & rhs) const
{
    return (mData == rhs.mData) && (mFirst == rhs.mFirst) &&
        (mLast == rhs.mLast);
}

template<class T>
DRLESegment<T> DRLESegment<T>::
tail(int64 first) const
{
    return DRLESegment<T>(*this, first, mLast);
}

template<class T>
DRLESegment<T> DRLESegment<T>::
head(int64 last) const
{
    return DRLESegment<T>(*this, mFirst, last);
}

template<class T>
void DRLESegment<T>::
absorb(const DRLESegment & segment)
{
    assert(mData == segment.mData);
    assert(mFirst <= segment.last() + 1 && mLast >= segment.first()-1);
    
    if (segment.first() < mFirst)
        mFirst = segment.first();
    if (segment.last() > mLast)
        mLast = segment.last();
}

template<class T>
bool DRLESegment<T>::
canEat(const DRLESegment & segment) const
{
    assert(mFirst < segment.first());
    assert(mLast > segment.last());
    return mData == segment.mData;
}

template<class T>
bool DRLESegment<T>::
canMerge(const DRLESegment & segment) const
{
    return (mData == segment.mData);
}

#pragma mark *** Operators ***


// Binary operations with scalars

template<class T, class S>
DynamicRLE<T> operator+(const DynamicRLE<T> & lhs, const S & rhs)
{
    DynamicRLE<T> result;
    RLE::scalarTransform(lhs, rhs, result, std::plus<T>());
    return result;
}
template<class T, class S>
DynamicRLE<T> operator+(const S & lhs, const DynamicRLE<T> & rhs)
{
    return rhs + lhs;
}

template<class T, class S>
DynamicRLE<T> operator-(const DynamicRLE<T> & lhs, const S & rhs)
{
    DynamicRLE<T> result;
    RLE::scalarTransform(lhs, rhs, result, std::minus<T>());
    return result;
}

template<class T, class S>
struct AntiMinus
{
    T operator()(const T & lhs, const S & rhs) { return rhs - lhs; }
};

template<class T, class S>
DynamicRLE<T> operator-(const S & lhs, const DynamicRLE<T> & rhs)
{
    DynamicRLE<T> result;
    RLE::scalarTransform(rhs, lhs, result, AntiMinus<T,S>());
    return result;
}

template<class T, class S>
DynamicRLE<T> operator*(const DynamicRLE<T> & lhs, const S & rhs)
{
    DynamicRLE<T> result;
    RLE::scalarTransform(lhs, rhs, result, std::multiplies<T>());
    return result;
}
template<class T, class S>
DynamicRLE<T> operator*(const S & lhs, const DynamicRLE<T> & rhs)
{
    return rhs * lhs;
}

template<class T, class S>
DynamicRLE<T> operator/(const DynamicRLE<T> & lhs, const S & rhs)
{
    DynamicRLE<T> result;
    RLE::scalarTransform(lhs, rhs, result, std::divides<T>());
    return result;
}
template<class T, class S>
struct AntiDivides
{
    T operator()(const T & lhs, const S & rhs) { return rhs / lhs; }
};
template<class T, class S>
DynamicRLE<T> operator/(const S & lhs, const DynamicRLE<T> & rhs)
{
    DynamicRLE<T> result;
    RLE::scalarTransform(rhs, lhs, result, AntiDivides<T,S>());
    return result;
}


template<class T, class S>
DynamicRLE<T> operator%(const DynamicRLE<T> & lhs, const S & rhs)
{
    DynamicRLE<T> result;
    RLE::scalarTransform(lhs, rhs, result, std::modulus<T>());
    return result;
}
template<class T, class S>
struct AntiModulus
{
    T operator()(const T & lhs, const S & rhs) { return rhs % lhs; }
};
template<class T, class S>
DynamicRLE<T> operator%(const S & lhs, const DynamicRLE<T> & rhs)
{
    DynamicRLE<T> result;
    RLE::scalarTransform(rhs, lhs, result, AntiModulus<T, S>());
    return result;
}

// Unary negate

template<class T>
DynamicRLE<T> operator-(const DynamicRLE<T> & lhs)
{
    DynamicRLE<T> result;
    RLE::transform(lhs, result, std::negate<T>());
    return result;
}

// Binary operations with other arrays

template<class T>
DynamicRLE<T> operator+(const DynamicRLE<T> & lhs, const DynamicRLE<T> & rhs)
{
    DynamicRLE<T> result;
    RLE::transform(lhs, rhs, result, std::plus<T>(), T(), T());
    return result;
}

template<class T>
DynamicRLE<T> operator-(const DynamicRLE<T> & lhs, const DynamicRLE<T> & rhs)
{
    DynamicRLE<T> result;
    RLE::transform(lhs, rhs, result, std::minus<T>(), T(), T());
    return result;
}

template<class T>
DynamicRLE<T> operator*(const DynamicRLE<T> & lhs, const DynamicRLE<T> & rhs)
{
    DynamicRLE<T> result;
    RLE::transform(lhs, rhs, result, std::multiplies<T>(), T(), T());
    return result;
}

template<class T>
DynamicRLE<T> operator/(const DynamicRLE<T> & lhs, const DynamicRLE<T> & rhs)
{
    DynamicRLE<T> result;
    RLE::transform(lhs, rhs, result, std::divides<T>(), T(), T());
    return result;
}

template<class T>
DynamicRLE<T> operator%(const DynamicRLE<T> & lhs, const DynamicRLE<T> & rhs)
{
    DynamicRLE<T> result;
    RLE::transform(lhs, rhs, result, std::modulus<T>(), T(), T());
    return result;
}

template<class Out, class S, class T>
DynamicRLE<Out> multiply(const DynamicRLE<S> & lhs, const DynamicRLE<T> & rhs)
{
    DynamicRLE<Out> result;
    RLE::transform(lhs, rhs, result, ScalarMultiply<Out,S,T>(), S(), T());
    return result;
}

}; // namespace RLE

#endif
