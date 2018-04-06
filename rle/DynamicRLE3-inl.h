/*
 *  DynamicRLE3-inl.h
 *  rle_index
 *
 *  Created by Paul Hansen on 2/2/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _DYNAMICRLE3_

#include <cmath>

namespace RLE
{

template<class T>
DynamicRLE3<T>::
DynamicRLE3() :
    RLEBase3<DRLESegment<T> >(true, true, true)
{
}

template<class T>
DynamicRLE3<T>::
DynamicRLE3(bool hasX, bool hasY, bool hasZ) :
    RLEBase3<DRLESegment<T> >(hasX, hasY, hasZ)
{
}

template<class T>
template<class Scalar>
DynamicRLE3<T>::
DynamicRLE3(const Scalar nonSymmetricDimensions[]) :
    RLEBase3<DRLESegment<T> >(nonSymmetricDimensions)
{
}

//template<class T>
//DynamicRLE3<T>::
//DynamicRLE3(const RLEBase3<DRLESegment<T> > & copyMe) :
//    RLEBase3<DRLESegment<T> >(copyMe)
//{
//}

//template<class T>
//DynamicRLE3<T>::
//DynamicRLE3(const DynamicRLE3<T> & copyMe) :
//    RLEBase3<DRLESegment<T> >(copyMe)
//{
//}

template<class T>
DynamicRLE3<T> & DynamicRLE3<T>::
operator=(const DynamicRLE3<T> & rhs)
{
//    *this = rhs;
    ((RLEBase3<DRLESegment<T> >&)*this) = rhs;
    return *this;
}


template<class T>
std::vector<T> DynamicRLE3<T>::
values() const
{
    std::vector<T> vec(RLEBase3<DRLESegment<T> >::numRuns());
    
    int64 curPos = 0;
    typename RLEBase3<DRLESegment<T> >::AVLConstIterator itr;
    itr = RLEBase3<DRLESegment<T> >::mSegmentTree.begin();
    while (itr.hasNext())
    {
        vec[curPos++] = itr.next().mark();
    }
    
    return vec;
}


#pragma mark *** Operators ***


// Binary operations with scalars

template<class T, class S>
struct Multiplies
{
    T operator()(const T & lhs, const S & rhs)
    {
        return lhs*rhs;
    }
};

template<class T, class S>
struct Divides
{
    T operator()(const T & lhs, const S & rhs)
    {
        return lhs / rhs;
    }
};

template<class T, class S>
struct AntiDivides3
{
    T operator()(const T & lhs, const S & rhs) { return rhs / lhs; }
};

template<class T, class S>
struct AntiMultiplies3
{
    S operator()(const T & lhs, const S & rhs) { return rhs*lhs; }
};

template<class T, class S>
DynamicRLE3<T> operator+(const DynamicRLE3<T> & lhs, const S & rhs)
{
    DynamicRLE3<T> result;
    RLE::scalarTransform(lhs, rhs, result, std::plus<T>());
    return result;
}
template<class T, class S>
DynamicRLE3<T> operator+(const S & lhs, const DynamicRLE3<T> & rhs)
{
    return rhs + lhs;
}

template<class T, class S>
DynamicRLE3<T> operator-(const DynamicRLE3<T> & lhs, const S & rhs)
{
    DynamicRLE3<T> result;
    RLE::scalarTransform(lhs, rhs, result, std::minus<T>());
    return result;
}

template<class T, class S>
struct AntiMinus3
{
    T operator()(const T & lhs, const S & rhs) { return rhs - lhs; }
};

template<class T, class S>
DynamicRLE3<T> operator-(const S & lhs, const DynamicRLE3<T> & rhs)
{
    DynamicRLE3<T> result;
    RLE::scalarTransform(rhs, lhs, result, AntiMinus3<T,S>());
    return result;
}

template<class T, class S>
DynamicRLE3<T> operator*(const DynamicRLE3<T> & lhs, const S & rhs)
{
    DynamicRLE3<T> result;
    RLE::scalarTransform(lhs, rhs, result, Multiplies<T, S>());
    return result;
}
template<class T, class S>
DynamicRLE3<T> operator*(const S & lhs, const DynamicRLE3<T> & rhs)
{
    DynamicRLE3<T> result;
    RLE::scalarTransform(rhs, lhs, result, AntiMultiplies3<S, T>());
    return result;
}

template<class T, class S>
DynamicRLE3<T> operator/(const DynamicRLE3<T> & lhs, const S & rhs)
{
    DynamicRLE3<T> result;
    RLE::scalarTransform(lhs, rhs, result, Divides<T, S>());
    return result;
}

template<class T, class S>
DynamicRLE3<T> operator/(const S & lhs, const DynamicRLE3<T> & rhs)
{
    DynamicRLE3<T> result;
    RLE::scalarTransform(rhs, lhs, result, AntiDivides3<T,S>());
    return result;
}


template<class T, class S>
DynamicRLE3<T> operator%(const DynamicRLE3<T> & lhs, const S & rhs)
{
    DynamicRLE3<T> result;
    RLE::scalarTransform(lhs, rhs, result, std::modulus<T>());
    return result;
}
template<class T, class S>
struct AntiModulus3
{
    T operator()(const T & lhs, const S & rhs) { return rhs % lhs; }
};
template<class T, class S>
DynamicRLE3<T> operator%(const S & lhs, const DynamicRLE3<T> & rhs)
{
    DynamicRLE3<T> result;
    RLE::scalarTransform(rhs, lhs, result, AntiModulus<T, S>());
    return result;
}



// Unary negate

template<class T>
DynamicRLE3<T> operator-(const DynamicRLE3<T> & lhs)
{
    DynamicRLE3<T> result;
    RLE::transform(lhs, result, std::negate<T>());
    return result;
}

// Binary operations with other arrays

template<class T>
DynamicRLE3<T> operator+(const DynamicRLE3<T> & lhs, const DynamicRLE3<T> & rhs)
{
    DynamicRLE3<T> result;
    
    if (lhs.numRuns() == 0)
        result = rhs;
    else if (rhs.numRuns() == 0)
        result = lhs;
    else
        RLE::transform(lhs, rhs, result, std::plus<T>(), T(), T());
    
    return result;
}

template<class T>
DynamicRLE3<T> operator-(const DynamicRLE3<T> & lhs, const DynamicRLE3<T> & rhs)
{
    DynamicRLE3<T> result;
    
    if (lhs.numRuns() == 0)
        result = -rhs;
    else if (rhs.numRuns() == 0)
        result = lhs;
    else
        RLE::transform(lhs, rhs, result, std::minus<T>(), T(), T());
    
    return result;
}

template<class T>
DynamicRLE3<T> operator*(const DynamicRLE3<T> & lhs, const DynamicRLE3<T> & rhs)
{
    DynamicRLE3<T> result;
    
    RLE::transform(lhs, rhs, result, std::multiplies<T>(), T(), T());
    return result;
}

template<class T>
DynamicRLE3<T> operator/(const DynamicRLE3<T> & lhs, const DynamicRLE3<T> & rhs)
{
    DynamicRLE3<T> result;
    RLE::transform(lhs, rhs, result, std::divides<T>(), T(), T());
    return result;
}

template<class T>
DynamicRLE3<T> operator%(const DynamicRLE3<T> & lhs, const DynamicRLE3<T> & rhs)
{
    DynamicRLE3<T> result;
    RLE::transform(lhs, rhs, result, std::modulus<T>(), T(), T());
    return result;
}



}; // namespace RLE



#endif
