/*
 *  DynamicRLE3.h
 *  rle_index
 *
 *  Created by Paul Hansen on 2/2/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _DYNAMICRLE3_
#define _DYNAMICRLE3_

#include "RLEBase3.h"
#include "DynamicRLE.h"

namespace RLE
{

template<class T>
class DynamicRLE3 : public RLEBase3<DRLESegment<T> >
{
public:
    DynamicRLE3();
    DynamicRLE3(bool hasX, bool hasY, bool hasZ);
    
    template<class Scalar>
    DynamicRLE3(const Scalar nonSymmetricDimensions[]);
    
//    DynamicRLE3(const DynamicRLE3<T> & copyMe);
    
//    DynamicRLE3(const RLEBase3<DRLESegment<T> > & copyMe);
//    DynamicRLE3<T> & operator=(const RLEBase3<DRLESegment<T> > & rhs);
    DynamicRLE3<T> & operator=(const DynamicRLE3<T> & rhs);
    
    std::vector<T> values() const;
};

#pragma mark *** Operators ***


// Binary operations with scalars

template<class T, class S>
DynamicRLE3<T> operator+(const DynamicRLE3<T> & lhs, const S & rhs);
template<class T, class S>
DynamicRLE3<T> operator+(const S & lhs, const DynamicRLE3<T> & rhs);
template<class T, class S>
DynamicRLE3<T> operator-(const DynamicRLE3<T> & lhs, const S & rhs);
template<class T, class S>
DynamicRLE3<T> operator-(const S & lhs, const DynamicRLE3<T> & rhs);
template<class T, class S>
DynamicRLE3<T> operator*(const DynamicRLE3<T> & lhs, const S & rhs);
template<class T, class S>
DynamicRLE3<T> operator*(const S & lhs, const DynamicRLE3<T> & rhs);
template<class T, class S>
DynamicRLE3<T> operator/(const DynamicRLE3<T> & lhs, const S & rhs);
template<class T, class S>
DynamicRLE3<T> operator/(const S & lhs, const DynamicRLE3<T> & rhs);
template<class T, class S>
DynamicRLE3<T> operator%(const DynamicRLE3<T> & lhs, const S & rhs);
template<class T, class S>
DynamicRLE3<T> operator%(const S & lhs, const DynamicRLE3<T> & rhs);

// Unary negate

template<class T>
DynamicRLE3<T> operator-(const DynamicRLE3<T> & lhs);

// Binary operations with other arrays

template<class T>
DynamicRLE3<T> operator+(const DynamicRLE3<T> & lhs, const DynamicRLE3<T> & rhs);
template<class T>
DynamicRLE3<T> operator-(const DynamicRLE3<T> & lhs, const DynamicRLE3<T> & rhs);
template<class T>
DynamicRLE3<T> operator*(const DynamicRLE3<T> & lhs, const DynamicRLE3<T> & rhs);
template<class T>
DynamicRLE3<T> operator/(const DynamicRLE3<T> & lhs, const DynamicRLE3<T> & rhs);


}; // namespace RLE

#include "DynamicRLE3-inl.h"

#endif
