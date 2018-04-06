/*
 *  VectorData.h
 *  support
 *
 *  Created by Paul Hansen on 12/10/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _VECTORDATA_
#define _VECTORDATA_

#include <iostream>
#include <vector>

namespace RLE
{

template<class T>
class VectorData
{
public:
    VectorData(int64 initialLength = 1);
    VectorData(int64 initialLength, const T & defaultMark);
    
    void mark(int64 first, int64 last, const T & mark);
    void erase(int64 first, int64 last);
    bool isMarked(int64 index) const;
    const T & markAt(int64 index) const;
    
    uint64 bytes() const;
    
    template<class U>
    friend std::ostream & operator<<(std::ostream & str,
        const VectorData<U> & s);
    
    std::vector<T> mMarks;
    T mDefaultMark;
};

template<class T>
std::ostream & operator<<(std::ostream & str, const VectorData<T> & r);

// Binary operations with scalars

template<class T, class S>
VectorData<T> operator+(const VectorData<T> & lhs, const S & rhs);
template<class T, class S>
VectorData<T> operator+(const S & lhs, const VectorData<T> & rhs);
template<class T, class S>
VectorData<T> operator-(const VectorData<T> & lhs, const S & rhs);
template<class T, class S>
VectorData<T> operator-(const S & lhs, const VectorData<T> & rhs);
template<class T, class S>
VectorData<T> operator*(const VectorData<T> & lhs, const S & rhs);
template<class T, class S>
VectorData<T> operator*(const S & lhs, const VectorData<T> & rhs);
template<class T, class S>
VectorData<T> operator/(const VectorData<T> & lhs, const S & rhs);
template<class T, class S>
VectorData<T> operator/(const S & lhs, const VectorData<T> & rhs);

// Unary negate

template<class T>
VectorData<T> operator-(const VectorData<T> & lhs);

// Binary operations with other arrays

template<class T>
VectorData<T> operator+(const VectorData<T> & lhs, const VectorData<T> & rhs);
template<class T>
VectorData<T> operator-(const VectorData<T> & lhs, const VectorData<T> & rhs);
template<class T>
VectorData<T> operator*(const VectorData<T> & lhs, const VectorData<T> & rhs);
template<class T>
VectorData<T> operator/(const VectorData<T> & lhs, const VectorData<T> & rhs);

}; // namespace RLE

#include "VectorData-inl.h"

#endif
