/*
 *  DynamicRLE.h
 *  refactorRLE
 *
 *  Created by Paul Hansen on 2/8/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _DYNAMICRLE_
#define _DYNAMICRLE_

#include "RLEBase.h"

namespace RLE
{

template<class T>
class DRLESegment : public BaseSegment
{
public:
    typedef T DataType; // argument for constructor
    typedef const T & DataReturnType;
    typedef T MarkType;
    typedef const T& MarkReturnType;
//    enum { hasDefault = true };
    
    DRLESegment(int64 first, int64 last);
    DRLESegment(int64 first, int64 last, const DataType & data);
    DRLESegment(const DRLESegment & copyMe, int64 first, int64 last); // truncation
    
    inline bool operator==(const DRLESegment & rhs) const;
    
    inline DataReturnType data() const { return mData; }
    inline MarkReturnType mark() const { return mData; }
    inline MarkReturnType markAt(int64 position) const { return mData; }
    static bool isDefaultMark(const DataType & mark) { return false; }
    
    inline DRLESegment tail(int64 first) const;
    inline DRLESegment head(int64 last) const;
    inline void absorb(const DRLESegment& segment);
    inline bool canEat(const DRLESegment & segment) const;
    inline bool canMerge(const DRLESegment & segment) const;
protected:
    DataType mData;
};

template<class T>
class DynamicRLE : public RLEBase<DRLESegment<T> >
{
public:
    DynamicRLE();
    DynamicRLE(const RLEBase<DRLESegment<T> > & copyMe);
    DynamicRLE<T> & operator=(const RLEBase<DRLESegment<T> > & rhs);
    
    std::vector<T> values() const;
};

#pragma mark *** Operators ***

// Binary operations with scalars

template<class T, class S>
DynamicRLE<T> operator+(const DynamicRLE<T> & lhs, const S & rhs);
template<class T, class S>
DynamicRLE<T> operator+(const S & lhs, const DynamicRLE<T> & rhs);
template<class T, class S>
DynamicRLE<T> operator-(const DynamicRLE<T> & lhs, const S & rhs);
template<class T, class S>
DynamicRLE<T> operator-(const S & lhs, const DynamicRLE<T> & rhs);
template<class T, class S>
DynamicRLE<T> operator*(const DynamicRLE<T> & lhs, const S & rhs);
template<class T, class S>
DynamicRLE<T> operator*(const S & lhs, const DynamicRLE<T> & rhs);
template<class T, class S>
DynamicRLE<T> operator/(const DynamicRLE<T> & lhs, const S & rhs);
template<class T, class S>
DynamicRLE<T> operator/(const S & lhs, const DynamicRLE<T> & rhs);
template<class T, class S>
DynamicRLE<T> operator%(const DynamicRLE<T> & lhs, const S & rhs);
template<class T, class S>
DynamicRLE<T> operator%(const S & lhs, const DynamicRLE<T> & rhs);

// Unary negate

template<class T>
DynamicRLE<T> operator-(const DynamicRLE<T> & lhs);

// Binary operations with other arrays

template<class T>
DynamicRLE<T> operator+(const DynamicRLE<T> & lhs, const DynamicRLE<T> & rhs);
template<class T>
DynamicRLE<T> operator-(const DynamicRLE<T> & lhs, const DynamicRLE<T> & rhs);
template<class T>
DynamicRLE<T> operator*(const DynamicRLE<T> & lhs, const DynamicRLE<T> & rhs);
template<class T>
DynamicRLE<T> operator/(const DynamicRLE<T> & lhs, const DynamicRLE<T> & rhs);
template<class T>
DynamicRLE<T> operator%(const DynamicRLE<T> & lhs, const DynamicRLE<T> & rhs);

// More flexible arithmetic
// If I want to define an arithmetic operation between types S and T, I can't
// assume that the return type is always S or always T.  I can templatize either
// "S operator*(S&, T&)" or "T operator*(S&, T&)", but not both.  I don't know
// of any way to actually overload the operator and get the result that I want.

template<class Out, class S, class T>
struct ScalarMultiply
{
    Out operator()(const S & lhs, const T & rhs)
    {
        return lhs*rhs;
    }
};

template<class Out, class S, class T>
DynamicRLE<Out> multiply(const DynamicRLE<S> & lhs, const DynamicRLE<T> & rhs);

}; // namespace RLE

#include "DynamicRLE-inl.h"

#endif
