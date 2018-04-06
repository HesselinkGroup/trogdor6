/*
 *  geometry2.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on Jan 6, 2010
 *  Copyright 2008 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _GEOMETRY2_
#define _GEOMETRY2_

#include "VectorMatrix2.h"


typedef Matrix2i Mat2i;
typedef Matrix2f Mat2f;
typedef Matrix2d Mat2d;
typedef Matrix2b Mat2b;

#pragma mark *** Vectors ***

/**
 * If relTol is positive, triangles that don't intersect might return true.
 * If relTol is negative, triangles with only coincident edges might return false.
 */
template<typename T>
bool segmentIntersects(const Vector2<T> & p0, const Vector2<T> & p1, const Vector2<T> & q0, const Vector2<T> & q1, const T & relTol);


#pragma mark *** Rect2 ***

template<typename T>
class Rect2
{
public:
	Rect2();
	Rect2(T i0, T j0, T i1, T j1);
	Rect2(const Vector2<T> & p0, const Vector2<T> & p1);
	Rect2(const Rect2<T> & copyMe);
	
    // size is calculated as p2-p1
	T size(unsigned int dim) const;
	Vector2<T> size() const;
    
    /**
     *  Return the signed volume of this rectangle.  The volume is the product
     *  of size() in each dimension.  If p1 and p2 represent the corners of
     *  the rectangle then this is the volume (for instance if this rectangle
     *  is floating-point-valued).  If p1 and p2 are inclusive indices of cells
     *  at corners, consider using count() instead.
     */
    T volume() const;
    
    // num is calculated as p2-p1+1 (so, "inclusive" size suitable for ints)
    T num(unsigned int dim) const;
    Vector2<T> num() const;
    
    /**
     *  Return the signed volume of this rectangle.  The volume is the product
     *  of num() in each dimension.  If p1 and p2 are inclusive indices of cells
     *  at corners, this will return the total number of cells in the rectangle.
     *  If p1 and p2 represent corners, consider using volume() instead.
     */
    T count() const;
	
	bool
	encloses(const Rect2<T> & inRect) const;
	
	bool
	encloses(const Vector2<T> & inPt) const;
	
	bool
	encloses(T x, T y) const;
	
	bool
	intersects(const Rect2<T> & inRect) const;
	
	int
	numNonSingularDims() const;
    
    bool operator==(const Rect2<T> & rhs) const
    {
        return (p1 == rhs.p1 && p2 == rhs.p2);
    }
    bool operator!=(const Rect2<T> & rhs) const
    {
        return (p1 != rhs.p1 || p2 != rhs.p2);
    }
	
	Vector2<T> p1;
	Vector2<T> p2; // must be indexwise < p1.
};

typedef Rect2<long> Rect2i;
typedef Rect2<double> Rect2d;

template<typename T>
Rect2<T>
operator * (const Rect2<T> & lhs, const T rhs);

template<typename T>
Rect2<T>
operator / (const Rect2<T> & lhs, const T rhs);

template<typename T>
Rect2<T>
operator * (const T scalar, const Rect2<T> & rhs);

template<typename T>
Rect2<T> &
operator *= (Rect2<T> & lhs, const T scalar);

template<typename T>
Rect2<T> &
operator /= (Rect2<T> & lhs, const T scalar);

template<typename T>
Rect2<T>
operator * (const Matrix2<T> & lhs, const Rect2<T> & rhs);

template<typename T>
Rect2<T>
operator + (const Rect2<T> & lhs, const Vector2<T> & rhs);

template<typename T>
Rect2<T>
operator - (const Rect2<T> & lhs, const Vector2<T> & rhs);

template<typename T>
Rect2<T>
operator + (const Vector2<T> & lhs, const Rect2<T> & rhs);

template<typename T>
Rect2<T>
operator - (const Vector2<T> & lhs, const Rect2<T> & rhs);

template<typename T>
Rect2<T>
inset( const Rect2<T> & inRect, T dx0, T dy0, T dx1, T dy1);

template<typename T>
Rect2<T>
inset( const Rect2<T> & inRect, T distance);

template<typename T>
Vector2<T>
clip( const Rect2<T> & clipRect, const Vector2<T> & v );

template<typename T>
Rect2<T>
clip( const Rect2<T> & clipRect, const Rect2<T> & rectToClip);

template<typename T>
Rect2<T>
cyclicPermute(const Rect2<T> & r, unsigned int nn);

template<typename T>
bool operator<(const Rect2<T> & lhs, const Rect2<T> & rhs);

template<typename T>
bool operator>(const Rect2<T> & lhs, const Rect2<T> & rhs);

#pragma mark *** Triangle ***

template<typename T>
class Triangle2
{
public:
    Triangle2()
    {
    }
    
    Triangle2(const Vector2<T> & p0, const Vector2<T> & p1,
        const Vector2<T> & p2)
    {
        mVertices[0] = p0;
        mVertices[1] = p1;
        mVertices[2] = p2;
    }
    
    Rect2<T> bounds() const
    {
        return Rect2<T>(
            vec_min(vec_min(mVertices[0], mVertices[1]), mVertices[2]),
            vec_max(vec_max(mVertices[0], mVertices[1]), mVertices[2]));
    }
    
    Vector2<T> & operator[](int vertex) { return mVertices[vertex]; }
    const Vector2<T> & operator[](int vertex) const
        { return mVertices[vertex]; }
    
    // Points on the boundary are considered to be inside.
    bool encloses(const Vector2<T> & point) const;
    
    Vector2<T> centroid() const { return (mVertices[0] + mVertices[1] + mVertices[2])/3.0; }
private:
    Vector2<T> mVertices[3];
};

template<class T>
bool intersects(const Triangle2<T> & t1, const Triangle2<T> & t2, const T & relTol);

typedef Triangle2<double> Triangle2d;

#pragma mark *** Input/Output ***

template<typename T>
std::ostream & operator<<(std::ostream & str, const Rect2<T> & rect);

template<typename T>
std::istream & operator>>(std::istream & str, Rect2<T> & rect);

template<typename T>
std::ostream & operator<<(std::ostream & str, const Triangle2<T> & triangle);

template<typename T>
std::istream & operator>>(std::istream & str, Triangle2<T> & triangle);

#include "geometry2-inl.h"

#endif
