/*
 *  geometry.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 2/8/08.
 *  Copyright 2008 Stanford University.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _GEOMETRY2_

#pragma mark *** Vectors ***


template<typename T>
bool segmentIntersects(const Vector2<T> & p0, const Vector2<T> & p1, const Vector2<T> & q0, const Vector2<T> & q1, const T & relTol)
{
    Vector2<T> dp = p1 - p0;
    Vector2<T> dq = q1 - q0;
    
    // If segment p is all on the same side of line q, reject intersection.
    if (determinant(p0-q0, dq) * determinant(p1-q0, dq) >= 0.0)
    {
        return false;
    }
    
    // Create two rays:
    //  p0 + t*(p1-p0)
    //  q0 + s*(q1-q0)
    // Find their intersection times (t,s).
    // If t and s are both in [0,1] then the edges intersect.
    //
    // Implement some tolerance in here later if that's desirable.
    
    Matrix2<T> m = Matrix2<T>::withColumns(dp, -dq);
    Vector2<T> ts = inverse(m) * (q0 - p0);
    
    if (ts[0] > -relTol && ts[0] < (1.0+relTol) && ts[1] > -relTol && ts[1] < (1.0+relTol))
    {
        return true;
    }
    else
    {
        return false;
    }
}

#pragma mark *** Rect ***


template<typename T>
Rect2<T>::Rect2() :
    p1(0,0,0),
    p2(0,0,0)
{
}

template<typename T>
Rect2<T>::Rect2(T i0, T j0, T i1, T j1) :
    p1(i0, j0),
    p2(i1, j1)
{
}

template<typename T>
Rect2<T>::Rect2( const Vector2<T> & inP1,
    const Vector2<T> & inP2 ) :
	p1(inP1), p2(inP2)
{
}

template<typename T>
Rect2<T>::Rect2(const Rect2<T> & copyMe) :
    p1(copyMe.p1),
    p2(copyMe.p2)
{
}

template<typename T>
T Rect2<T>::size(unsigned int dim) const
{
	return p2[dim] - p1[dim];
}

template<typename T>
Vector2<T> Rect2<T>::size() const
{
	return p2 - p1 ;
}

template<typename T>
T Rect2<T>::volume() const
{
    return size(0)*size(1);
}

template<typename T>
T Rect2<T>::num(unsigned int dim) const
{
	return p2[dim] - p1[dim] + 1;
}

template<typename T>
Vector2<T> Rect2<T>::num() const
{
	return p2 - p1 + 1;
}

template<typename T>
T Rect2<T>::count() const
{
    return num(0)*num(1);
}


template<typename T>
bool Rect2<T>::encloses(const Rect2<T> & inRect) const
{
    return (encloses(inRect.p1) && encloses(inRect.p2));
}

template<typename T>
bool Rect2<T>::encloses(const Vector2<T> & inPt) const
{
    if ( inPt[0] >= p1[0] && inPt[0] <= p2[0] &&
         inPt[1] >= p1[1] && inPt[1] <= p2[1])
         return 1;
    return 0;
}

template<typename T>
bool Rect2<T>::encloses(T x, T y) const
{
    if ( x >= p1[0] && x <= p2[0] &&
         y >= p1[1] && y <= p2[1])
         return 1;
    return 0;
}

template<typename T>
bool Rect2<T>::intersects(const Rect2<T> & inRect) const
{
	return (vec_ge(p2, inRect.p1) && vec_le(p1, inRect.p2));
}

template<typename T>
int Rect2<T>::numNonSingularDims() const
{
    int ndims = 0;
    if (p1[0] != p2[0])
        ndims++;
    if (p1[1] != p2[1])
        ndims++;
    
    return ndims;
}


#pragma mark *** Other functions ***


template<typename T>
Rect2<T>
operator * (const Rect2<T> & lhs, const T rhs)
{
    return Rect2<T>(lhs.p1*rhs, lhs.p2*rhs);
}

template<typename T>
Rect2<T>
operator / (const Rect2<T> & lhs, const T rhs)
{
    return Rect2<T>(lhs.p1/rhs, lhs.p2/rhs);
}

template<typename T>
Rect2<T>
operator * (const T scalar, const Rect2<T> & rhs)
{
    return Rect2<T>(rhs.p1*scalar, rhs.p2*scalar);
}

template<typename T>
Rect2<T> &
operator *= (Rect2<T> & lhs, const T scalar)
{
    lhs.p1 *= scalar;
    lhs.p2 *= scalar;
    return lhs;
}

template<typename T>
Rect2<T> &
operator /= (Rect2<T> & lhs, const T scalar)
{
    lhs.p1 /= scalar;
    lhs.p2 /= scalar;
    return lhs;
}


template<typename T>
Rect2<T>
operator * (const Matrix2<T> & lhs, const Rect2<T> & rhs)
{
    return Rect2<T>(lhs*rhs.p1, lhs*rhs.p2);
}


template<typename T>
Rect2<T>
operator + (const Rect2<T> & lhs, const Vector2<T> & rhs)
{
	return Rect2<T>(lhs.p1 + rhs, lhs.p2 + rhs);
}

template<typename T>
Rect2<T>
operator - (const Rect2<T> & lhs, const Vector2<T> & rhs)
{
	return Rect2<T>(lhs.p1 - rhs, lhs.p2 - rhs);
}

template<typename T>
Rect2<T>
operator + (const Vector2<T> & lhs, const Rect2<T> & rhs)
{
	return Rect2<T>(rhs.p1 + lhs, rhs.p2 + lhs);
}

template<typename T>
Rect2<T>
operator - (const Vector2<T> & lhs, const Rect2<T> & rhs)
{
	return Rect2<T>(lhs - rhs.p1, lhs - rhs.p2);
}



template<typename T>
Rect2<T>
inset( const Rect2<T> & inRect, T dx0, T dx1, T dy0, T dy1)
{
    Rect2<T> out(inRect);
    out.p1[0] += dx0;
    out.p1[1] += dy0;
    out.p2[0] -= dx1;
    out.p2[1] -= dy1;
    return out;
}


template<typename T, typename S>
Rect2<T>
inset( const Rect2<T> & inRect, S distance )
{
    return inset(inRect, distance, distance, distance, distance);
}
                   
template<typename T>
Vector2<T>  clip(const Rect2<T> & clipRect, const Vector2<T> & v)
{
    Vector2<T> out(v);
    if (out[0] < clipRect.p1[0])
        out[0] = clipRect.p1[0];
    else if (out[0] > clipRect.p2[0])
        out[0] = clipRect.p2[0];
    
    if (out[1] < clipRect.p1[1])
        out[1] = clipRect.p1[1];
    else if (out[1] > clipRect.p2[1])
        out[1] = clipRect.p2[1];
    
    return out;
}


template<typename T>
Rect2<T> clip( const Rect2<T> & clipRect, const Rect2<T> & rectToClip)
{
	Rect2<T> out( clip(clipRect, rectToClip.p1), clip(clipRect, rectToClip.p2) );
	return out;
}

template<typename T>
Rect2<T> cyclicPermute(const Rect2<T> & r, unsigned int nn)
{
    return Rect2<T>(cyclicPermute(r.p1,nn), cyclicPermute(r.p2,nn));
}

template<typename T>
bool operator<(const Rect2<T> & lhs, const Rect2<T> & rhs)
{
    if (lhs.p1 < rhs.p1)
        return 1;
    else if (lhs.p1 > rhs.p1)
        return 0;
    else if (lhs.p2 < rhs.p2)
        return 1;
    return 0;
}


template<typename T>
bool operator>(const Rect2<T> & lhs, const Rect2<T> & rhs)
{
    if (lhs.p1 > rhs.p1)
        return 1;
    else if (lhs.p1 < rhs.p1)
        return 0;
    else if (lhs.p2 > rhs.p2)
        return 1;
    return 0;
}

#pragma mark *** Triangle ***

template<typename T>
bool Triangle2<T>::encloses(const Vector2<T> & point) const
{
    // Triangle is (a,b,c), point is p.
    // side_ab = sign(cross(b-a, p-a))
    // side_bc = sign(cross(c-b, p-b))
    // side_ca = sign(cross(a-c, p-c))
    // If side_ab == side_bc == side_ca then the point is inside the triangle
    
    T ab = cross(mVertices[1]-mVertices[0], point-mVertices[0]);
    T bc = cross(mVertices[2]-mVertices[1], point-mVertices[1]);
    T ca = cross(mVertices[0]-mVertices[2], point-mVertices[2]);
    
    return (ab*bc >= 0 && bc*ca >= 0);
}

template<class T>
bool intersects(const Triangle2<T> & t1, const Triangle2<T> & t2, const T & relTol)
{
    // 1. If any edges intersect, the triangles intersect.
    for (int ii = 0; ii < 3; ii++)
    {
        const Vector2<T> & p0 = t1[ii];
        const Vector2<T> & p1 = t1[(ii+1)%3];
        for (int jj = 0; jj < 3; jj++)
        {
            const Vector2<T> & q0 = t2[jj];
            const Vector2<T> & q1 = t2[(jj+1)%3];
            
            if (segmentIntersects(p0, p1, q0, q1, relTol))
            {
                return true;
            }
        }
    }
    
    // 2. If one triangle encloses a point of the other, they intersect
    
    if (t1.encloses(t2.centroid()))
    {
        return true;
    }
    else
    {
        return false;
    }
}

template<typename T>
std::ostream & operator<<(std::ostream & str, const Rect2<T> & rect)
{
	str << "[" << rect.p1 << ", " << rect.p2 << "]";
	return str;
}

template<typename T>
std::istream & operator>>(std::istream & str, Rect2<T> & rect)
{
	str >> rect.p1[0] >> rect.p1[1]
		>> rect.p2[0] >> rect.p2[1];
	return str;
}

template<typename T>
std::ostream & operator<<(std::ostream & str, const Triangle2<T> & triangle)
{
	str << "[" << triangle[0] << ", " << triangle[1] << ", "
        << triangle[2] << "]";
	return str;
}

template<typename T>
std::istream & operator>>(std::istream & str, Triangle2<T> & triangle)
{
	str >> triangle[0] >> triangle[1] >> triangle[2];
    
	return str;
}


#endif
