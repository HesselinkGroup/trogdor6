/*
 *  VectorData3.h
 *  refactorRLE
 *
 *  Created by Paul Hansen on 2/10/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _VECTORDATA3_
#define _VECTORDATA3_

#include "VectorData.h"

namespace RLE
{

template<class T>
class VectorData3 : public VectorData<T>
{
public:
    VectorData3() :
        VectorData<T>(1),
        m_nx(1), m_ny(1), m_nz(1), m_nxny(1)
    {
    }
    
    VectorData3(int64 nx, int64 ny, int64 nz, const T & defaultMark) :
        VectorData<T>(nx*ny*nz, defaultMark),
        m_nx(nx), m_ny(ny), m_nz(nz), m_nxny(nx*ny)
    {
    }
    
    int64 nx() const { return m_nx; }
    int64 ny() const { return m_ny; }
    int64 nz() const { return m_nz; }
    
    void mark(int64 x0, int64 y0, int64 z0, int64 x1, int64 y1, int64 z1, const T & val)
    {
        for (int64 zz = z0; zz <= z1; zz++)
        for (int64 yy = y0; yy <= y1; yy++)
        {
            int64 index = linearIndex(x0, yy, zz);
            VectorData<T>::mark(index, index + x1 - x0, val);
        }
    }
    
    void erase(int64 x0, int64 y0, int64 z0, int64 x1, int64 y1, int64 z1)
    {
        for (int64 zz = z0; zz <= z1; zz++)
        for (int64 yy = y0; yy <= y1; yy++)
        {
            int64 index = linearIndex(x0, yy, zz);
            VectorData<T>::erase(index, index + x1 - x0);
        }
    }
    
    bool isMarked(int64 x0, int64 y0, int64 z0) const
    {
        return VectorData<T>::isMarked(linearIndex(x0, y0, z0));
    }
    
    const T & markAt(int64 ii, int64 jj, int64 kk) const
    {
        return VectorData<T>::markAt(linearIndex(ii, jj, kk));
    }
    
    int64 linearIndex(int64 ii, int64 jj, int64 kk) const
    {
        return ii + jj*m_nx + kk*m_nxny;
    }
    
//    void cartesianCoordinates(int64 index, int64 & ii, int64 & jj, int64 & kk) const
//    {
//    }
protected:
    int64 m_nx, m_ny, m_nz;
    int64 m_nxny;
};

// Binary operations with scalars

template<class T, class S>
VectorData3<T> operator+(const VectorData3<T> & lhs, const S & rhs);
template<class T, class S>
VectorData3<T> operator+(const S & lhs, const VectorData3<T> & rhs);
template<class T, class S>
VectorData3<T> operator-(const VectorData3<T> & lhs, const S & rhs);
template<class T, class S>
VectorData3<T> operator-(const S & lhs, const VectorData3<T> & rhs);
template<class T, class S>
VectorData3<T> operator*(const VectorData3<T> & lhs, const S & rhs);
template<class T, class S>
VectorData3<T> operator*(const S & lhs, const VectorData3<T> & rhs);
template<class T, class S>
VectorData3<T> operator/(const VectorData3<T> & lhs, const S & rhs);
template<class T, class S>
VectorData3<T> operator/(const S & lhs, const VectorData3<T> & rhs);

// Unary negate

template<class T>
VectorData3<T> operator-(const VectorData3<T> & lhs);

// Binary operations with other arrays

template<class T>
VectorData3<T> operator+(const VectorData3<T> & lhs, const VectorData3<T> & rhs);
template<class T>
VectorData3<T> operator-(const VectorData3<T> & lhs, const VectorData3<T> & rhs);
template<class T>
VectorData3<T> operator*(const VectorData3<T> & lhs, const VectorData3<T> & rhs);
template<class T>
VectorData3<T> operator/(const VectorData3<T> & lhs, const VectorData3<T> & rhs);

}; // namespace RLE


#include "VectorData3-inl.h"

#endif
