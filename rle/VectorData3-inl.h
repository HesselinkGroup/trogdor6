/*
 *  VectorData3-inl.h
 *  refactorRLE
 *
 *  Created by Paul Hansen on 2/10/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _VECTORDATA3_

namespace RLE
{

#pragma mark *** Operators ***

// Binary operations with scalars

template<class T, class S>
VectorData3<T> operator+(const VectorData3<T> & lhs, const S & rhs)
{
    VectorData3<T> result(lhs.nx(), lhs.ny(), lhs.nz(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < lhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] + rhs;
    return result;
}
template<class T, class S>
VectorData3<T> operator+(const S & lhs, const VectorData3<T> & rhs)
{
    VectorData3<T> result(rhs.nx(), rhs.ny(), rhs.nz(), rhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = rhs.mMarks[nn] + lhs;
    return result;
}

template<class T, class S>
VectorData3<T> operator-(const VectorData3<T> & lhs, const S & rhs)
{
    VectorData3<T> result(lhs.nx(), lhs.ny(), lhs.nz(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < lhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] - rhs;
    return result;
}
template<class T, class S>
VectorData3<T> operator-(const S & lhs, const VectorData3<T> & rhs)
{
    VectorData3<T> result(rhs.nx(), rhs.ny(), rhs.nz(), rhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs - rhs.mMarks[nn];
    return result;
}

template<class T, class S>
VectorData3<T> operator*(const VectorData3<T> & lhs, const S & rhs)
{
    VectorData3<T> result(lhs.nx(), lhs.ny(), lhs.nz(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < lhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] * rhs;
    return result;
}

template<class T, class S>
VectorData3<T> operator*(const S & lhs, const VectorData3<T> & rhs)
{
    VectorData3<T> result(rhs.nx(), rhs.ny(), rhs.nz(), rhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs * rhs.mMarks[nn];
    return result;
}

template<class T, class S>
VectorData3<T> operator/(const VectorData3<T> & lhs, const S & rhs)
{
    VectorData3<T> result(lhs.nx(), lhs.ny(), lhs.nz(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < lhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] / rhs;
    return result;
}

template<class T, class S>
VectorData3<T> operator/(const S & lhs, const VectorData3<T> & rhs)
{
    VectorData3<T> result(rhs.nx(), rhs.ny(), rhs.nz(), rhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs / rhs.mMarks[nn];
    return result;
}

// Unary negate

template<class T>
VectorData3<T> operator-(const VectorData3<T> & lhs)
{
    VectorData3<T> result(lhs.nx(), lhs.ny(), lhs.nz(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < lhs.mMarks.size(); nn++)
        result.mMarks[nn] = -lhs.mMarks[nn];
    return result;
}

// Binary operations with other arrays

template<class T>
VectorData3<T> operator+(const VectorData3<T> & lhs, const VectorData3<T> & rhs)
{
    assert(lhs.mMarks.size() == rhs.mMarks.size());
    VectorData3<T> result(lhs.nx(), lhs.ny(), lhs.nz(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] + rhs.mMarks[nn];
    return result;
}

template<class T>
VectorData3<T> operator-(const VectorData3<T> & lhs, const VectorData3<T> & rhs)
{
    assert(lhs.mMarks.size() == rhs.mMarks.size());
    VectorData3<T> result(lhs.nx(), lhs.ny(), lhs.nz(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] - rhs.mMarks[nn];
    return result;
}

template<class T>
VectorData3<T> operator*(const VectorData3<T> & lhs, const VectorData3<T> & rhs)
{
    assert(lhs.mMarks.size() == rhs.mMarks.size());
    VectorData3<T> result(lhs.nx(), lhs.ny(), lhs.nz(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] * rhs.mMarks[nn];
    return result;
}

template<class T>
VectorData3<T> operator/(const VectorData3<T> & lhs, const VectorData3<T> & rhs)
{
    assert(lhs.mMarks.size() == rhs.mMarks.size());
    VectorData3<T> result(lhs.nx(), lhs.ny(), lhs.nz(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] / rhs.mMarks[nn];
    return result;
}


}; // namespace RLE

#endif
