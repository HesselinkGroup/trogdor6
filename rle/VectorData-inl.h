/*
 *  VectorData.cpp
 *  support
 *
 *  Created by Paul Hansen on 12/10/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _VECTORDATA_

namespace RLE
{

template<class T>
VectorData<T>::
VectorData(int64 initialLength) :
    mMarks(initialLength, 0),
    mDefaultMark(0)
{
}

template<class T>
VectorData<T>::
VectorData(int64 initialLength, const T & defaultMark) :
    mMarks(initialLength, defaultMark),
    mDefaultMark(defaultMark)
{
}

template<class T>
void VectorData<T>::
mark(int64 first, int64 last, const T & mark)
{
    assert(first >= 0);
    if (mark == mDefaultMark)
    {
        erase(first, last);
        return;
    }
    
    if (last >= mMarks.size())
        mMarks.resize(last+1, mDefaultMark);
    
    for (int64 nn = first; nn <= last; nn++)
        mMarks[nn] = mark;
}

template<class T>
void VectorData<T>::
erase(int64 first, int64 last)
{
    assert(first >= 0);
    
    if (last >= mMarks.size())
        last = mMarks.size()-1;
    
    for (int64 nn = first; nn <= last; nn++)
        mMarks[nn] = mDefaultMark;
}

template<class T>
bool VectorData<T>::
isMarked(int64 index) const
{
    assert(index >= 0);
    if (index < mMarks.size())
    {
        //std::cout << "Mark is " << mMarks[index] << endl;
        return (mMarks[index] != mDefaultMark);
    }
    return 0;
}

template<class T>
const T & VectorData<T>::
markAt(int64 index) const
{
    assert(index >= 0);
    if (index < mMarks.size())
    {
        return mMarks[index];
    }
    return mDefaultMark;
}

template<class T>
uint64 VectorData<T>::
bytes() const
{
    return sizeof(T)*mMarks.size();
}

template<class T>
std::ostream &
operator<<(std::ostream & str, const VectorData<T> & r)
{
    bool inMarkedSegment = 0;
    for (int64 nn = 0; nn < r.mMarks.size(); nn++)
    {
        if (inMarkedSegment)
        {
            if (r.mMarks[nn] == r.mDefaultMark) // End a segment
            {
                str << nn-1 << "] ";
                inMarkedSegment = 0;
            }
        }
        else
        {
            if (r.mMarks[nn] == 1) // start a segment
            {
                str << "[" << nn << ", ";
                inMarkedSegment = 1;
            }
        }
    }
    if (inMarkedSegment)
        str << r.mMarks.size()-1 << "] ";
    return str;
}


#pragma mark *** Operators ***


// Binary operations with scalars

template<class T, class S>
VectorData<T> operator+(const VectorData<T> & lhs, const S & rhs)
{
    VectorData<T> result(lhs.mMarks.size(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < lhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] + rhs;
    return result;
}
template<class T, class S>
VectorData<T> operator+(const S & lhs, const VectorData<T> & rhs)
{
    VectorData<T> result(rhs.mMarks.size(), rhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = rhs.mMarks[nn] + lhs;
    return result;
}

template<class T, class S>
VectorData<T> operator-(const VectorData<T> & lhs, const S & rhs)
{
    VectorData<T> result(lhs.mMarks.size(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < lhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] - rhs;
    return result;
}
template<class T, class S>
VectorData<T> operator-(const S & lhs, const VectorData<T> & rhs)
{
    VectorData<T> result(rhs.mMarks.size(), rhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs - rhs.mMarks[nn];
    return result;
}

template<class T, class S>
VectorData<T> operator*(const VectorData<T> & lhs, const S & rhs)
{
    VectorData<T> result(lhs.mMarks.size(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < lhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] * rhs;
    return result;
}

template<class T, class S>
VectorData<T> operator*(const S & lhs, const VectorData<T> & rhs)
{
    VectorData<T> result(rhs.mMarks.size(), rhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs * rhs.mMarks[nn];
    return result;
}

template<class T, class S>
VectorData<T> operator/(const VectorData<T> & lhs, const S & rhs)
{
    VectorData<T> result(lhs.mMarks.size(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < lhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] / rhs;
    return result;
}

template<class T, class S>
VectorData<T> operator/(const S & lhs, const VectorData<T> & rhs)
{
    VectorData<T> result(rhs.mMarks.size(), rhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs / rhs.mMarks[nn];
    return result;
}

// Unary negate

template<class T>
VectorData<T> operator-(const VectorData<T> & lhs)
{
    VectorData<T> result(lhs.mMarks.size(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < lhs.mMarks.size(); nn++)
        result.mMarks[nn] = -lhs.mMarks[nn];
    return result;
}

// Binary operations with other arrays

template<class T>
VectorData<T> operator+(const VectorData<T> & lhs, const VectorData<T> & rhs)
{
    assert(lhs.size() == rhs.size());
    VectorData<T> result(lhs.mMarks.size(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] + rhs.mMarks[nn];
    return result;
}

template<class T>
VectorData<T> operator-(const VectorData<T> & lhs, const VectorData<T> & rhs)
{
    assert(lhs.size() == rhs.size());
    VectorData<T> result(lhs.mMarks.size(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] - rhs.mMarks[nn];
    return result;
}

template<class T>
VectorData<T> operator*(const VectorData<T> & lhs, const VectorData<T> & rhs)
{
    assert(lhs.size() == rhs.size());
    VectorData<T> result(lhs.mMarks.size(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] * rhs.mMarks[nn];
    return result;
}

template<class T>
VectorData<T> operator/(const VectorData<T> & lhs, const VectorData<T> & rhs)
{
    assert(lhs.size() == rhs.size());
    VectorData<T> result(lhs.mMarks.size(), lhs.mDefaultMark);
    for (int64 nn = 0; nn < rhs.mMarks.size(); nn++)
        result.mMarks[nn] = lhs.mMarks[nn] / rhs.mMarks[nn];
    return result;
}



}; // namespace RLE

#endif
