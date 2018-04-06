/*
 *  FieldEnum.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/23/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef T6_FIELDENUM
#define T6_FIELDENUM

enum FieldName
{
    kBeginElectricFields,
    kD,
    kE,
    kEE,
    kJ,
    kJE,
    kEndElectricFields,
    kBeginMagneticFields,
    kB,
    kH,
    kHH,
    kM,
    kMH,
    kEndMagneticFields
};

struct Field
{
public:
    Field() : mName(kD), m_ii(0), m_jj(0)
    {
        
    }
    
    Field(FieldName nn, int xyz) : mName(nn), m_ii(xyz), m_jj(-1)
    {
        
    }
    
    Field(FieldName nn, int row, int col) : mName(nn),
        m_ii(row), m_jj(col)
    {
        
    }
    
    static Field d(int xyz) { return Field(kD, xyz); }
    static Field e(int xyz) { return Field(kE, xyz); }
    static Field ee(int ii, int jj) { return Field(kEE, ii, jj); }
    static Field j(int xyz) { return Field(kJ, xyz); }
    static Field je(int xyz) { return Field(kJE, xyz); }
    static Field b(int xyz) { return Field(kB, xyz); }
    static Field h(int xyz) { return Field(kH, xyz); }
    static Field hh(int ii, int jj) { return Field(kHH, ii, jj); }
    static Field m(int xyz) { return Field(kM, xyz); }
    static Field mh(int xyz) { return Field(kMH, xyz); }
    
    int xyz() const { return m_ii; }
    int i() const { return m_ii; }
    int j() const { return m_jj; }
    FieldName name() const { return mName; }
    
    bool isElectric() const
    {
        return mName > kBeginElectricFields && mName < kEndElectricFields;
    }
    bool isMagnetic() const
    {
        return mName > kBeginMagneticFields && mName < kEndMagneticFields;
    }
    
private:
    FieldName mName;
    int m_ii, m_jj;
};

#endif
