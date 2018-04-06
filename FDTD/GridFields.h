/*
 *  GridFields.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/18/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _GRIDFIELDS_
#define _GRIDFIELDS_

#include <vector>
#include "Precision.h"
#include "FieldEnum.h"

class ElectromagneticFieldIndices;

class GridFields
{
public:
    GridFields(const ElectromagneticFieldIndices & fieldIndices);
    
    Precision::Float* e(int xyz) { return &(mE[xyz][0]); }
    Precision::Float* h(int xyz) { return &(mH[xyz][0]); }
    Precision::Float* ee(int i, int j) { return &(mEE[i][j][0]); }
    Precision::Float* hh(int i, int j) { return &(mHH[i][j][0]); }
    Precision::Float* d(int xyz) { return &(mD[xyz][0]); }
    Precision::Float* b(int xyz) { return &(mB[xyz][0]); }
    Precision::Float* j(int xyz) { return &(mJ[xyz][0]); }
    Precision::Float* m(int xyz) { return &(mM[xyz][0]); }
    Precision::Float* je(int xyz) { return &(mJ_E[xyz][0]); }
    Precision::Float* mh(int xyz) { return &(mM_H[xyz][0]); }
    Precision::Float* field(Field f);
    
    // two deprecated-ish accessors
    Precision::Float* field(FieldName fName, int xyz);
    Precision::Float* field(FieldName fName, int i, int j);
    
    void clearJ()
    {
        for (int xyz = 0; xyz < 3; xyz++)
        {
            for (int nn = 0; nn < mJ[xyz].size(); nn++)
                mJ[xyz][nn] = 0;
            for (int nn = 0; nn < mJ_E[xyz].size(); nn++)
                mJ_E[xyz][nn] = 0;
        }
    }
    
    void clearM()
    {
        for (int xyz = 0; xyz < 3; xyz++)
        {
            for (int nn = 0; nn < mM[xyz].size(); nn++)
                mM[xyz][nn] = 0;
            for (int nn = 0; nn < mM_H[xyz].size(); nn++)
                mM_H[xyz][nn] = 0;
        }
    }
    
    long bytes() const;
    
    
private:
    std::vector<Precision::Float> mE[3];
    std::vector<Precision::Float> mH[3];
    std::vector<Precision::Float> mD[3];
    std::vector<Precision::Float> mB[3];
    std::vector<Precision::Float> mEE[3][3]; // anisotropic subcomponents
    std::vector<Precision::Float> mHH[3][3];
    std::vector<Precision::Float> mJ[3];
    std::vector<Precision::Float> mM[3];
    std::vector<Precision::Float> mJ_E[3];
    std::vector<Precision::Float> mM_H[3];
};





#endif
