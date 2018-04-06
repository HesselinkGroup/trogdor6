/*
 *  HuygensSurface.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 11/21/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef T6_HUYGENS_SURFACE
#define T6_HUYGENS_SURFACE

#include "../SimulationDescription.h"
#include "../FieldEnum.h"
#include "../GridFields.h"
#include "../Operation.h"
#include <cassert>
#include <iostream>

struct RunlineHuygens
{
    RunlineHuygens(unsigned long headJM, unsigned long headHE,
        unsigned long inDestStride, unsigned long inSrcStride,
        long length) :
        jm(headJM),
        he(headHE),
        srcStride(inSrcStride),
        destStride(inDestStride),
        mLength(length)
    {
        assert(length > 0);
    }
    
    long length() const { return mLength; }
    
    unsigned long jm;
    unsigned long he;
    unsigned long srcStride;
    unsigned long destStride;
    long mLength;
};
inline std::ostream & operator<<(std::ostream & str, const RunlineHuygens & rl)
{
    str << "[eh " << rl.he << " jm " << rl.jm << " length " << rl.length() << "];";
    return str;
}

class HuygensSurface : public Operation
{
public:
    HuygensSurface() {}
    HuygensSurface(GridDescPtr sourceGrid, double factor,
        Field writeField, Field readField,
        const std::vector<RunlineHuygens> & runlines);
    
    GridDescPtr sourceGrid() const { return mSourceGrid; }
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields);
    unsigned long bytes() const;
    void apply(long timestep, Precision::Float dt);
    void printRunlines(std::ostream & str) const;
    void allocate() {}
private:
    GridDescPtr mSourceGrid;
    Field mWriteField;
    Field mReadField;
    std::vector<RunlineHuygens> mRunlines;
    Precision::Float mFactor;
    Precision::Float* mHeadHE;
    Precision::Float* mHeadJM;
};


#endif
