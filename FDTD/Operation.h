/*
 *  Operation.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/31/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef T6_OPERATION
#define T6_OPERATION

#include <string>
#include <iostream>
#include <map>
#include <string>

#include "Precision.h"
#include "GridFields.h"
#include "Pointer.h"

class GridFields;

class Operation
{
public:
    Operation() {}
    virtual ~Operation() {}
    virtual void apply(long timestep, Precision::Float dt) = 0;
    virtual void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields) = 0;
    virtual void allocate() = 0;
    virtual unsigned long bytes() const = 0;
    
    virtual void printRunlines(std::ostream & str) const {  }
    virtual long numRunlines() const { return 0; }
    virtual long numCells() const { return 0; }
    
    void name(const std::string & nom) { mName = nom; }
    const std::string & name() const { return mName; }
private:
    std::string mName;
};

#endif


