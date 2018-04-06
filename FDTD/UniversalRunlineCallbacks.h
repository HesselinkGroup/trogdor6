/*
 *  UniversalRunlineCallbacks.h
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 2/29/12.
 *  Copyright 2012 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef UNIVERSALRUNLINECALLBACKS_H
#define UNIVERSALRUNLINECALLBACKS_H

#include "operations/JMPML.h"
#include "operations/OutputEHDB.h"
#include "operations/Source.h"
#include "operations/HuygensSurface.h"


// Generic callback that does nothing.
struct RLECallback
{
    RLECallback(int* numRuns = 0L) : mNumRuns(numRuns)
    {
        if (mNumRuns != 0L)
            *mNumRuns = 0;
    }
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
//        LOG << " length " << length << "\n";
        if (mNumRuns != 0L)
            (*mNumRuns)++;
    }
    
    int* mNumRuns;
};



struct CallbackHuygensSurface
{
    CallbackHuygensSurface(std::vector<RunlineHuygens> & runlines):
        mRunlines(runlines)
    {
    }
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        RunlineHuygens runline(indices[0], indices[1], strides[0], strides[1],
            length);
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineHuygens> & mRunlines;
};

struct CallbackJMPML
{
    CallbackJMPML(std::vector<RunlineJMPML> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        RunlineJMPML runline(indices[0], indices[1], indices[2],
            indices[3], indices[4], strides[4], length);
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineJMPML> & mRunlines;
};

struct CallbackSource
{
    CallbackSource(std::vector<RunlineOutput> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        RunlineOutput runline(indices[0], length);
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineOutput> & mRunlines;
};

struct CallbackPeriodicBC
{
    CallbackPeriodicBC(std::vector<RunlinePeriodicBC> & runlines):
        mRunlines(runlines)
    {}
    
    //
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // read, write, length
        RunlinePeriodicBC runline(indices[0], indices[1], length);
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlinePeriodicBC> & mRunlines;
};


#endif


