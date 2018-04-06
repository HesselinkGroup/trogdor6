/*
 *  AdjointRunlineCallbacks.h
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 2/29/12.
 *  Copyright 2012 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef ADJOINTRUNLINECALLBACKS_H
#define ADJOINTRUNLINECALLBACKS_H


#include <vector>


struct CallbackAdjointDB
{
    CallbackAdjointDB(std::vector<RunlineAdjointDB> & runlines) :
        mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // db, eh, aux, coeff, coeff stride, length
        RunlineAdjointDB runline(indices[0], indices[1], indices[2],
            indices[3], strides[3], length);
        //cout << "made " << runline << "\n";
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineAdjointDB> & mRunlines;
};


struct CallbackAdjointDB_JM
{
    CallbackAdjointDB_JM(std::vector<RunlineAdjointDB_JM> & runlines) :
        mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // db, eh, aux, coeff, coeff stride, current, length
        RunlineAdjointDB_JM runline(indices[0], indices[1], indices[2],
            indices[3], strides[3], indices[4], length);
        //std::cerr << runline << "\n";
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineAdjointDB_JM> & mRunlines;
};

struct CallbackAdjointEH
{
    CallbackAdjointEH(std::vector<RunlineAdjointEH> & runlines) :
        mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // eh, gj1, gj2, gk1, gk2, aux, coeff, coeffStride, length
        RunlineAdjointEH runline(indices[0], indices[1], indices[2],
            indices[3], indices[4], indices[5], indices[6], strides[6],
            length);
//        cerr << "RL: " << runline << "\n";
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineAdjointEH> & mRunlines;
};

struct CallbackAdjointEH_JM
{
    CallbackAdjointEH_JM(std::vector<RunlineAdjointEH_JM> & runlines) :
        mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // eh, gj1, gj2, gk1, gk2, aux, coeff, coeffStride, current, length
        RunlineAdjointEH_JM runline(indices[0], indices[1], indices[2],
            indices[3], indices[4], indices[5], indices[6], strides[6],
            indices[7], length);
//        cerr << "RL: " << runline << "\n";
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineAdjointEH_JM> & mRunlines;
};





#endif
