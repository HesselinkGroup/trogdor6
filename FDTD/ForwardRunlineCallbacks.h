/*
 *  ForwardRunlineCallbacks.h
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 2/28/12.
 *  Copyright 2012 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef FORWARD_RUNLINE_CALLBACKS_H
#define FORWARD_RUNLINE_CALLBACKS_H

#include <vector>

struct CallbackDB
{
    CallbackDB(std::vector<RunlineDB> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        RunlineDB runline(indices[0], indices[1], indices[2],
            indices[3], indices[4], length);
        //cout << "made " << runline << "\n";
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineDB> & mRunlines;
};

struct CallbackDB_JM
{
    CallbackDB_JM(std::vector<RunlineDB_JM> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        RunlineDB_JM runline(indices[0], indices[1], indices[2],
            indices[3], indices[4], indices[5], length);
        //std::cerr << runline << "\n";
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineDB_JM> & mRunlines;
};

struct CallbackNondispersiveEH
{
    CallbackNondispersiveEH(std::vector<RunlineNondispersiveEH> & runlines) :
        mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // eh, hej1, hej2, hek1, hek2, coeff, coeff stride, length
        RunlineNondispersiveEH runline(indices[0], indices[1], indices[2],
            indices[3], indices[4], indices[5], strides[5], length);
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineNondispersiveEH> & mRunlines;
};

struct CallbackNondispersiveEH_JM
{
    CallbackNondispersiveEH_JM(std::vector<RunlineNondispersiveEH_JM> & runlines) :
        mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // eh, hej1, hej2, hek1, hek2, coeff, coeff stride, current, length
        RunlineNondispersiveEH_JM runline(indices[0], indices[1], indices[2],
            indices[3], indices[4], indices[5], strides[5], indices[6], length);
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineNondispersiveEH_JM> & mRunlines;
};

struct CallbackEH
{
    CallbackEH(std::vector<RunlineEH> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // eh, db, aux, coeff, coeff stride, length
        RunlineEH runline(indices[0], indices[1], indices[2], indices[3],
            strides[3], length);
//        cerr << "RL: " << runline << "\n";
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineEH> & mRunlines;
};


struct CallbackSumEH_0
{
    CallbackSumEH_0(std::vector<RunlineSumEH_0> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        RunlineSumEH_0 runline(indices[0], indices[1], indices[2], indices[3],
            indices[4], indices[5], length);
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineSumEH_0> & mRunlines;
};

struct CallbackSumEH_i
{
    CallbackSumEH_i(std::vector<RunlineSumEH_i> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // input: ei, eii, eij, eik
        // ei, eii, eij, eik
        RunlineSumEH_i runline(indices[0], indices[1], indices[2], indices[3],
            length);
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineSumEH_i> & mRunlines;
};

struct CallbackSumEH_j
{
    CallbackSumEH_j(std::vector<RunlineSumEH_j> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // input: ei, eii, eij0-3, eik0-3
        RunlineSumEH_j runline(indices[0], indices[1],
            indices[2], indices[3], indices[4], indices[5], // eij
            indices[6], indices[7], indices[8], indices[9], // eik
            length);
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineSumEH_j> & mRunlines;
};

struct CallbackAnisotropicEH_0
{
    CallbackAnisotropicEH_0(std::vector<RunlineAnisotropicEH_0> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // eh, db, aux, coeff, coeff stride, length
        RunlineAnisotropicEH_0 runline(indices[0], indices[1], indices[2],
            indices[3], indices[4], strides[4], length);
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineAnisotropicEH_0> & mRunlines;
};

struct CallbackAnisotropicEH_i
{
    CallbackAnisotropicEH_i(std::vector<RunlineAnisotropicEH_i> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // input: eh, db0, db1, db2, db3, aux, eps
        // eh, db, aux, coeff, coeff stride, length
        RunlineAnisotropicEH_i runline(indices[0], indices[1], indices[2],
            indices[3], indices[4], indices[5], indices[6],
            strides[6], length);
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineAnisotropicEH_i> & mRunlines;
};

struct CallbackAnisotropicEH_j
{
    CallbackAnisotropicEH_j(std::vector<RunlineAnisotropicEH_j> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        // input: e, d, aux, eps
        // eh, db, aux, coeff, coeff stride, length
        RunlineAnisotropicEH_j runline(indices[0], indices[1], indices[2],
            indices[3], strides[3], length);
        mRunlines.push_back(runline);
    }
    
    std::vector<RunlineAnisotropicEH_j> & mRunlines;
};

#endif
