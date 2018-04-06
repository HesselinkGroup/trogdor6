/*
 *  MergeRunlines-inl.h
 *  rle
 *
 *  Created by Paul Hansen on 3/18/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _RUNLINE_

namespace RLE
{


template<class Callback>
void merge(std::vector<int64> indexArray[],
    std::vector<uint64> strideArray[],
    std::vector<uint64> lengthArray[],
    uint64 numStreams,
    Callback callback)
{
    std::vector<uint64> ii(numStreams, 0); // iterators for data streams
    std::vector<uint64> runPosition(numStreams, 0); // 0 to lengths(i)-1
    
    std::vector<int64> currentIndex(numStreams);
    std::vector<uint64> currentStride(numStreams);
    bool done = 0;
    
    for (int64 mm = 0; mm < numStreams; mm++)
    {
        assert(indexArray[mm].size() == strideArray[mm].size());
        assert(indexArray[mm].size() == lengthArray[mm].size());
    }
    if (indexArray[0].size() == 0)
        return;
    
    while (!done)
    {
        // Step 1: calculate the end of the current run, and store the indices
        // and strides in each array to send to the callback.
        uint64 length = std::numeric_limits<uint64>::max();
        for (int64 mm = 0; mm < numStreams; mm++)
        if (ii[mm] < indexArray[mm].size())
        {
            uint64 lengthRemaining = lengthArray[mm][ii[mm]] - runPosition[mm];
            length = std::min(length, lengthRemaining);
            currentIndex[mm] = indexArray[mm][ii[mm]] +
                strideArray[mm][ii[mm]]*runPosition[mm];
            currentStride[mm] = strideArray[mm][ii[mm]];
        }
        assert(length < std::numeric_limits<uint64>::max());
        
        // Step 2: report the merged run, using the callback
        callback(&(currentIndex[0]), &(currentStride[0]), length, numStreams);
        
        // Step 3: update the current position
        done = 1;
        for (int64 mm = 0; mm < numStreams; mm++)
        {
            runPosition[mm] += length;
            if (ii[mm] < indexArray[mm].size())
            {
                if (runPosition[mm] >= lengthArray[mm][ii[mm]])
                {
                    runPosition[mm] = 0;
                    ii[mm]++;
                }
            }
            if (ii[mm] < indexArray[mm].size())
                done = 0;
        }
    }
}

template<class Callback, class IndArray>
void merge(IndArray const* arrays[], uint64 numStreams, Callback callback)
{
    if (arrays[0]->numRuns() == 0)
        return;
    
    std::vector<typename IndArray::AVLConstIterator> itrs(numStreams);
    for (int mm = 0; mm < numStreams; mm++)
        itrs[mm] = arrays[mm]->avlBegin();
        
    int64 curPos = std::numeric_limits<int64>::max();
    for (int64 mm = 0; mm < numStreams; mm++)
    if (itrs[mm].hasNext())
        curPos = std::min(curPos, itrs[mm].peekNext().first());
    
    std::vector<uint64> ii(numStreams, 0); // iterators for data streams
    std::vector<uint64> runPosition(numStreams, 0); // 0 to lengths(i)-1
    
    std::vector<int64> currentIndex(numStreams);
    std::vector<uint64> currentStride(numStreams);
    bool done = 0;
    
    while (!done)
    {
        // Step 1: calculate the end of the current run, and store the indices
        // and strides in each array to send to the callback.
        uint64 length = std::numeric_limits<uint64>::max();
        for (int64 mm = 0; mm < numStreams; mm++)
        if (itrs[mm].hasNext())
        {
            uint64 lengthRemaining = itrs[mm].peekNext().length() - runPosition[mm];
            length = std::min(length, lengthRemaining);
            
            currentIndex[mm] = itrs[mm].peekNext().markAt(
                itrs[mm].peekNext().first() + runPosition[mm]);
            currentStride[mm] = itrs[mm].peekNext().data().stride();
        }
        assert(length < std::numeric_limits<uint64>::max());
        
        // Step 2: report the merged run, using the callback
        callback(&(currentIndex[0]), &(currentStride[0]), length, numStreams);
        
        // Step 3: update the current position
        curPos += length;
        done = 1;
        for (int64 mm = 0; mm < numStreams; mm++)
        {
            runPosition[mm] += length;
            if (itrs[mm].hasNext())
            {
                if (runPosition[mm] >= itrs[mm].peekNext().length())
                {
                    runPosition[mm] = 0;
                    itrs[mm].next();
                }
            }
            if (itrs[mm].hasNext())
                done = 0;
        }
    }
}

template<class Callback>
void mergeOverlapping(std::vector<int64> indexArray[],
    std::vector<uint64> strideArray[],
    std::vector<uint64> lengthArray[],
    std::vector<int64> startArray[],
    uint64 numStreams,
    Callback callback)
{
    std::vector<uint64> ii(numStreams, 0); // iterators for data streams
    uint64 currentPosition = 0; // ranges from 0 to numCells - 1
    std::vector<uint64> runPosition(numStreams, 0); // 0 to lengths(i)-1
    std::vector<uint64> lengthRemaining(numStreams, 0);
    
    std::vector<int64> currentIndex(numStreams);
    std::vector<uint64> currentStride(numStreams);
    bool done = 0;
    
    for (int64 mm = 0; mm < numStreams; mm++)
    {
        assert(indexArray[mm].size() == strideArray[mm].size());
        assert(indexArray[mm].size() == lengthArray[mm].size());
        assert(indexArray[mm].size() == startArray[mm].size());
    }
//    if (indexArray[0].size() == 0)
//        return;
    
    int64 curPos = std::numeric_limits<int64>::max();
    for (int64 mm = 0; mm < numStreams; mm++)
    if (startArray[mm].size() > 0)
        curPos = std::min(curPos, startArray[mm][0]);
    
//    std::cerr << "First cur pos = " << curPos << "\n";
    
    while (!done)
    {
        // Step 1: calculate the end of the current run, and store the indices
        // and strides in each array to send to the callback.
        bool somethingIsMarked = false;
        int64 nextPos = std::numeric_limits<int64>::max();
        for (int mm = 0; mm < numStreams; mm++)
        {
            std::vector<int64> & indices = indexArray[mm]; // the data
            std::vector<int64> & starts = startArray[mm]; // run start
            std::vector<uint64> & lengths = lengthArray[mm]; // run length
            std::vector<uint64> & strides = strideArray[mm]; // run stride
            uint64 & iii = ii[mm];
            
            int64 myNextPos;
            if (iii < indices.size()) // if curRun < numRuns
            // if itr != rle.end()
            {
                if (curPos < starts[iii])
                // if curPos < itr.runStart()
                {
//                    std::cerr << curPos << " < " << startArray[mm][ii[mm]] << "\n";

                    myNextPos = starts[iii];
                    // myNextPos = itr.runStart()
                    currentIndex[mm] = -1;
                    currentStride[mm] = -1;
                }
                else
                {
//                    std::cerr << curPos << " >= " << startArray[mm][ii[mm]] << "\n";
                    myNextPos = starts[iii] + lengths[iii];
                    // myNextPos = itr.runEnd()+1
                    
//                    std::cerr << "\tmy next is " << startArray[mm][ii[mm]] << " + "
//                        << lengthArray[mm][ii[mm]] << "\n";
                    currentStride[mm] = strides[iii];
                    // currentStride[mm] = itr->data().stride()
                    currentIndex[mm] = indices[iii] + 
                        strides[iii]*(curPos - starts[iii]);
                    // currentIndex[mmm] = itr->markAt(curPos)
                    somethingIsMarked = true;
                }
                nextPos = std::min(nextPos, myNextPos);
//                std::cerr << "\tnew next = " << nextPos << "\n";
            }
            else
            {
//                std::cerr << "Stream " << mm << " is empty.\n";
                currentIndex[mm] = -1;
                currentStride[mm] = -1;
            }
        }
        
        // Step 2: check if we're done; if not, callback and step!
        if (nextPos == std::numeric_limits<int64>::max())
        {
//            std::cerr << "Done!\n";
            done = 1;
        }
        else
        {
//            std::cerr << "Not done... ";
            if (somethingIsMarked) // we never call the callback with NO indices
            {
                callback(&(currentIndex[0]), &(currentStride[0]),
                    nextPos - curPos, numStreams);
            }
            curPos = nextPos;
            
//            std::cerr << "cur pos is now " << curPos << "\n";
            
            for (int mm = 0; mm < numStreams; mm++)
            {
                std::vector<int64> & indices = indexArray[mm];
                std::vector<int64> & starts = startArray[mm];
                std::vector<uint64> & lengths = lengthArray[mm];
                std::vector<uint64> & strides = strideArray[mm];
                uint64 & iii = ii[mm];
                
                if (iii < indices.size())
                // if itr != rle.end()
                {
//                    std::cerr << "Stream " << mm << " not done (ii = "
//                        << ii[mm] << ")\n";
                    if (curPos >= starts[iii] + lengths[iii])
                    // if curPos > itr.runEnd()
                    {
//                        std::cerr << "\tcurPos " << curPos << " >= "
//                            << startArray[mm][ii[mm]] << " + "
//                            << lengthArray[mm][ii[mm]] << "\n";
                        iii++;
                        // itr.nextMarkedRun()
//                        std::cerr << "\tii[" << mm << "] = " << ii[mm] << "\n";
                    }
                    else
                    {
//                        std::cerr << "\tcurPos " << curPos << " < "
//                            << startArray[mm][ii[mm]] << " + "
//                            << lengthArray[mm][ii[mm]] << "\n";
                    }
                }
            }
        }
    }
}

// Originally: 1.2x slower than non-iterator version
// Caching IndexArray::end() for each array: 1.07x slower
//    xx Saving itr = itrs[mm] in loops: 1.07x slower (not keeping this one)
// Using AVLConstIterators: 0.93 the runtime (huzzah, it's faster)
// 
// This also should use less memory.

template<class Callback, class IndArray>
void mergeOverlapping(IndArray const* arrays[], uint64 numStreams, Callback callback)
{
    typedef typename IndArray::SegmentType seg;
    
    std::vector<typename IndArray::AVLConstIterator> itrs(numStreams);
    for (int mm = 0; mm < numStreams; mm++)
    {
        itrs[mm] = arrays[mm]->avlBegin();
    }
    
    std::vector<int64> currentIndex(numStreams);
    std::vector<uint64> currentStride(numStreams);
    bool done = 0;
    
    int64 curPos = std::numeric_limits<int64>::max();
    for (int64 mm = 0; mm < numStreams; mm++)
    if (itrs[mm].hasNext())
        curPos = std::min(curPos, itrs[mm].peekNext().first());
    
    while (!done)
    {
        // Step 1: calculate the end of the current run, and store the indices
        // and strides in each array to send to the callback.
        bool somethingIsMarked = false;
        int64 nextPos = std::numeric_limits<int64>::max();
        for (int mm = 0; mm < numStreams; mm++)
        {   
            int64 myNextPos;
            
            if (itrs[mm].hasNext())
            {
                if (curPos < itrs[mm].peekNext().first())
                {
                    myNextPos = itrs[mm].peekNext().first();
                    currentIndex[mm] = -1;
                    currentStride[mm] = -1;
                }
                else
                {
                    myNextPos = itrs[mm].peekNext().last()+1;
                    
                    currentStride[mm] = itrs[mm].peekNext().data().stride();
                    currentIndex[mm] = itrs[mm].peekNext().markAt(curPos);
                    somethingIsMarked = true;
                }
                nextPos = std::min(nextPos, myNextPos);
            }
            else
            {
                currentIndex[mm] = -1;
                currentStride[mm] = -1;
            }
        }
        
        // Step 2: check if we're done; if not, callback and step!
        if (nextPos == std::numeric_limits<int64>::max())
        {
            done = 1;
        }
        else
        {
            if (somethingIsMarked) // we never call the callback with NO indices
            {
                callback(&(currentIndex[0]), &(currentStride[0]),
                    nextPos - curPos, numStreams);
            }
            curPos = nextPos;
                        
            for (int mm = 0; mm < numStreams; mm++)
            {
                if (itrs[mm].hasNext())
                {
                    if (curPos > itrs[mm].peekNext().last())
                    {
                        itrs[mm].next();
                    }
                }
            }
        }
    }
}





};

#endif
