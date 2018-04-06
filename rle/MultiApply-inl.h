/*
 *  MultiTransform.h
 *  rle
 *
 *  Created by Paul C Hansen on 10/12/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef MULTIAPPLY_H

namespace RLE
{

template<class Array, class Callback>
void multiApply(Array const* arrays[], int numArrays, Callback callback,
    int gapPolicy)
{
    std::vector<typename Array::AVLConstIterator> itrs(numArrays);
    for (int mm = 0; mm < numArrays; mm++)
    {
        itrs[mm] = arrays[mm]->avlBegin();
    }
    
    std::vector<typename Array::SegmentType::MarkType> currentData(numArrays);
    std::vector<int> markedArrays(numArrays, 0);
    bool done = 0;
    
    int64 curPos = std::numeric_limits<int64>::max();
    for (int64 mm = 0; mm < numArrays; mm++)
    if (itrs[mm].hasNext())
        curPos = std::min(curPos, itrs[mm].peekNext().first());
    
    while (!done)
    {
        // Step 1: calculate the end of the current run, and store the indices
        // and strides in each array to send to the callback.
        bool somethingIsMarked = false;
        int64 nextPos = std::numeric_limits<int64>::max();
        for (int mm = 0; mm < numArrays; mm++)
        {   
            int64 myNextPos;
            
            if (itrs[mm].hasNext())
            {
                if (curPos < itrs[mm].peekNext().first())
                {
                    myNextPos = itrs[mm].peekNext().first();
                    markedArrays[mm] = false;
                }
                else
                {
                    myNextPos = itrs[mm].peekNext().last()+1;
                    markedArrays[mm] = true;
                    currentData[mm] = itrs[mm].peekNext().markAt(curPos);
                    somethingIsMarked = true;
                }
                nextPos = std::min(nextPos, myNextPos);
            }
            else
            {
                markedArrays[mm] = false;
            }
        }
        
        // Step 2: check if we're done; if not, callback and step!
        if (nextPos == std::numeric_limits<int64>::max())
        {
            done = 1;
        }
        else
        {
            if (gapPolicy == RLE::IncludeGaps || somethingIsMarked)
            {
                callback(curPos, nextPos-1, &(currentData[0]),
                    &(markedArrays[0]), numArrays);
            }
            curPos = nextPos;
                        
            for (int mm = 0; mm < numArrays; mm++)
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

}; // namespace RLE

#endif
