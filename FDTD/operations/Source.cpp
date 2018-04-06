/*
 *  SourceJM.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/24/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "Source.h"

#include <cmath>

using namespace std;

/*
    TIME FILE:
    On each timestep, read one float for each field.  Thus it will take more
    than one float per timestep to handle circularly-polarized light, say.
    
    SPACE-TIME FILE:
    On each timestep, read one float for each field, for each cell.
*/

SourceJM::
SourceJM(const CurrentSourceDescription & description) :
    mDurations(description.durations()),
    mCurrentDuration(0)
{
    if (description.isSpaceVarying())
    {
        mType = kSpaceTimeFile;
        
        mFileStream.open(description.spaceTimeFile().c_str());
        if (!mFileStream.good())
            throw(std::logic_error("Cannot open source space-time file"));
    }
    else
    {
        mType = kTimeFile;
        
        mFileStream.open(description.timeFile().c_str());
        if (!mFileStream.good())
            throw(std::logic_error("Cannot open source time file"));
    }
}

void SourceJM::
calcJ(long timestep)
{
    if (mCurrentDuration >= mDurations.size())
        return;
    if (timestep < mDurations[mCurrentDuration].first())
        return;
    while (mCurrentDuration < mDurations.size() && 
        timestep > mDurations[mCurrentDuration].last())
        mCurrentDuration++;
    if (mCurrentDuration >= mDurations.size())
        return;
    
    Precision::Float val;
//    Precision::Float t0 = 40.0;
//    Precision::Float sigma = 10.0;
//    Precision::Float val = 1e6*exp(-(timestep-t0)*(timestep-t0)/sigma/sigma)*(timestep-t0);
    assert(mType == kTimeFile || mType == kSpaceTimeFile);
    assert(mBuffer.size() > 0);
    if (mType == kTimeFile)
    {
        for (int xyz = 0; xyz < 3; xyz++)
        if (mRunlinesJ[xyz].size() > 0)
        {
            mFileStream.read((char*)&val, sizeof(Precision::Float));
            for (int rr = 0; rr < mRunlinesJ[xyz].size(); rr++)
            {
                for (int ll = 0; ll < mRunlinesJ[xyz][rr].length(); ll++)
                    *(mHeadJ[xyz] + mRunlinesJ[xyz][rr].head() + ll) += val;
            }
        }
        
        for (int xyz = 0; xyz < 3; xyz++)
        if (mRunlinesJE[xyz].size() > 0)
        {
            mFileStream.read((char*)&val, sizeof(Precision::Float));
            for (int rr = 0; rr < mRunlinesJE[xyz].size(); rr++)
            {
                for (int ll = 0; ll < mRunlinesJE[xyz][rr].length(); ll++)
                    *(mHeadJE[xyz] + mRunlinesJE[xyz][rr].head() + ll) += val;
            }
        }
    }
    else if (mType == kSpaceTimeFile)
    {   
        for (int xyz = 0; xyz < 3; xyz++)
        if (mRunlinesJ[xyz].size() > 0)
        {
            long totalValues = totalLength(mRunlinesJ[xyz]);
            assert(totalValues <= mBuffer.size());
            mFileStream.read((char*)&(mBuffer[0]),
                totalValues*sizeof(Precision::Float));
            
            unsigned int bb = 0;
            for (int rr = 0; rr < mRunlinesJ[xyz].size(); rr++)
            {
                for (int ll = 0; ll < mRunlinesJ[xyz][rr].length(); ll++)
                {
                    *(mHeadJ[xyz] + mRunlinesJ[xyz][rr].head() + ll) += 
                        mBuffer[bb++];
                }
            }
        }
        
        for (int xyz = 0; xyz < 3; xyz++)
        if (mRunlinesJE[xyz].size() > 0)
        {
            long totalValues = totalLength(mRunlinesJE[xyz]);
            assert(totalValues <= mBuffer.size());
            mFileStream.read((char*)&(mBuffer[0]),
                totalValues*sizeof(Precision::Float));
            
            unsigned int bb = 0;
            for (int rr = 0; rr < mRunlinesJE[xyz].size(); rr++)
            {
                for (int ll = 0; ll < mRunlinesJE[xyz][rr].length(); ll++)
                {
                    *(mHeadJE[xyz] + mRunlinesJE[xyz][rr].head() + ll) += 
                        mBuffer[bb++];
                }
            }
        }
    }
}

void SourceJM::
calcM(long timestep)
{
    if (mCurrentDuration >= mDurations.size())
        return;
    if (timestep < mDurations[mCurrentDuration].first())
        return;
    while (mCurrentDuration < mDurations.size() && 
        timestep > mDurations[mCurrentDuration].last())
        mCurrentDuration++;
    if (mCurrentDuration >= mDurations.size())
        return;
    
    Precision::Float val;
    assert(mType == kTimeFile || mType == kSpaceTimeFile);
    assert(mBuffer.size() > 0);
    
    if (mType == kTimeFile)
    {
        for (int xyz = 0; xyz < 3; xyz++)
        if (mRunlinesM[xyz].size() > 0)
        {
            mFileStream.read((char*)&val, sizeof(Precision::Float));
            for (int rr = 0; rr < mRunlinesM[xyz].size(); rr++)
            {
                for (int ll = 0; ll < mRunlinesM[xyz][rr].length(); ll++)
                    *(mHeadM[xyz] + mRunlinesM[xyz][rr].head() + ll) += val;
            }
        }
        
        for (int xyz = 0; xyz < 3; xyz++)
        if (mRunlinesMH[xyz].size() > 0)
        {
            mFileStream.read((char*)&val, sizeof(Precision::Float));
            for (int rr = 0; rr < mRunlinesMH[xyz].size(); rr++)
            {
                for (int ll = 0; ll < mRunlinesMH[xyz][rr].length(); ll++)
                    *(mHeadMH[xyz] + mRunlinesMH[xyz][rr].head() + ll) += val;
            }
        }
    }
    else if (mType == kSpaceTimeFile)
    {
        if (sizeof(Precision::Float) != sizeof(val))
            throw(std::logic_error("WRONG FLOATING POINT SIZE"));
        
        for (int xyz = 0; xyz < 3; xyz++)
        if (mRunlinesM[xyz].size() > 0)
        {
            long totalValues = totalLength(mRunlinesM[xyz]);
            assert(totalValues <= mBuffer.size());
            mFileStream.read((char*)&(mBuffer[0]),
                totalValues*sizeof(Precision::Float));
            
            unsigned long bb = 0;
            for (int rr = 0; rr < mRunlinesM[xyz].size(); rr++)
            {
                for (int ll = 0; ll < mRunlinesM[xyz][rr].length(); ll++)
                {
                    *(mHeadM[xyz] + mRunlinesM[xyz][rr].head() + ll) += 
                        mBuffer[bb++];
                }
            }
        }
        
        for (int xyz = 0; xyz < 3; xyz++)
        if (mRunlinesMH[xyz].size() > 0)
        {
            long totalValues = totalLength(mRunlinesMH[xyz]);
            assert(totalValues <= mBuffer.size());
            mFileStream.read((char*)&(mBuffer[0]),
                totalValues*sizeof(Precision::Float));
            
            unsigned long bb = 0;
            for (int rr = 0; rr < mRunlinesMH[xyz].size(); rr++)
            {
                for (int ll = 0; ll < mRunlinesMH[xyz][rr].length(); ll++)
                {
                    *(mHeadMH[xyz] + mRunlinesMH[xyz][rr].head() + ll) += 
                        mBuffer[bb++];
                }
            }
        }
    }
}

void SourceJM::
allocate()
{
    long maxLength = 0;
    for (int xyz = 0; xyz < 3; xyz++)
    {
        long len = totalLength(mRunlinesJ[xyz]);
        if (len > maxLength)
            maxLength = len;
        len = totalLength(mRunlinesM[xyz]);
        if (len > maxLength)
            maxLength = len;
        len = totalLength(mRunlinesJE[xyz]);
        if (len > maxLength)
            maxLength = len;
        len = totalLength(mRunlinesMH[xyz]);
        if (len > maxLength)
            maxLength = len;
    }
    
    mBuffer = std::vector<Precision::Float>(maxLength);
}

unsigned long SourceJM::
totalLength(const std::vector<RunlineSource> & rl) const
{
    unsigned long total = 0;
    for (int rr = 0; rr < rl.size(); rr++)
        total += rl[rr].length();
    
    return total;
}

SourceEH::
SourceEH(const SourceDescription & description) :
    mDurations(description.durations()),
    mCurrentDuration(0)
{
    if (description.isSpaceVarying())
    {
        mType = kSpaceTimeFile;
        
        mFileStream.open(description.spaceTimeFile().c_str());
        if (!mFileStream.good())
            throw(std::logic_error("Cannot open source space-time file"));
    }
    else
    {
        mType = kTimeFile;
        
        mFileStream.open(description.timeFile().c_str());
        if (!mFileStream.good())
            throw(std::logic_error("Cannot open source time file"));
    }
}

void SourceEH::
calcE(long timestep)
{
    if (mCurrentDuration >= mDurations.size())
        return;
    if (timestep < mDurations[mCurrentDuration].first())
        return;
    while (mCurrentDuration < mDurations.size() && 
        timestep > mDurations[mCurrentDuration].last())
        mCurrentDuration++;
    if (mCurrentDuration >= mDurations.size())
        return;
    
    Precision::Float val;
    assert(mType == kTimeFile || mType == kSpaceTimeFile);
    
    if (mType == kTimeFile)
    {
        for (int xyz = 0; xyz < 3; xyz++)
        if (mRunlinesE[xyz].size() > 0)
        {
            mFileStream.read((char*)&val, sizeof(Precision::Float));
            for (int rr = 0; rr < mRunlinesE[xyz].size(); rr++)
            {
                for (int ll = 0; ll < mRunlinesE[xyz][rr].length(); ll++)
                    *(mHeadE[xyz] + mRunlinesE[xyz][rr].head() + ll) = val;
            }
        }
    }
    else if (mType == kSpaceTimeFile)
    {   
        for (int xyz = 0; xyz < 3; xyz++)
        if (mRunlinesE[xyz].size() > 0)
        {
            for (int rr = 0; rr < mRunlinesE[xyz].size(); rr++)
            {
                mFileStream.read((char*)(mHeadE[xyz]+mRunlinesE[xyz][rr].head()),
                    sizeof(Precision::Float)*mRunlinesE[xyz][rr].length());
            }
        }
    }
}

void SourceEH::
calcH(long timestep)
{
    if (mCurrentDuration >= mDurations.size())
        return;
    if (timestep < mDurations[mCurrentDuration].first())
        return;
    while (mCurrentDuration < mDurations.size() && 
        timestep > mDurations[mCurrentDuration].last())
        mCurrentDuration++;
    if (mCurrentDuration >= mDurations.size())
        return;
    
    Precision::Float val;
    assert(mType == kTimeFile || mType == kSpaceTimeFile);
    
    if (mType == kTimeFile)
    {
        for (int xyz = 0; xyz < 3; xyz++)
        if (mRunlinesH[xyz].size() > 0)
        {
            mFileStream.read((char*)&val, sizeof(Precision::Float));
            for (int rr = 0; rr < mRunlinesH[xyz].size(); rr++)
            {
                for (int ll = 0; ll < mRunlinesH[xyz][rr].length(); ll++)
                    *(mHeadH[xyz] + mRunlinesH[xyz][rr].head() + ll) = val;
            }
        }
    }
    else if (mType == kSpaceTimeFile)
    {
        if (sizeof(Precision::Float) != sizeof(val))
            throw(std::logic_error("WRONG FLOATING POINT SIZE"));
        
        for (int xyz = 0; xyz < 3; xyz++)
        if (mRunlinesH[xyz].size() > 0)
        {
            for (int rr = 0; rr < mRunlinesH[xyz].size(); rr++)
            {
                mFileStream.read((char*)(mHeadH[xyz]+mRunlinesH[xyz][rr].head()),
                    sizeof(Precision::Float)*mRunlinesH[xyz][rr].length());
            }
        }
    }
}


