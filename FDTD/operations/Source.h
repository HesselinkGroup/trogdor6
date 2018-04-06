/*
 *  Source.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/24/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _SOURCEJM_
#define _SOURCEJM_

#include "OutputEHDB.h"
#include "../SimulationDescription.h"
#include "../FieldEnum.h"
#include "../GridFields.h"
#include "../Operation.h"

typedef RunlineOutput RunlineSource;

class SourceJM : public Operation
{
public:
    SourceJM() {}
    SourceJM(const CurrentSourceDescription & description);
    
    void runlinesJ(int xyz, std::vector<RunlineSource> & runlines)
    { mRunlinesJ[xyz] = runlines; }
    void runlinesM(int xyz, std::vector<RunlineSource> & runlines)
    { mRunlinesM[xyz] = runlines; }
    void runlinesJE(int xyz, std::vector<RunlineSource> & runlines)
    { mRunlinesJE[xyz] = runlines; }
    void runlinesMH(int xyz, std::vector<RunlineSource> & runlines)
    { mRunlinesMH[xyz] = runlines; }
    
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields)
    {
        for (int xyz = 0; xyz < 3; xyz++)
        {
            mHeadJ[xyz] = currentGridFields.j(xyz);
            mHeadJE[xyz] = currentGridFields.je(xyz);
            mHeadM[xyz] = currentGridFields.m(xyz);
            mHeadMH[xyz] = currentGridFields.mh(xyz);
        }
    }
    
    unsigned long bytes() const
    {
        long totalBytes = 0;
        for (int xyz = 0; xyz < 3; xyz++)
        {
            totalBytes += mRunlinesJ[xyz].size()*sizeof(RunlineSource);
            totalBytes += mRunlinesM[xyz].size()*sizeof(RunlineSource);
            totalBytes += mRunlinesJE[xyz].size()*sizeof(RunlineSource);
            totalBytes += mRunlinesMH[xyz].size()*sizeof(RunlineSource);
        }
        totalBytes += mBuffer.size()*sizeof(Precision::Float);
        return totalBytes;
    }
    
    void apply(long timestep, Precision::Float dt) {}
    
    void calcJ(long timestep); // and J_E
    void calcM(long timestep); // and M_H
    
    long numCells() const
    {
        long cells = 0;
        for (int xyz = 0; xyz < 3; xyz++)
        {
            for (int nn = 0; nn < mRunlinesJ[xyz].size(); nn++)
                cells += mRunlinesJ[xyz][nn].length();
            for (int nn = 0; nn < mRunlinesM[xyz].size(); nn++)
                cells += mRunlinesM[xyz][nn].length();
            for (int nn = 0; nn < mRunlinesJE[xyz].size(); nn++)
                cells += mRunlinesJE[xyz][nn].length();
            for (int nn = 0; nn < mRunlinesMH[xyz].size(); nn++)
                cells += mRunlinesMH[xyz][nn].length();
        }
        return cells;
    }
    
    void allocate();
private:
    enum SourceType { kTimeFile, kSpaceTimeFile };
    SourceType mType;
    std::ifstream mFileStream;
    std::vector<RunlineSource> mRunlinesJ[3];
    std::vector<RunlineSource> mRunlinesM[3];
    std::vector<RunlineSource> mRunlinesJE[3];
    std::vector<RunlineSource> mRunlinesMH[3];
    Precision::Float* mHeadJ[3];
    Precision::Float* mHeadM[3];
    Precision::Float* mHeadJE[3];
    Precision::Float* mHeadMH[3];
    std::vector<Precision::Float> mBuffer;
    
    std::vector<Duration> mDurations;
    int mCurrentDuration;
    
    unsigned long totalLength(const std::vector<RunlineSource> & rl) const;
};


class SourceEH : public Operation
{
public:
    SourceEH() {}
    SourceEH(const SourceDescription & description);
    
    void runlinesE(int xyz, std::vector<RunlineSource> & runlines)
    { mRunlinesE[xyz] = runlines; }
    void runlinesH(int xyz, std::vector<RunlineSource> & runlines)
    { mRunlinesH[xyz] = runlines; }
    
    
    void setPointers(GridFields & currentGridFields,
        std::map<int, Pointer<GridFields> > & allGridFields)
    {
        for (int xyz = 0; xyz < 3; xyz++)
        {
            mHeadE[xyz] = currentGridFields.e(xyz);
            mHeadH[xyz] = currentGridFields.h(xyz);
        }
    }
    
    unsigned long bytes() const
    {
        long totalBytes = 0;
        for (int xyz = 0; xyz < 3; xyz++)
        {
            totalBytes += mRunlinesE[xyz].size()*sizeof(RunlineSource);
            totalBytes += mRunlinesH[xyz].size()*sizeof(RunlineSource);
        }
        return totalBytes;
    }
    
    void apply(long timestep, Precision::Float dt) {}
    
    void calcE(long timestep);
    void calcH(long timestep);
    
    void allocate() {}
private:
    enum SourceType { kTimeFile, kSpaceTimeFile };
    SourceType mType;
    std::ifstream mFileStream;
    std::vector<RunlineSource> mRunlinesE[3];
    std::vector<RunlineSource> mRunlinesH[3];
    Precision::Float* mHeadE[3];
    Precision::Float* mHeadH[3];
    
    std::vector<Duration> mDurations;
    int mCurrentDuration;
};















#endif
