/*
 *  OutputEHDB.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/23/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _OUTPUTEHDB_
#define _OUTPUTEHDB_

#include "../SimulationDescription.h"
#include "../FieldEnum.h"
#include "../GridFields.h"

#include <iostream>
#include <fstream>

class GridFields;

struct RunlineOutput
{
    RunlineOutput() {}
    RunlineOutput(long head, long length) :
        mHead(head),
        mLength(length)
    {
        assert(length > 0);
    }
    
    long head() const { return mHead; }
    long length() const { return mLength; }
    
    long mHead;
    long mLength;
};
inline std::ostream & operator<<(std::ostream & str, const RunlineOutput & rl)
{
    str << "[head " << rl.head() << " length " << rl.length() << "]";
    return str;
}

struct RunlineTensorialOutput
{
    RunlineTensorialOutput() {}
    RunlineTensorialOutput(long head_ii, long head_i, long length) :
        mHead_ii(head_ii), mHead_i(head_i), mLength(length)
    {
        assert(length > 0);
    }
    
    long head_i() const { return mHead_i; }
    long head_ii() const { return mHead_ii; }
    long length() const { return mLength; }
    
    long mHead_ii, mHead_i;
    long mLength;
};
inline std::ostream & operator<<(std::ostream & str, const RunlineTensorialOutput & rl)
{
    str << "[head_ii " << rl.head_ii() << " head_i " << rl.head_i()
        << " length " << rl.length() << "]";
    return str;
}

struct RunlineAverageTwo
{
    RunlineAverageTwo() {}
    RunlineAverageTwo(long head0, long head1, long length) :
        mHead0(head0), mHead1(head1), mLength(length)
    {
        assert(length > 0);
    }
    
    long head0() const { return mHead0; }
    long head1() const { return mHead1; }
    long length() const { return mLength; }
    
    long mHead0, mHead1, mLength;
};
inline std::ostream & operator<<(std::ostream & str, const RunlineAverageTwo & rl)
{
    str << "[head0 " << rl.head0() << " head1 " << rl.head1()
        << " length " << rl.length() << "]";
    return str;
}

class HandleOutputBase
{
public:
    HandleOutputBase(const OutputField & f) { mField = f; }
    virtual ~HandleOutputBase() {}
    virtual long calc(std::ostream & binaryStream) = 0;
    virtual void setPointers(GridFields & fields) = 0;
    virtual long bytes() const = 0;
    virtual void printRunlines(std::ostream & str) const = 0;
    
    const OutputField & field() const { return mField; }
    
private:
    OutputField mField;
};

template<class RunlineType>
class HandleOutput : public HandleOutputBase
{
public:
    HandleOutput(const OutputField & field,
        const std::vector<RunlineType> & runlines) :
        HandleOutputBase(field),
        mRunlines(runlines)
    {
    }
    
    const std::vector<RunlineType> & runlines() const { return mRunlines; }
    std::vector<RunlineType> & runlines() { return mRunlines; }
    virtual void printRunlines(std::ostream & str) const
    {
        str << field().fieldAsString() << ":\n";
        for (int rr = 0; rr < mRunlines.size(); rr++)
            str << mRunlines[rr] << "\n";
    }
    virtual long bytes() const
    {
        return mRunlines.size() * sizeof(RunlineType);
    }
protected:
    std::vector<RunlineType> mRunlines;
};

class HandleRawField : public HandleOutput<RunlineOutput>
{
public:
    HandleRawField(const OutputField & field,
        const std::vector<RunlineOutput> & runlines);
    virtual long calc(std::ostream & binaryStream);
    virtual void setPointers(GridFields & fields);
private:
    Precision::Float* mHead;
};

class HandleTensorialDiagonalField : public HandleOutput<RunlineTensorialOutput>
{
public:
    HandleTensorialDiagonalField(const OutputField & field,
        const std::vector<RunlineTensorialOutput> & runlines);
    virtual long calc(std::ostream & binaryStream);
    virtual void setPointers(GridFields & fields);
private:
    Precision::Float *mHead0, *mHead1;
};


class HandleTensorialOffDiagonalField : public HandleOutput<RunlineTensorialOutput>
{
public:
    HandleTensorialOffDiagonalField(const OutputField & field,
        const std::vector<RunlineTensorialOutput> & runlines);
    virtual long calc(std::ostream & binaryStream);
    virtual void setPointers(GridFields & fields);
private:
    Precision::Float* mHead;
};

class HandleAveragedField : public HandleOutput<RunlineAverageTwo>
{
public:
    HandleAveragedField(const OutputField & field,
        const std::vector<RunlineAverageTwo> & runlines);
    virtual long calc(std::ostream & binaryStream);
    virtual void setPointers(GridFields & fields);
private:
    Precision::Float* mHead;
};


class OutputEHDB
{
public:
    OutputEHDB() {}
    OutputEHDB(const OutputDescription & description,
        Precision::Vec3 dxyz, Precision::Float dt, Precision::Vec3 origin);
    
    ~OutputEHDB();
    
    void appendField(const OutputField & field,
        const std::vector<RunlineOutput> & runlines);
    void appendField(const OutputField & field,
        const std::vector<RunlineAverageTwo> & runlines);
    void appendTensorialDiagonalField(const OutputField & field,
        const std::vector<RunlineTensorialOutput> & runlines);
    void appendTensorialOffDiagonalField(const OutputField & field,
        const std::vector<RunlineTensorialOutput> & runlines);
    
    void setPointers(GridFields & fields);
    
    long bytes() const
    {
        long totalBytes = 0;
        for (int ff = 0; ff < mFields.size(); ff++)
        {
            totalBytes += mFields.at(ff)->bytes();
        }
        return totalBytes;
    }
    
    void calc(int timestep);
    
    void printRunlines(std::ostream & str) const
    {
        for (int ff = 0; ff < mFields.size(); ff++)
            mFields[ff]->printRunlines(str);
    }
    
    std::string fileName() const
    {
        return mFileName;
    }
        
private:
    void writeOutputHeader(std::ostream & file, std::string fileName,
        std::string dataFileName,
        Precision::Vec3 dxyz, Precision::Float dt, Precision::Vec3 origin) const;
    void writeOutputFields(std::ostream & file,
        const OutputDescription & description) const;
    void writeOutputRegions(std::ostream & file,
        const std::vector<Region> & regions) const;
    void writeOutputDurations(std::ostream & file,
        const std::vector<Duration> & durations) const;
    
    std::ofstream mDataFile;
    std::string mFileName;
    
    std::vector<HandleOutputBase*> mFields;
    
    std::vector<Duration> mDurations;
    int mCurrentDuration;
};






#endif
