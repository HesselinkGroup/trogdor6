/*
 *  OutputEHDB.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/24/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "OutputEHDB.h"

#include "TimeWrapper.h"
#include "../YeeUtilities.h"
#include "../GridFields.h"
#include "../OutputDirectory.h"

//#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>

using namespace std;
//using namespace boost::posix_time;
using namespace YeeUtilities;

HandleRawField::
HandleRawField(const OutputField & field, const vector<RunlineOutput> & runlines) :
    HandleOutput<RunlineOutput>(field, runlines),
    mHead(0L)
{
}

long HandleRawField::
calc(std::ostream & binaryStream)
{
    if (mRunlines.size() > 0)
    {
        assert(mHead != 0L);
    }
    
    long numCells = 0;
    for (int rr = 0; rr < mRunlines.size(); rr++)
    {
        numCells += runlines()[rr].length();
        if (runlines()[rr].head() == -1)
        {
            const Precision::Float zed = 0.0;
            for (int nn = 0; nn < runlines()[rr].length(); nn++)
                binaryStream.write((char*) &zed, sizeof(Precision::Float));
        }
        else
        {
            binaryStream.write((char*) (mHead + runlines()[rr].head()),
                (std::streamsize)(runlines()[rr].length()*sizeof(Precision::Float)));
        }
    }
    return numCells;
}

void HandleRawField::
setPointers(GridFields & fields)
{
    std::cout << "Field name " << field().fieldAsString() << "\n";
    assert(field().ii() == field().jj());
    mHead = fields.field(field().field(), field().ii());
}

//HandleSparseField::
//HandleSparseField(const OutputField & field, const vector<RunlineOutput> & runlines) :
//    HandleOutput<RunlineOutput>(field, runlines),
//    mHead(0L)
//{
//}
//
//long HandleSparseField::
//calc(std::ostream & binaryStream)
//{
//    assert(mHead != 0L);
//    long numCells = 0;
//    for (int rr = 0; rr < mRunlines.size(); rr++)
//    {
//        numCells += runlines()[rr].length();
//        if (runlines()[rr].head() != -1)
//        {
//            binaryStream.write((char*) (mHead + runlines()[rr].head()),
//                (std::streamsize)(runlines()[rr].length()*sizeof(Precision::Float)));
//        }
//        else
//        {
//            const Precision::Float zed = 0.0;
//            for (int nn = 0; nn < runlines()[rr].length(); nn++)
//                binaryStream.write((char*) &zed, sizeof(Precision::Float));
//        }
//    }
//}

//void HandleSparseField::
//setPointers(GridFields & fields)
//{
//    mHead = fields.field(field().field(), field().ii());
//}


HandleTensorialDiagonalField::
HandleTensorialDiagonalField(const OutputField & field,
    const vector<RunlineTensorialOutput> & runlines) :
    HandleOutput<RunlineTensorialOutput>(field, runlines),
    mHead0(0L),
    mHead1(0L)
{
}

long HandleTensorialDiagonalField::
calc(ostream & binaryStream)
{
    //assert(mHead0 != 0L); // might be null!
    assert(mHead1 != 0L); // should never be null.
    long numCells = 0;
    for (int rr = 0; rr < mRunlines.size(); rr++)
    {
        Precision::Float* headPointer;
        unsigned long offset;
        if (runlines()[rr].head_ii() != -1)
        {
            headPointer = mHead0;
            offset = runlines()[rr].head_ii();
        }
        else
        {
            headPointer = mHead1;
            offset = runlines()[rr].head_i();
        }
        numCells += runlines()[rr].length();
        binaryStream.write((char*) (headPointer + offset),
            (std::streamsize)(runlines()[rr].length()*sizeof(Precision::Float)));
    }
    return numCells;
}

void HandleTensorialDiagonalField::
setPointers(GridFields & fields)
{
    assert(field().ii() == field().jj());
    mHead0 = fields.field(field().field(), field().ii(), field().ii()); // may be null
    
    if (field().field() == kEE)
        mHead1 = fields.field(kE, field().ii());
    else if (field().field() == kHH)
        mHead1 = fields.field(kH, field().ii());
    else
        throw(std::logic_error("Bad field"));
    
    assert(mHead1 != 0L);
}


HandleTensorialOffDiagonalField::
HandleTensorialOffDiagonalField(const OutputField & field,
    const vector<RunlineTensorialOutput> & runlines) :
    HandleOutput<RunlineTensorialOutput>(field, runlines),
    mHead(0L)
{
}

long HandleTensorialOffDiagonalField::
calc(ostream & binaryStream)
{
    long numCells = 0;
    for (int rr = 0; rr < mRunlines.size(); rr++)
    {
        if (runlines()[rr].head_ii() != -1)
        {
            assert(mHead != 0L);
            binaryStream.write((char*) (mHead + runlines()[rr].head_ii()),
                (std::streamsize)(runlines()[rr].length()*sizeof(Precision::Float)));
            numCells += runlines()[rr].length();
        }
        else
        {
            const Precision::Float ZED(0.0);
            for (int ll = 0; ll < runlines()[rr].length(); ll++)
            {
                binaryStream.write((char*) &ZED,
                    (std::streamsize)sizeof(Precision::Float));
                numCells += 1;
            }
        }
    }
    return numCells;
}

void HandleTensorialOffDiagonalField::
setPointers(GridFields & fields)
{
    assert(field().ii() != field().jj());
    mHead = fields.field(field().field(), field().ii(), field().jj());
}

HandleAveragedField::
HandleAveragedField(const OutputField & field,
    const vector<RunlineAverageTwo> & runlines) :
    HandleOutput<RunlineAverageTwo>(field, runlines),
    mHead(0L)
{
}

long HandleAveragedField::
calc(ostream & binaryStream)
{
    assert(mHead != 0L);
    long numCells = 0;
    Precision::Float meanVal;
    for (int rr = 0; rr < mRunlines.size(); rr++)
    for (int ll = 0; ll < mRunlines[rr].length(); ll++)
    {
        meanVal = 0.5*(*(mHead + runlines()[rr].head0() + ll) +
            *(mHead + runlines()[rr].head1() + ll));
        
        binaryStream.write((char*) &meanVal,
            (std::streamsize)sizeof(Precision::Float));
        numCells += 1;
    }
    return numCells;
}

void HandleAveragedField::
setPointers(GridFields & fields)
{
    assert(field().ii() == field().jj());
    mHead = fields.field(field().field(), field().ii());
}

OutputEHDB::
OutputEHDB(const OutputDescription & description,
    Precision::Vec3 dxyz, Precision::Float dt, Precision::Vec3 origin) :
    mFileName(description.file()),
    mDurations(description.durations()),
    mCurrentDuration(0)
{
    string fileName = description.file() + ".txt";
    ofstream file(OutputDirectory::path(fileName).c_str());
    writeOutputHeader(file, fileName, description.file(), dxyz, dt, origin);
    writeOutputFields(file, description);
    writeOutputRegions(file, description.regions());
    writeOutputDurations(file, description.durations());
    file.close();
    
    mDataFile.open(OutputDirectory::path(description.file()).c_str());
    
    if (!mDataFile)
    {
        cerr << "Error opening file " << OutputDirectory::path(description.file()) << "\n";
        throw(std::runtime_error("Could not open output file!"));
    }
    
    mDataFile << setprecision(10);
}

OutputEHDB::
~OutputEHDB()
{
    mDataFile.close();
    
    for (int ff = 0; ff < mFields.size(); ff++)
        delete mFields.at(ff);
}


void OutputEHDB::
appendField(const OutputField & field, const vector<RunlineOutput> & runlines)
{
    mFields.push_back((HandleOutputBase*) new HandleRawField(field, runlines));
}

void OutputEHDB::
appendField(const OutputField & field,
    const vector<RunlineAverageTwo> & runlines)
{
    mFields.push_back((HandleOutputBase*) new HandleAveragedField(field, runlines));
}

void OutputEHDB::
appendTensorialDiagonalField(const OutputField & field,
    const vector<RunlineTensorialOutput> & runlines)
{
    mFields.push_back((HandleOutputBase*)
        new HandleTensorialDiagonalField(field, runlines));
}

void OutputEHDB::
appendTensorialOffDiagonalField(const OutputField & field,
    const vector<RunlineTensorialOutput> & runlines)
{
    mFields.push_back((HandleOutputBase*)
        new HandleTensorialOffDiagonalField(field, runlines));
}

void OutputEHDB::
setPointers(GridFields & fields)
{
    for (int ff = 0; ff < mFields.size(); ff++)
        mFields[ff]->setPointers(fields);
}


void OutputEHDB::
calc(int timestep)
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
    
    if ( (timestep - mDurations[mCurrentDuration].first()) %
        (mDurations[mCurrentDuration].period()) != 0)
        return;
    
    long numCells = 0;
    for (int ff = 0; ff < mFields.size(); ff++)
        numCells += mFields[ff]->calc(mDataFile);
    
//    LOG << "Timestep " << timestep << " wrote " << numCells << " cells.\n";
}

void OutputEHDB::
writeOutputHeader(ostream & file,
    string fileName,
    string dataFileName,
    Precision::Vec3 dxyz,
    Precision::Float dt,
    Precision::Vec3 origin) const
{
//    ptime now(second_clock::local_time());
    // Header: all IO description files will have this.
    file << "trogdor5data\n";
//    file << "trogdorMajorVersion " << TROGDOR_MAJOR_VERSION << "\n";
//    file << "trogdorSVNVersion NOTUSED\n";
    file << "trogdorBuildDate " << __DATE__ << "\n";
    if (sizeof(Precision::Float) == sizeof(float))
        file << "precision float32\n";
    else if (sizeof(Precision::Float) == sizeof(double))
        file << "precision float64\n";
    else
        throw(std::logic_error("Precision is not 32-bit or 64-bit."));
    file << "date " << TimeWrapper::timestampNow() << "\n";
    file << "dxyz " << dxyz << "\n";
    file << "dt " << dt << "\n";
    file << "origin " << setprecision(10) << origin << "\n";
//    file << "runlineDirection " << vp.lattice().runlineDirection() << "\n";
    file << "specfile " << OutputDirectory::path(fileName) << "\n";
    file << "datafile " << OutputDirectory::path(dataFileName) << "\n";
}

void OutputEHDB::
writeOutputFields(ostream & file, const OutputDescription & description) const
{
    for (int nn = 0; nn < description.fields().size(); nn++)
    {
        file << "field " << description.fields().at(nn).fieldAsString() << " "
            << description.fields().at(nn).cellOffset() << " "
            << description.fields().at(nn).timeOffset() << "\n";
    }
}

void OutputEHDB::
writeOutputRegions(ostream & file, 
    const vector<Region> & regions) const
{
    for (int nn = 0; nn < regions.size(); nn++)
    {
        file << "region "
            << regions[nn].yeeCells()
            << " stride "
            << regions[nn].stride();
        
        if (regions[nn].hasBounds())
        {
            file << " bounds " << regions[nn].bounds();
        }
        
        file << "\n";
    }
}

void OutputEHDB::
writeOutputDurations(ostream & file,
    const vector<Duration> & durations) const
{
    for (int nn = 0; nn < durations.size(); nn++)
    {
        file << "duration from "
            << durations[nn].first() << " to "
            << durations[nn].last() << " period "
            << durations[nn].period();
            
        if (durations[nn].hasInterval())
            file << " interval " << durations[nn].interval() << "\n";
        
        file << "\n";
    }
}






