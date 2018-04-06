/*
 *  OrientationReport.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 11/12/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "OrientationReport.h"
#include "OrientedVoxels.h"
#include "WriteSomeRLE.h"
#include "OutputDirectory.h"

#include <sstream>


using namespace std;
using namespace RLE;

OrientationReport::
OrientationReport() :
    mNumArrays(0)
{
}

OrientationReport::
OrientationReport(string fileName, Rect3i cells) :
    mNumArrays(0),
    mCells(cells),
    mFileName(fileName),
    mDataFile(fileName.c_str(), ios::binary)
{
}

void OrientationReport::
open(string fileName, Rect3i cells)
{
    mNumArrays = 0;
    mCells = cells;
    mFileName = fileName;
    mDataFile.open(OutputDirectory::path(fileName).c_str(), ios::binary);
}

void OrientationReport::
write(const OrientedVoxels & vox)
{
    mNumArrays++;
    
    for (int jj = 0; jj < 3; jj++) // tensor elements
    for (int ii = 0; ii < 3; ii++)
    for (int zz = mCells.p1[2]; zz <= mCells.p2[2]; zz++) // spatial indices
    for (int yy = mCells.p1[1]; yy <= mCells.p2[1]; yy++)
    for (int xx = mCells.p1[0]; xx <= mCells.p2[0]; xx++)
    {
        Precision::Float val = vox.orientations().at(xx,yy,zz)(ii,jj);
        mDataFile.write((char*)&val, sizeof(Precision::Float));
    }
}

void OrientationReport::
close()
{
    mDataFile.close();
    
    ofstream txt(OutputDirectory::path(mFileName + ".txt").c_str());
    
    txt << "trogdor6data\n";
    txt << "trogdorBuildDate " << __DATE__ << "\n";
    if (sizeof(Precision::Float) == sizeof(float))
        txt << "precision float32\n";
    else if (sizeof(Precision::Float) == sizeof(double))
        txt << "precision float64\n";
    else
        throw(std::logic_error("Precision is not 32-bit or 64-bit."));
    txt << "dims " << mCells.num(0) << " " << mCells.num(1)
        << " " << mCells.num(2) << " 3 3";
    
    if (mNumArrays == 3)
        txt << " 3\n";
    else if (mNumArrays == 9)
        txt << " 3 3\n";
    else
        throw(std::logic_error("Should be 3 or 9"));
    
    txt.close();
}


OrientationSensitivityReport::
OrientationSensitivityReport()
{
}

OrientationSensitivityReport::
OrientationSensitivityReport(string fileName, Rect3i cells)
{
    open(fileName, cells);
}

void OrientationSensitivityReport::
open(string fileName, Rect3i cells)
{
    mCells = cells;
    mFile.open(OutputDirectory::path(fileName).c_str());
    
    mFile << "function f = orientationSensitivity()\n";
    mFile << "% f{vertex, epsilon row, epsilon column, octant} = "
        "[d/dvx, d/dvy, d/dvz]\n";
    
    mFile << setprecision(16);
}

void OrientationSensitivityReport::
write(const OrientedVoxels & vox, int row, int col, int octant)
{
    map<unsigned int, DynamicRLE3<Vector3d> > oj = 
        vox.orientationJacobians(row, col);
    
    map<unsigned int, DynamicRLE3<Vector3d> >::const_iterator itr;
    for (itr = oj.begin(); itr != oj.end(); itr++)
    {
        ostringstream varName;
        varName << "f{" << itr->first+1 << ","
            << row+1 << "," << col+1 << "," << octant+1 << "}";
        matlabRLE(itr->second, mCells, varName.str(), mFile);
    }
}

void OrientationSensitivityReport::
close()
{
    mFile.close();
}


