/*
 *  FillFactorReport.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 2/8/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "FillFactorReport.h"
#include "OrientedVoxels.h"
#include "WriteSomeRLE.h"
#include "OutputDirectory.h"

#include <sstream>

using namespace std;
using namespace RLE;

FillFactorReport::
FillFactorReport() :
    mNumArrays(0),
    mNumMaterials(0)
{
}

FillFactorReport::
FillFactorReport(string fileName, Rect3i cells) :
    mNumArrays(0),
    mNumMaterials(0),
    mCells(cells),
    mFileName(fileName),
    mDataFile(fileName.c_str(), ios::binary)
{
}

void FillFactorReport::
open(string fileName, Rect3i cells)
{
    mNumArrays = 0;
    mNumMaterials = 0;
    mCells = cells;
    mFileName = fileName;
    mDataFile.open(OutputDirectory::path(fileName).c_str(), ios::binary);
}

void FillFactorReport::
write(const OrientedVoxels & vox)
{
    mNumArrays++;
    mNumMaterials = vox.fillFactors().size();
    
    for (int ff = 0; ff < vox.fillFactors().size(); ff++)
    for (int zz = mCells.p1[2]; zz <= mCells.p2[2]; zz++)
    for (int yy = mCells.p1[1]; yy <= mCells.p2[1]; yy++)
    for (int xx = mCells.p1[0]; xx <= mCells.p2[0]; xx++)
    {
        Precision::Float val = vox.fillFactors()[ff].at(xx,yy,zz);
        mDataFile.write((char*)&val, sizeof(Precision::Float));
    }
}

void FillFactorReport::
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
        << " " << mCells.num(2) << " " << mNumMaterials;
    
    if (mNumArrays == 3)
        txt << " 3\n";
    else if (mNumArrays == 9)
        txt << " 3 3\n";
    else
        throw(std::logic_error("Strange number of arrays, nether xyz nor ij."));
    
    txt.close();
}

FillFactorSensitivityReport::
FillFactorSensitivityReport()
{
}

FillFactorSensitivityReport::
FillFactorSensitivityReport(string fileName, Rect3i cells)
{
    open(fileName, cells);
}

void FillFactorSensitivityReport::
open(string fileName, Rect3i cells)
{
    mCells = cells;
    mFile.open(OutputDirectory::path(fileName).c_str());
    
    mFile << "function f = fillFactorSensitivity()\n";
    mFile << "% f{octant, material, vertex} = [d/dvx, d/dvy, d/dvz]\n";
    mFile << setprecision(16);
}

void FillFactorSensitivityReport::
write(const OrientedVoxels & vox, int octant)
{
    map<unsigned int, RLE::DynamicRLE3<Vector3d> >::const_iterator itr;
    
    for (int nn = 0; nn < vox.fillFactorJacobians().size(); nn++)
    for (itr = vox.fillFactorJacobians().at(nn).begin();
        itr != vox.fillFactorJacobians().at(nn).end(); itr++)
    {
        ostringstream varName;
        varName << "f{" << octant+1 << "," << nn+1 << "," << itr->first+1 << "}";
        matlabRLE(itr->second, mCells, varName.str(), mFile);
    }
}

void FillFactorSensitivityReport::
close()
{
    mFile.close();
}


