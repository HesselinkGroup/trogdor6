/*
 *  FillFactorReport.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 2/8/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef T6_FILLFACTORREPORT
#define T6_FILLFACTORREPORT

#include <fstream>
#include <vector>
#include "geometry.h"

class OrientedVoxels;

class FillFactorReport
{
public:
    FillFactorReport();
    FillFactorReport(std::string fileName, Rect3i cells);
    
    void open(std::string fileName, Rect3i cells);
    void write(const OrientedVoxels & vox);
    void close();

private:
    int mNumArrays;
    int mNumMaterials;
    Rect3i mCells;
    std::string mFileName;
    std::ofstream mDataFile;
};

class FillFactorSensitivityReport
{
public:
    FillFactorSensitivityReport();
    FillFactorSensitivityReport(std::string fileName, Rect3i cells);
    
    void open(std::string fileName, Rect3i cells);
    void write(const OrientedVoxels & vox, int octant);
    void close();

private:
    Rect3i mCells;
    std::ofstream mFile;
};



#endif
