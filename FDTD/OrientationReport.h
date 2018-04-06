/*
 *  OrientationReport.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 11/12/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _ORIENTATIONREPORT_
#define _ORIENTATIONREPORT_

#include <fstream>
#include <vector>
#include "geometry.h"

class OrientedVoxels;

class OrientationReport
{
public:
    OrientationReport();
    OrientationReport(std::string fileName, Rect3i cells);
    
    void open(std::string fileName, Rect3i cells);
    void write(const OrientedVoxels & vox);
    void close();

private:
    int mNumArrays;
    Rect3i mCells;
    std::string mFileName;
    std::ofstream mDataFile;
};

class OrientationSensitivityReport
{
public:
    OrientationSensitivityReport();
    OrientationSensitivityReport(std::string fileName, Rect3i cells);
    
    void open(std::string fileName, Rect3i cells);
    void write(const OrientedVoxels & vox, int tensorRow, int tensorColumn,
        int octant);
    void close();

private:
    Rect3i mCells;
    std::ofstream mFile;
};



#endif
