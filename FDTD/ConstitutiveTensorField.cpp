/*
 *  ConstitutiveTensorField.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/23/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "ConstitutiveTensorField.h"
#include "YeeUtilities.h"
#include "Log.h"
//#include "calc.hh"
#include "RLEOperations.h"
#include "UserPreferences.h"
#include "OutputDirectory.h"

#include <map>
#include <fstream>

using namespace YeeUtilities;
using namespace RLE;
using namespace std;
namespace P = Precision;

ConstitutiveTensorField::
ConstitutiveTensorField() :
    mOctant(Octant::electric())
{
}

ConstitutiveTensorField::
ConstitutiveTensorField(const Octant & o, Rect3i bounds) :
    mOctant(o),
    mBounds(bounds)
{
//    setDimensions(bounds);
}

//void ConstitutiveTensorField::
//setDimensions(Rect3i yeeCells)
//{
//    bool hasX = (yeeCells.num(0) != 1);
//    bool hasY = (yeeCells.num(1) != 1);
//    bool hasZ = (yeeCells.num(2) != 1);
//    for (int xyz = 0; xyz < 3; xyz++)
//    for (int dd = 0; dd < 3; dd++)
//        (*this)(xyz,dd).dimensions(hasX, hasY, hasZ);
//}

void ConstitutiveTensorField::
extrude(Rect3i innerHalfCells, Rect3i outerHalfCells)
{
    // Diagonal terms
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Rect3i inner = YeeUtilities::halfToYee(innerHalfCells, tensorOctant(xyz,xyz));
        Rect3i outer = YeeUtilities::halfToYee(outerHalfCells, tensorOctant(xyz,xyz));
        extrude_slow((*this)(xyz, xyz), inner, outer);
    }
    
    // Off-diagonal terms
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Rect3i inner = YeeUtilities::halfToYee(innerHalfCells,
            tensorOctant(xyz,(xyz+1)%3));
        Rect3i outer = YeeUtilities::halfToYee(outerHalfCells,
            tensorOctant(xyz,(xyz+1)%3));
        extrude_slow((*this)(xyz, (xyz+1)%3), inner, outer);
    }
}

void ConstitutiveTensorField::
copy(const ConstitutiveTensorField & source, Rect3i copyCells,
    Rect3i pasteCells)
{
    assert(copyCells.size() == pasteCells.size());
    Vector3i shift = pasteCells.p1 - copyCells.p1;
    
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    {
        for (int zz = copyCells.p1[2]; zz <= copyCells.p2[2]; zz++)
        for (int yy = copyCells.p1[1]; yy <= copyCells.p2[1]; yy++)
        for (int xx = copyCells.p1[0]; xx <= copyCells.p2[0]; xx++)
        if (source.mComponents[ii][jj].isMarked(xx,yy,zz))
        {
            mComponents[ii][jj].mark(xx+shift[0], yy+shift[1], zz+shift[2],
                source.mComponents[ii][jj](xx,yy,zz));
        }
    }
}

SupportRegion3 ConstitutiveTensorField::
support(int i, int j, int numeratorOrder, int denominatorOrder) const
{
    const DynamicRLE3<P::RationalFunction> & rle = (*this)(i,j);
    
    SupportRegion3 out(rle.hasDimension(0),
        rle.hasDimension(1),
        rle.hasDimension(2));
    transform(rle, out, IsOrder(numeratorOrder, denominatorOrder));
    return out;
}


DynamicRLE3<P::RationalFunction> ConstitutiveTensorField::
filter(int i, int j, int numeratorOrder, int denominatorOrder) const
{
    const DynamicRLE3<P::RationalFunction> & rle = (*this)(i,j);
    
    DynamicRLE3<P::RationalFunction> out;
    transform(rle, out, SelectOrder(numeratorOrder, denominatorOrder));
    return out;
}

DynamicRLE3<P::RationalFunction> & ConstitutiveTensorField::
operator()(int i, int j)
{
    assert(i >= 0 && i < 3 && j >= 0 && j < 3);
    return mComponents[i][j];
//    if (i == j)
//        return mDiagonal[i];
//    else
//        return mOffDiagonal[(i+j)%3];
}

const DynamicRLE3<P::RationalFunction> & ConstitutiveTensorField::
operator()(int i, int j) const
{
    assert(i >= 0 && i < 3 && j >= 0 && j < 3);
    return mComponents[i][j];
}

int ConstitutiveTensorField::
tensorOctant(int i, int j) const
{
    assert(i >= 0 && i < 3 && j >= 0 && j < 3);
    
    if (i == j)
        return mOctant(i, j);
    else
    {
        string oct = UserPreferences::valueAs<string>("offDiagonalOctant");
        
        if (oct == "i")
            return mOctant(i, i);
        else if (oct == "j")
            return mOctant(j, j);
        else if (oct == "0")
            return 0;
        else if (oct == "7")
            return 7;
        else
            throw(std::logic_error(string("Invalid tensor position ")+oct));
    }
}

double ConstitutiveTensorField::
bytes() const
{
    double total = 0.0;
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
        total += mComponents[ii][jj].bytes();
    return total;
}

void ConstitutiveTensorField::
writeBinary(std::string fileName, Rect3i cells) const
{
    ofstream description(OutputDirectory::path(fileName + ".txt").c_str());
    description << "trogdor6data\n";
    description << "trogdorBuildDate " <<__DATE__ << "\n";
    if (sizeof(Precision::Float) == sizeof(float))
        description << "precision float32\n";
    else if (sizeof(Precision::Float) == sizeof(double))
        description << "precision float64\n";
    else
        throw(std::logic_error("Precision is not 32-bit or 64-bit."));
    description << "specfile " << OutputDirectory::path(fileName + ".txt") << "\n";
    description << "datafile " << OutputDirectory::path(fileName) << "\n";
    description << "region " << cells << " stride "
        << Precision::Vec3(1,1,1) << "\n";
    
    //const int MAX_ORDER = Precision::RationalFunction::maximumOrder();
    const int MAX_ORDER = 10;
    
    int maxNumerOrder, maxDenomOrder;
    
    ofstream data(fileName.c_str(), ios::binary);
    
    // Permittivity numerator, which is the denominator of 1/eps:
    for (int lag = 0; lag <= MAX_ORDER; lag++)
    {
        //for (int xyz = 0; xyz < 3; xyz++)
        for (int jj = 0; jj < 3; jj++)
        for (int ii = 0; ii < 3; ii++)
        {
            DynamicRLE3<Precision::Float> term;
            //transform(mDiagonal[xyz], term, DenominatorTerm(lag)); // careful!
            transform((*this)(ii,jj), term, DenominatorTerm(lag));
            
//            LOG << "My runs: " << (*this)(ii,jj).numRuns() << ", now "
//                << term.numRuns() << " runs to write.\n";
            
            bool theFlag = 0;
            //if (term.numRuns() > 0)
            {
                //maxNumerOrder = lag;
                for (int zz = cells.p1[2]; zz <= cells.p2[2]; zz++)
                for (int yy = cells.p1[1]; yy <= cells.p2[1]; yy++)
                for (int xx = cells.p1[0]; xx <= cells.p2[0]; xx++)
                {
                    Precision::Float val = term.at(xx,yy,zz);
                    data.write((char*)&val, sizeof(Precision::Float));
                    
                    if (val != 0.0)
                        theFlag = 1;
                }
            }
            
            if (theFlag)
            {
                LOG << "Denom nonzed: (" << ii << ", " << jj << ") lag "
                    << lag << "\n";
            }
        }
    }
    
    // Permittivity denominator, which is the numerator of 1/eps:
    for (int lag = 0; lag <= MAX_ORDER; lag++)
    {
        //for (int xyz = 0; xyz < 3; xyz++)
        for (int jj = 0; jj < 3; jj++)
        for (int ii = 0; ii < 3; ii++)
        {
            DynamicRLE3<Precision::Float> term;
            //transform(mDiagonal[xyz], term, NumeratorTerm(lag)); // careful!
            transform((*this)(ii,jj), term, NumeratorTerm(lag));
            
            bool theFlag = 0;
            //if (term.numRuns() > 0)
            {
                //maxDenomOrder = lag;
                for (int zz = cells.p1[2]; zz <= cells.p2[2]; zz++)
                for (int yy = cells.p1[1]; yy <= cells.p2[1]; yy++)
                for (int xx = cells.p1[0]; xx <= cells.p2[0]; xx++)
                {
                    Precision::Float val = term.at(xx,yy,zz);
                    data.write((char*)&val, sizeof(Precision::Float));
                    if (val != 0.0)
                        theFlag = 1;
                }
            }
            if (theFlag)
            {
                LOG << "Numer nonzed: (" << ii << ", " << jj << ") lag "
                    << lag << "\n";
            }
        }
    }
    data.close();
    
    maxNumerOrder = MAX_ORDER;
    maxDenomOrder = MAX_ORDER;
    description << "numerator " << maxNumerOrder+1 << " lags\n";
    description << "denominator " << maxDenomOrder+1 << " lags\n";
    
    description.close();
}





