/*
 *  OrientedVoxels.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 6/26/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "OrientedVoxels.h"
#include "RLEOperations.h"
#include "VoxelizePolyhedra.h"
#include "WriteSomeRLE.h"
#include "UserPreferences.h"
#include "Log.h"

using namespace std;
using namespace RLE;
using namespace SimpleMesh;
using namespace YeeUtilities;

template<class T>
static inline Matrix3<T> sDivide(const Matrix3<T> & matrix,
    double scalar)
{
    return (1.0 / scalar) * matrix;
}

static inline Vector3d sTimes(const Vector3d & v1, const Vector3d & v2)
{
    return v1*v2;
}

struct GetMatrixElement2
{
    GetMatrixElement2(int i, int j) : m_i(i), m_j(j) {}
    
    double operator()(const Matrix3d & mat)
    {
        return mat(m_i, m_j);
    }
    
    int m_i, m_j;
};

struct GetElement
{
    GetElement(int xyz) : mXYZ(xyz) {}
    double operator()(const Vector3d & v) const { return v[mXYZ]; }
    int mXYZ;
};

struct GetMatrixElement
{
    GetMatrixElement(int i, int j, int xyz) : m_i(i), m_j(j), m_xyz(xyz) {}
    double operator()(const Matrix3<Vector3d> & m) const
    {
        return m(m_i, m_j)[m_xyz];
    }
    int m_i, m_j, m_xyz;
};

struct NonzeroElement
{
    NonzeroElement(int xyz) { mXYZ = xyz; }
    bool operator()(const Vector3d & v) { return fabs(v[mXYZ]) > 1e-100; }
    int mXYZ;
};

bool sNonzero(const double & d)
{
    return (fabs(d) > 1e-100);
}



#pragma mark *** New stuff ***

OrientedVoxels::
OrientedVoxels() :
    mOrientationSmoothing(0),
    mFillFactorSmoothing(0)
{
}

OrientedVoxels::
OrientedVoxels(Rect3i voxelBounds, Rect3d realBounds) :
    mRealBounds(realBounds),
    mVoxelBounds(voxelBounds),
    mClipBounds(realBounds),
    mNonSymmetricDimensions(true,true,true),
    mOrientationSmoothing(0),
    mFillFactorSmoothing(0)
{
}

OrientedVoxels::
OrientedVoxels(Rect3i voxelBounds, Rect3d realBounds, Rect3d clipBounds) :
    mRealBounds(realBounds),
    mVoxelBounds(voxelBounds),
    mClipBounds(clipBounds),
    mNonSymmetricDimensions(true,true,true),
    mOrientationSmoothing(0),
    mFillFactorSmoothing(0)
{
}

static void checkFillFactorsBounded(const DynamicRLE3<double> & rle)
{
    DynamicRLE3<double>::ConstIterator itr;
    
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        const double SMALL = 1e-6;
        if (itr.mark() > 1 + SMALL || itr.mark() < -SMALL)
        {
            cerr << "Out-of-bounds fill factor (" << itr.mark() << ") at ";
            long x, y, z;
            rle.cartesianCoordinates(itr.position(), x, y, z);
            cerr << x << ", " << y << ", " << z << "\n";
            throw(std::logic_error("Fill factor out of bounds."));
        }
    }
}

void OrientedVoxels::
downsampleFillFactors(const vector<DynamicRLE3<double> > & halfCellFillFactors,
    int octant)
{
    mFillFactors.resize(halfCellFillFactors.size());
    
    for (int nn = 0; nn < halfCellFillFactors.size(); nn++)
    {
//        checkFillFactorsBounded(halfCellFillFactors.at(nn));
        DynamicRLE3<double> ff = averageDownsample2Clamp(
            halfCellFillFactors.at(nn),
            yeeToHalf(mVoxelBounds) + (halfCellOffset(octant) - 1),
            yeeToHalf(mVoxelBounds),
            mVoxelBounds);
        
        checkFillFactorsBounded(ff);
        
        if (fillFactorSmoothing() > 0)
        {
            cerr << "Smoothing fill factors (" << fillFactorSmoothing() << ")\n";
            mFillFactors.at(nn) = smoothClampN(ff,
                voxelBounds(), voxelBounds(), fillFactorSmoothing());
        }
        else
            mFillFactors.at(nn) = ff;
    }
    

}

void OrientedVoxels::
downsampleFillJacobians(
    const vector<map<unsigned int, DynamicRLE3<Vector3d> > > & halfCellJacobians,
    int octant)
{
    mFillFactorJacobians.resize(halfCellJacobians.size());
    
    map<unsigned int, DynamicRLE3<Vector3d> >::const_iterator itr;
    
    for (int nn = 0; nn < halfCellJacobians.size(); nn++)
    for (itr = halfCellJacobians[nn].begin();
        itr != halfCellJacobians[nn].end(); itr++)
    {
        DynamicRLE3<Vector3d> tooBigTemporary = sumDownsample2Clamp(
            itr->second,
            yeeToHalf(mVoxelBounds) + (halfCellOffset(octant) - 1),
            yeeToHalf(mVoxelBounds),
            mVoxelBounds);
        
        // weight after downsampling
        RLE::transform(tooBigTemporary, Vector3d(0.125, 0.125, 0.125),
            mFillFactorJacobians.at(nn)[itr->first], sTimes);
        
        if (fillFactorSmoothing() > 0)
        {
            mFillFactorJacobians.at(nn)[itr->first] = smoothClampN(
                mFillFactorJacobians.at(nn)[itr->first],
                voxelBounds(), voxelBounds(), fillFactorSmoothing());
        }
    }
}

Matrix3<Vector3d> sOuterProduct(const Matrix3d & m,
    const Vector3d & threeScalars)
{
    return Matrix3<Vector3d>(
        m[0]*threeScalars, m[1]*threeScalars, m[2]*threeScalars,
        m[3]*threeScalars, m[4]*threeScalars, m[5]*threeScalars,
        m[6]*threeScalars, m[7]*threeScalars, m[8]*threeScalars);
}

static Matrix3<Vector3d> sDivideBySquare(const Matrix3<Vector3d> & m,
    double d)
{
    double d2 = 1.0/(d*d);
    return Matrix3<Vector3d>(d2*m[0], d2*m[1], d2*m[2],
        d2*m[3], d2*m[4], d2*m[5], d2*m[6], d2*m[7], d2*m[8]);
}

// this horrid thing is needed for calling smoothClampN(), below.  Where oh
// where can I get a nice, clean result_of for C++????? :-( :-( :-(
template<class T>
Matrix3<Vector3<T> > operator*(double w, const Matrix3<Vector3<T> > & m)
{
    return Matrix3<Vector3<T> >(
        w*m[0], w*m[1], w*m[2],
        w*m[3], w*m[4], w*m[5],
        w*m[6], w*m[7], w*m[8]);
}


void OrientedVoxels::
initOrientations(const vector<std::vector<SimpleMesh::Triangle> > & meshTriangles)
{
//    cerr << "initOrientations()\n";
    DynamicRLE3<Matrix3d> unnormalizedMatrices;
    
    VoxelizePolyhedra vox(voxelBounds(), realBounds(), nonSymmetricDimensions(), clipBounds());
    
    vox.orientations(meshTriangles, unnormalizedMatrices);
    assert(hasCorrectDimensions(unnormalizedMatrices));
    
    if (!UserPreferences::defines("noNormalizedOrientation"))
    {
        DynamicRLE3<double> traces;
        
        if (orientationSmoothing() > 0)
        {
            cerr << "Smoothing orientations (" << orientationSmoothing() << ")\n";
            unnormalizedMatrices = smoothClampN(unnormalizedMatrices,
                voxelBounds(), voxelBounds(), orientationSmoothing());
            assert(hasCorrectDimensions(unnormalizedMatrices));
        }
        
        transform(unnormalizedMatrices, traces, trace<double>);
        assert(hasCorrectDimensions(traces));
        
        transform(unnormalizedMatrices, traces, mOrientation,
            sDivide<double>, Matrix3d(), 1.0);
        assert(hasCorrectDimensions(mOrientation));
        
        map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > > temp;
        vox.orientationJacobians(meshTriangles, temp);
        
        map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > >::const_iterator itr;
        if (orientationSmoothing() > 0)
        {
            for (itr = temp.begin(); itr != temp.end(); itr++)
            {
                assert(hasCorrectDimensions(itr->second));
                temp[itr->first] = smoothClampN(itr->second,
                    voxelBounds(), voxelBounds(), orientationSmoothing());
            }
        }
        
//        cerr << "Orientation Jacobian: " << temp.size() << " entries.\n";
        for (itr = temp.begin(); itr != temp.end(); itr++)
        {
            assert(hasCorrectDimensions(itr->second));
            DynamicRLE3<Vector3d> Dtraces;
            transform(itr->second, Dtraces, trace<Vector3d>);
            assert(hasCorrectDimensions(Dtraces));
//            cerr << "Dtraces: " << Dtraces.numRuns() << " runs.\n";
            
            // DP = DM/trM - MtrDM/trM/trM
            DynamicRLE3<Matrix3<Vector3d> > DM_trM;
            
            DynamicRLE3<double> tracesMasked;
            restriction(traces, itr->second, tracesMasked);
            assert(hasCorrectDimensions(tracesMasked));
            
            transform(itr->second, tracesMasked, DM_trM,
                sDivide<Vector3d>, Matrix3<Vector3d>(), 1.0);
            assert(hasCorrectDimensions(DM_trM));
//            cerr << "DM_trM: " << DM_trM.numRuns() << " runs.\n";
            
            DynamicRLE3<Matrix3d> unnormalizedMatricesMasked;
            restriction(unnormalizedMatrices, Dtraces, unnormalizedMatricesMasked);
            assert(hasCorrectDimensions(unnormalizedMatricesMasked));
            
            DynamicRLE3<Matrix3<Vector3d> > MtrDM;
            transform(unnormalizedMatricesMasked, Dtraces, MtrDM, sOuterProduct,
                Matrix3d(), Vector3d());
            assert(hasCorrectDimensions(MtrDM));
//            cerr << "MtrDM: " << MtrDM.numRuns() << " runs.\n";
            
            DynamicRLE3<Matrix3<Vector3d> > MtrDM_trM_trM;
            transform(MtrDM, tracesMasked, MtrDM_trM_trM, sDivideBySquare,
                Matrix3<Vector3d>(), 1.0);
            assert(hasCorrectDimensions(MtrDM_trM_trM));
            
//            cerr << "MtrDM_trM_trM: " << MtrDM_trM_trM.numRuns() << " runs.\n";
            
            mOrientationJacobians[itr->first] = DM_trM - MtrDM_trM_trM;
            assert(hasCorrectDimensions(mOrientationJacobians[itr->first]));
            assert(hasCorrectDimensions(mOrientationJacobians[itr->first]));
        }
    }
    else // this little block is really just for debugging purposes
    {
        LOG << "Not normalizing orientation.\n";
        vox.orientationJacobians(meshTriangles, mOrientationJacobians);
        mOrientation = unnormalizedMatrices;
    }
}



//void OrientedVoxels::
//initOrientations(const vector<TrackedPolyhedron*> & polyhedra)
//{
////    cerr << "initOrientations()\n";
//    DynamicRLE3<Matrix3d> unnormalizedMatrices;
//
//    VoxelizePolyhedra vox(voxelBounds(), realBounds(), nonSymmetricDimensions(), clipBounds());
//
//    vox.orientations(polyhedra, unnormalizedMatrices);
//
//    if (!UserPreferences::defines("noNormalizedOrientation"))
//    {
//        DynamicRLE3<double> traces;
//
//        if (orientationSmoothing() > 0)
//        {
//            cerr << "Smoothing orientations (" << orientationSmoothing() << ")\n";
//            unnormalizedMatrices = smoothClampN(unnormalizedMatrices,
//                voxelBounds(), voxelBounds(), orientationSmoothing());
//        }
//
//        transform(unnormalizedMatrices, traces, trace<double>);
//
//        transform(unnormalizedMatrices, traces, mOrientation,
//            sDivide<double>, Matrix3d(), 1.0);
//
//        map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > > temp;
//        vox.orientationJacobians(polyhedra, temp);
//
//        map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > >::const_iterator itr;
//        if (orientationSmoothing() > 0)
//        {
//            for (itr = temp.begin(); itr != temp.end(); itr++)
//            {
//                temp[itr->first] = smoothClampN(itr->second,
//                    voxelBounds(), voxelBounds(), orientationSmoothing());
//            }
//        }
//
////        cerr << "Orientation Jacobian: " << temp.size() << " entries.\n";
//        for (itr = temp.begin(); itr != temp.end(); itr++)
//        {
//            DynamicRLE3<Vector3d> Dtraces;
//            transform(itr->second, Dtraces, trace<Vector3d>);
////            cerr << "Dtraces: " << Dtraces.numRuns() << " runs.\n";
//
//            // DP = DM/trM - MtrDM/trM/trM
//            DynamicRLE3<Matrix3<Vector3d> > DM_trM;
//
//            DynamicRLE3<double> tracesMasked;
//            restriction(traces, itr->second, tracesMasked);
//
//            transform(itr->second, tracesMasked, DM_trM,
//                sDivide<Vector3d>, Matrix3<Vector3d>(), 1.0);
////            cerr << "DM_trM: " << DM_trM.numRuns() << " runs.\n";
//
//            DynamicRLE3<Matrix3d> unnormalizedMatricesMasked;
//            restriction(unnormalizedMatrices, Dtraces, unnormalizedMatricesMasked);
//
//            DynamicRLE3<Matrix3<Vector3d> > MtrDM;
//            transform(unnormalizedMatricesMasked, Dtraces, MtrDM, sOuterProduct,
//                Matrix3d(), Vector3d());
////            cerr << "MtrDM: " << MtrDM.numRuns() << " runs.\n";
//
//            DynamicRLE3<Matrix3<Vector3d> > MtrDM_trM_trM;
//            transform(MtrDM, tracesMasked, MtrDM_trM_trM, sDivideBySquare,
//                Matrix3<Vector3d>(), 1.0);
//
////            cerr << "MtrDM_trM_trM: " << MtrDM_trM_trM.numRuns() << " runs.\n";
//
//            mOrientationJacobians[itr->first] = DM_trM - MtrDM_trM_trM;
//        }
//    }
//    else // this little block is really just for debugging purposes
//    {
//        LOG << "Not normalizing orientation.\n";
//        vox.orientationJacobians(polyhedra, mOrientationJacobians);
//        mOrientation = unnormalizedMatrices;
//    }
//}


// TODO: see if I can stuff this into the sensitiveCells() function.
struct HasNonzeroElement
{
    HasNonzeroElement(int i, int j) : m_i(i), m_j(j) {}
    
    bool operator()(const Matrix3<Vector3d> & matrix) const
    {
        return matrix(m_i,m_j)[0] != 0 ||
            matrix(m_i,m_j)[1] != 0 ||
            matrix(m_i,m_j)[2] != 0;
    }
    
    int m_i, m_j;
};

RLE::SupportRegion3 OrientedVoxels::
sensitiveCells(int i, int j) const
{
//    SupportRegion3 out(mVoxelBounds.size().asArray());
    SupportRegion3 out(mNonSymmetricDimensions.asArray());
    
    map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > >::const_iterator itr;
    for (itr = mOrientationJacobians.begin();
        itr != mOrientationJacobians.end(); itr++)
    {
        assert(hasCorrectDimensions(itr->second));
        SupportRegion3 temp;
        transform(itr->second, temp, HasNonzeroElement(i,j));
        assert(hasCorrectDimensions(temp));
        out += temp;
    }
//    cerr << "----- Orientation:\n" << out << "\n";
    
    SupportRegion3 fillTemp(mVoxelBounds.size().asArray());
    
    map<unsigned int, DynamicRLE3<Vector3d> >::const_iterator itr2;
    for (int nn = 0; nn < mFillFactorJacobians.size(); nn++)
    for (itr2 = mFillFactorJacobians[nn].begin();
        itr2 != mFillFactorJacobians[nn].end(); itr2++)
    {
        assert(hasCorrectDimensions(itr2->second));
        SupportRegion3 temp;
        transform(itr2->second, temp, Vector3d::IsNonzero());
        assert(hasCorrectDimensions(temp));
        out += temp;
        fillTemp += temp;
    }
    
//    cerr << "----- Fill:\n" << fillTemp << "\n";
    
    assert(hasCorrectDimensions(out));
    return out;
}

const std::vector<RLE::DynamicRLE3<double> > & OrientedVoxels::
fillFactors() const
{
    for (int ii = 0; ii < mFillFactors.size(); ii++)
        assert(hasCorrectDimensions(mFillFactors.at(ii)));
    return mFillFactors;
}

const RLE::DynamicRLE3<double> & OrientedVoxels::
fillFactors(int material) const
{
    assert(hasCorrectDimensions(mFillFactors.at(material)));
    return mFillFactors.at(material);
}

const RLE::DynamicRLE3<Matrix3d> & OrientedVoxels::
orientations() const
{
    assert(hasCorrectDimensions(mOrientation));
    return mOrientation;
}

static bool sIsNonzero(const double & val)
{
    return (val != 0);
}

RLE::DynamicRLE3<double> OrientedVoxels::
orientations(int i, int j) const
{
    RLE::DynamicRLE3<double> out1;
    transform(mOrientation, out1, Matrix3d::GetElement(i,j));
    
    RLE::DynamicRLE3<double> outNoNonzeros;
    filter(out1, outNoNonzeros, sIsNonzero);
    
    assert(hasCorrectDimensions(outNoNonzeros));
    
    return outNoNonzeros;
}

const std::vector<std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > > &
OrientedVoxels::
fillFactorJacobians() const
{
    return mFillFactorJacobians;
}

const map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > > & OrientedVoxels::
orientationJacobians() const
{
    map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > >::const_iterator itr;
    for (itr = mOrientationJacobians.begin(); itr != mOrientationJacobians.end();
        itr++)
    {
        if (!hasCorrectDimensions(itr->second))
        {
            throw std::runtime_error("Wrong dimensions!");
        }
//        assert(hasCorrectDimensions(itr->second));
    }
    
    return mOrientationJacobians;
}

std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > OrientedVoxels::
orientationJacobians(int i, int j) const
{
    map<unsigned int, RLE::DynamicRLE3<Vector3d> > out;
    
    map<unsigned int, RLE::DynamicRLE3<Matrix3<Vector3d> > >::const_iterator itr;
    for (itr = mOrientationJacobians.begin();
        itr != mOrientationJacobians.end(); itr++)
    {
        assert(hasCorrectDimensions(itr->second));
        RLE::DynamicRLE3<Vector3d> temp;
        transform(itr->second, temp, Matrix3<Vector3d>::GetElement(i,j));
        assert(hasCorrectDimensions(temp));
        filter(temp, out[itr->first], Vector3d::IsNonzero());
        
        assert(hasCorrectDimensions(out[itr->first]));
    }
    
    return out;
}

long OrientedVoxels::
bytes() const
{
    long fillBytes = 0;
    long orientBytes = 0;
    long dFillBytes = 0;
    long dOrientBytes = 0;
    
    for (int nn = 0; nn < mFillFactors.size(); nn++)
        fillBytes += mFillFactors.at(nn).bytes();
    
    orientBytes += mOrientation.bytes();
    
    for (int nn = 0; nn < mFillFactorJacobians.size(); nn++)
    {
        map<unsigned int, DynamicRLE3<Vector3d> >::const_iterator itr;
        for (itr = mFillFactorJacobians[nn].begin();
            itr != mFillFactorJacobians[nn].end(); itr++)
        {
            dFillBytes += itr->second.bytes();
        }
    }
    
    map<unsigned int, DynamicRLE3<Matrix3<Vector3d> > >::const_iterator itr;
    for (itr = mOrientationJacobians.begin();
        itr != mOrientationJacobians.end(); itr++)
    {
        dOrientBytes += itr->second.bytes();
    }
    
    LOG << "Fill factors " << fillBytes/1024.0 << " kB\n";
    LOG << "Orientation " << orientBytes/1024.0 << " kB\n";
    LOG << "Fill factor Jacobians " << dFillBytes/1024.0 << " kB\n";
    LOG << "Orientation Jacobians " << dOrientBytes/1024.0 << " kB\n";
    
    return fillBytes + orientBytes + dFillBytes + dOrientBytes;
}



