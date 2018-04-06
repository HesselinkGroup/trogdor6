/*
 *  RowVoxelizer.cpp
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 7/5/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "RowVoxelizer.h"
#include "ClipTriangle.h"

using namespace std;
using namespace SimpleMesh;

ClippedTriangle RowVoxelizer::
clipToRow(const Triangle & tri, int jj, int kk) const
{
    int u0 = mRowAxis;
    int u1 = (u0+1)%3;
    int u2 = (u1+1)%3;
    
    Rect2d pb = pixelBounds(jj, kk);
    
//    cerr << "pb = " << pb << "\n";
    
    ClippedTriangle clipped(tri);
    clipped = clipped.partAbove(u1, pb.p1[0]);
    clipped = clipped.partBelow(u1, pb.p2[0]);
    clipped = clipped.partAbove(u2, pb.p1[1]);
    clipped = clipped.partBelow(u2, pb.p2[1]);
    
    return clipped;
}

ClippedTriangle RowVoxelizer::
clipToRow(const Triangle & tri, int jj, int kk,
    vector<SensitiveEdge> & edgeBuffer) const
{
    int u0 = mRowAxis;
    int u1 = (u0+1)%3;
    int u2 = (u1+1)%3;
    
    Rect3d triBounds = tri.triangle().bounds();
    
    Rect2d pb = pixelBounds(jj, kk);
    
//    cerr << "pb = " << pb << "\n";
    
    Plane3d planeLow1(Vector3d::unit(u1), -pb.p1[0]),
        planeHigh1(Vector3d::unit(u1), -pb.p2[0]),
        planeLow2(Vector3d::unit(u2), -pb.p1[1]),
        planeHigh2(Vector3d::unit(u2), -pb.p2[1]);
    
//    cerr << "Clip " << tri << "\n";
//    cerr << "Planes:\n"
//        << planeLow1 << "\n" << planeHigh1 << "\n"
//        << planeLow2 << "\n" << planeHigh2 << "\n";
        
    ClippedTriangle clipped(tri);
    
    
    if (fabs(clipped.triangle()->unitNormal()[u1]) < 0.999999)
    {
        if (triBounds.p1[u1] < pb.p1[0] && triBounds.p2[u1] > pb.p1[0])
        {
            edgeBuffer[0] = SensitiveEdge(tri, planeLow1);
            clipped = clipped.partAbove(u1, pb.p1[0], &edgeBuffer[0]);
        }
        else
            clipped = clipped.partAbove(u1, pb.p1[0]);

        if (triBounds.p1[u1] < pb.p2[0] && triBounds.p2[u1] > pb.p2[0])
        {
            edgeBuffer[1] = SensitiveEdge(tri, planeHigh1);
            clipped = clipped.partBelow(u1, pb.p2[0], &edgeBuffer[1]);
        }
        else
            clipped = clipped.partBelow(u1, pb.p2[0]);
    }
    else
    {
        clipped = clipped.partAbove(u1, pb.p1[0]);
        clipped = clipped.partBelow(u1, pb.p2[0]);
    }
    
    if (fabs(clipped.triangle()->unitNormal()[u2]) < 0.999999)
    {
        if (triBounds.p1[u2] < pb.p1[1] && triBounds.p2[u2] > pb.p1[1])
        {
            edgeBuffer[2] = SensitiveEdge(tri, planeLow2);
            clipped = clipped.partAbove(u2, pb.p1[1], &edgeBuffer[2]);
        }
        else
            clipped = clipped.partAbove(u2, pb.p1[1]);
        
        if (triBounds.p1[u2] < pb.p2[1] && triBounds.p2[u2] > pb.p2[1])
        {
            edgeBuffer[3] = SensitiveEdge(tri, planeHigh2);
            clipped = clipped.partBelow(u2, pb.p2[1], &edgeBuffer[3]);
        }
        else
            clipped = clipped.partBelow(u2, pb.p2[1]);
    }
    else
    {
        clipped = clipped.partAbove(u2, pb.p1[1]);
        clipped = clipped.partBelow(u2, pb.p2[1]);
    }
    
    return clipped;
}

#pragma mark *** FillFactorRowVoxelizer ***


static void checkFillFactorsBounded(const RLE::DynamicRLE3<double> & rle, double vol)
{
    RLE::DynamicRLE3<double>::ConstIterator itr;
    
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        const double SMALL = 1e-6;
        if (itr.mark()/vol > 1 + SMALL || itr.mark()/vol < -SMALL)
        {
            cerr << "Out-of-bounds fill factor (" << itr.mark() << ") at ";
            long x, y, z;
            rle.cartesianCoordinates(itr.position(), x, y, z);
            cerr << x << ", " << y << ", " << z << "\n";
            //throw(std::logic_error("Fill factor out of bounds."));
        }
    }
}

FillFactorRowVoxelizer::
FillFactorRowVoxelizer(int rowAxis, Rect3d bounds, Rect3i voxelBounds, Vector3b nonSymmetricDimensions,
    RLE::DynamicRLE3<double> & outFillFactors) :
    RowVoxelizer(rowAxis, bounds, voxelBounds, nonSymmetricDimensions),
    mFillFactors(outFillFactors)
{
}

void FillFactorRowVoxelizer::
voxelize(int jj, int kk, const vector<RayIntersection> & intersections,
    int numHits)
{
    // 1.  Process all the triangles and make a hit list
    
    vector<FillRecord> fillRecords;
    
    for (int ii = 0; ii < numHits; ii++)
        voxelizeOneTri(intersections[ii].triangle(), jj, kk, fillRecords);
    
    sort(fillRecords.begin(), fillRecords.end());
    
    handleSortedFillRecords(jj, kk, fillRecords);
    
}


void FillFactorRowVoxelizer::
voxelizeOneTri(const Triangle* tri, int jj, int kk,
    std::vector<FillRecord> & fillRecords)
{
    ClippedTriangle clipped = clipToRow(*tri, jj, kk);
    if (clipped.vertices().size() < 3)
        return;
    
    Rect3d clippedBounds = clipped.bounds();
    
    double xmin = clippedBounds.p1[rowAxis()];
    double xmax = clippedBounds.p2[rowAxis()];
    
    int minVoxel = voxel(xmin);
    int maxVoxel = voxel(xmax);
    
    for (int vv = minVoxel; vv <= maxVoxel; vv++)
    {
        ClippedTriangle polyInVoxel = clipped.partAbove(rowAxis(), voxelMin(vv));
        polyInVoxel = polyInVoxel.partBelow(rowAxis(), voxelMin(vv+1));
        
//        if (polyInVoxel.vertices().size() < 3 && jj == -1 && kk == -1)
//            int blarg = 3;
        
        FillRecord fr(vv, polyInVoxel.area2d(rowAxis()),
            polyInVoxel.volume(rowAxis(), voxelMin(vv)));
        
//        if (jj == -1 && kk == -1)
//        {
//            cerr << "Poly " << polyInVoxel << "\n";
//            cerr << "FR " << fr << "\n";
//        }
        
        fillRecords.push_back(fr);
    }
}

void FillFactorRowVoxelizer::
handleSortedFillRecords(int jj, int kk, const vector<FillRecord> & fillRecords)
{
    if (fillRecords.empty())
        return;
    
    const double cellHeight = voxelSize()[rowAxis()];
    
    int currentCell;
    int lastCell = fillRecords[0].voxel-1;
    double floorArea = 0.0;
    int currentRecord = 0;
    
    
    while (currentRecord < fillRecords.size())
    {
        double volume = 0.0, newArea = 0.0;
        currentCell = fillRecords[currentRecord].voxel;
        fill(lastCell+1, currentCell-1, jj, kk, floorArea*cellHeight);
        lastCell = currentCell;
        
        while( currentRecord < fillRecords.size() &&
            fillRecords.at(currentRecord).voxel == currentCell )
        {
            volume += fillRecords[currentRecord].volume;
            newArea += fillRecords[currentRecord].area;
            currentRecord++;
        }
        
        double cellVol = volume + (floorArea-newArea)*cellHeight;
        fill(currentCell, jj, kk, cellVol);
        floorArea -= newArea;
        
//        assert(cellVol <= voxelVol*(1+SMALL_NUM));
//        assert(cellVol >= voxelVol*(-SMALL_NUM));
    }
    
    double vol = floorArea*cellHeight;
//    double voxelVol = voxelSize()[0]*voxelSize()[1]*voxelSize()[2];
//    const double SMALL_NUM = 1e-8;
//    assert(vol <= voxelVol*(1 + SMALL_NUM));
//    assert(vol >= voxelVol*(-SMALL_NUM));
    fill(lastCell+1, voxelBounds().p2[rowAxis()], jj, kk, vol);
    
//    checkFillFactorsBounded(mFillFactors, 
//        voxelSize()[0]*voxelSize()[1]*voxelSize()[2]);
}


#pragma mark *** OrientationRowVoxelizer ***

OrientationRowVoxelizer::
OrientationRowVoxelizer(int rowAxis, Rect3d bounds, Rect3i voxelBounds, Vector3b nonSymmetricDimensions,
    RLE::DynamicRLE3<Matrix3d> & outOrientations) :
    RowVoxelizer(rowAxis, bounds, voxelBounds, nonSymmetricDimensions),
    mOrientations(outOrientations)
{
}

void OrientationRowVoxelizer::
voxelize(int jj, int kk, const vector<RayIntersection> & intersections,
    int numHits)
{
    // 1.  Process all the triangles and make a hit list
    
    vector<ClippedTriangle> clipped(numHits);
    vector<OrientationRecord> orientationRecords;
    
    for (int ii = 0; ii < numHits; ii++)
    {
        //cerr << "Hit " << intersections[ii].triangle()->triangle() << "\n";
        clipped[ii] = clipToRow(*intersections[ii].triangle(), jj, kk);
        //cerr << "\t clip to pb(" << jj << ", " << kk << ") = "
//            << pixelBounds(jj, kk) << "\n";
//        cerr << "\t" << clipped[ii] << "\n";
        voxelizeOneTri(&clipped[ii], jj, kk, orientationRecords);
    }
    
    sort(orientationRecords.begin(), orientationRecords.end());
    handleSortedOrientationRecords(jj, kk, orientationRecords);
}


void OrientationRowVoxelizer::
voxelizeOneTri(const ClippedTriangle* clipped, int jj, int kk,
    std::vector<OrientationRecord> & orientationRecords)
{
    if (clipped->vertices().size() < 3)
        return;
    
    Rect3d clippedBounds = clipped->bounds();
    
    double xmin = clippedBounds.p1[rowAxis()];
    double xmax = clippedBounds.p2[rowAxis()];
    
    int minVoxel = voxel(xmin);
    int maxVoxel = voxel(xmax);
    
    for (int vv = minVoxel; vv <= maxVoxel; vv++)
    {
        orientationRecords.push_back(OrientationRecord(vv, clipped));
    }
}

void OrientationRowVoxelizer::
handleSortedOrientationRecords(int jj, int kk,
    const vector<OrientationRecord> & records)
{
    if (records.empty())
        return;
    
    int cell;
    int currentRecord = 0;
    
    while (currentRecord < records.size())
    {
        cell = records[currentRecord].voxel;
        
        Matrix3d totalOrientation = Matrix3d::zero();
        
        while( currentRecord < records.size() &&
            records.at(currentRecord).voxel == cell )
        {
            ClippedTriangle polyInVoxel = records[currentRecord].polygon->
                partAbove(rowAxis(), voxelMin(cell));
            polyInVoxel = polyInVoxel.partBelow(rowAxis(), voxelMin(cell+1));
            
            if (polyInVoxel.vertices().size() > 0) // just to speed it up...
            {
                Vector3d normal = polyInVoxel.arealNormal();
                totalOrientation += outerProduct(normal, unit(normal));
            }
            
            currentRecord++;
        }
        
        if (trace(totalOrientation) > 0) // this part culls empty matrices.
            fill(cell, jj, kk, totalOrientation);
    }
}

#pragma mark *** DFillFactorRowVoxelizer ***

DFillFactorRowVoxelizer::
DFillFactorRowVoxelizer(int rowAxis, Rect3d bounds, Rect3i voxelBounds, Vector3b nonSymmetricDimensions,
    std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > & outDFillFactors) :
    RowVoxelizer(rowAxis, bounds, voxelBounds, nonSymmetricDimensions),
    mDFillFactors(outDFillFactors)
{
}

void DFillFactorRowVoxelizer::
voxelize(int jj, int kk, const vector<RayIntersection> & intersections,
    int numHits)
{
    // 1.  Process all the triangles and make a hit list
    
    vector<ClippedTriangle> clipped(numHits);
    vector<DFillRecord> records;
    
    for (int ii = 0; ii < numHits; ii++)
    {
        clipped[ii] = clipToRow(*intersections[ii].triangle(), jj, kk);
        voxelizeOneTri(&clipped[ii], jj, kk, records);
    }
    
    sort(records.begin(), records.end());
    handleSortedRecords(jj, kk, records);
}

void DFillFactorRowVoxelizer::
voxelizeOneTri(const ClippedTriangle* clipped, int jj, int kk,
    std::vector<DFillRecord> & records)
{
    if (clipped->vertices().size() < 3)
        return;
    
    Rect3d clippedBounds = clipped->bounds();
    
    double xmin = clippedBounds.p1[rowAxis()];
    double xmax = clippedBounds.p2[rowAxis()];
    
    int minVoxel = voxel(xmin);
    int maxVoxel = voxel(xmax);
    
    for (int vv = minVoxel; vv <= maxVoxel; vv++)
    {   
        records.push_back(DFillRecord(vv, clipped));
    }
}

static inline void sAdd(map<unsigned int, Vector3d> & sensitivities,
    const map<unsigned int, Vector3d> & addend)
{
    map<unsigned int, Vector3d>::const_iterator itr;
    for (itr = addend.begin(); itr != addend.end(); itr++)
    {
        if (0 == sensitivities.count(itr->first))
        {
            sensitivities[itr->first] = itr->second;
        }
        else
            sensitivities[itr->first] += itr->second;
    }
}

void DFillFactorRowVoxelizer::
handleSortedRecords(int jj, int kk, const vector<DFillRecord> & records)
{
    if (records.empty())
        return;
    
    int cell;
    int currentRecord = 0;
    
//    if (rowAxis() != 2)
//        throw(std::logic_error("ClippedTriangle::DV() assumes row axis is Z."));
    
    while (currentRecord < records.size())
    {
        cell = records[currentRecord].voxel;
        
        map<unsigned int, Vector3d> sensitivities;
        
        while( currentRecord < records.size() &&
            records.at(currentRecord).voxel == cell )
        {
            const ClippedTriangle* clipped = records[currentRecord].polygon;
            
//            cerr << "Handle " << *clipped << " from " << voxelMin(cell)
//                << " to " << voxelMin(cell+1) << "\n";
            
            ClippedTriangle polyInVoxel;
            SensitiveEdge edgeAbove, edgeBelow;
            
            if (fabs(clipped->triangle()->unitNormal()[rowAxis()]) < 0.999999)
            {
                Plane3d topPlane(Vector3d::unit(rowAxis()), -voxelMin(cell+1)),
                    bottomPlane(Vector3d::unit(rowAxis()), -voxelMin(cell));
                
                edgeAbove = SensitiveEdge(*clipped->triangle(), topPlane);
                edgeBelow = SensitiveEdge(*clipped->triangle(), bottomPlane);
                
                //
//                cerr << "Above: " << topPlane << " edge " << edgeAbove.line()
//                    << "\n";
//                cerr << "Below: " << bottomPlane << " edge " << edgeBelow.line() << "\n";
//                
                polyInVoxel = clipped->partAbove(rowAxis(), voxelMin(cell),
                    &edgeBelow);
                polyInVoxel = polyInVoxel.partBelow(rowAxis(), voxelMin(cell+1),
                    &edgeAbove);
            }
            else
                polyInVoxel = *clipped;
            
            if (polyInVoxel.vertices().size() > 2)
            {
                sAdd(sensitivities, polyInVoxel.DV());
            }
            
            currentRecord++;
        }
        
        fill(cell, jj, kk, sensitivities);
    }
}


#pragma mark *** DOrientationRowVoxelizer ***

Matrix3<Vector3d> DOrientationRowVoxelizer::mZedMatrix = Matrix3<Vector3d>();

DOrientationRowVoxelizer::
DOrientationRowVoxelizer(int rowAxis, Rect3d bounds, Rect3i voxelBounds, Vector3b nonSymmetricDimensions,
    std::map<unsigned int, RLE::DynamicRLE3<Matrix3<Vector3d> > > & outDOrientation) :
    RowVoxelizer(rowAxis, bounds, voxelBounds, nonSymmetricDimensions),
    mDOrientation(outDOrientation)
{
}

// This version had the issue with insufficient sensitivity information.
//// this is THE SAME as the fill factor voxelizer.
//void DOrientationRowVoxelizer::
//voxelize(int jj, int kk, const vector<RayIntersection> & intersections,
//    int numHits)
//{
//    // 1.  Process all the triangles and make a hit list
//    
//    vector<ClippedTriangle> clipped(numHits);
//    vector<DOrientationRecord> records;
//    int numFills = 0;
//    
//    for (int ii = 0; ii < numHits; ii++)
//    {
//        clipped[ii] = clipToRow(*intersections[ii].triangle(), jj, kk);
//        
////        cerr << "Clipped " << *intersections[ii].triangle() << " to "
////            << clipped[ii] << "\n";
//        
//        voxelizeOneTri(&clipped[ii], jj, kk, records);
//    }
//    
//    sort(records.begin(), records.end());
//    handleSortedRecords(jj, kk, records);
//}

// Old way.
// This is the same as the fill factor function.
//void DOrientationRowVoxelizer::
//voxelizeOneTri(const ClippedTriangle* clipped, int jj, int kk,
//    std::vector<DOrientationRecord> & records)
//{
//    if (clipped->vertices().size() < 3)
//        return;
//    
//    Rect3d clippedBounds = clipped->bounds();
//    
//    double xmin = clippedBounds.p1[rowAxis()];
//    double xmax = clippedBounds.p2[rowAxis()];
//    
//    int minVoxel = voxel(xmin);
//    int maxVoxel = voxel(xmax);
//    
//    int xyz = dominantDirection(clipped->triangle()->unitNormal());
//    
//    for (int vv = minVoxel; vv <= maxVoxel; vv++)
//    {   
//        records.push_back(DOrientationRecord(vv, clipped, xyz));
//    }
//}

// new way.
void DOrientationRowVoxelizer::
voxelize(int jj, int kk, const vector<RayIntersection> & intersections,
    int numHits)
{
    // 1.  Process all the triangles and make a hit list
    
    vector<SensitiveEdge> extraEdges(6);
    
    for (int ii = 0; ii < numHits; ii++)
    {
        ClippedTriangle clipped = clipToRow(*intersections[ii].triangle(),
            jj, kk, extraEdges);
        
//        cerr << "Clipped " << *intersections[ii].triangle() << " to "
//            << clipped << "\n";
        
        voxelizeOneTri(&clipped, jj, kk);
    }
}

static inline void sAdd(map<unsigned int, Matrix3<Vector3d> > & sensitivities,
    const map<unsigned int, Matrix3<Vector3d> > & addend)
{
    map<unsigned int, Matrix3<Vector3d> >::const_iterator itr;
    for (itr = addend.begin(); itr != addend.end(); itr++)
    {
        if (0 == sensitivities.count(itr->first))
        {
            sensitivities[itr->first] = itr->second;
        }
        else
            sensitivities[itr->first] += itr->second;
    }
}


// new way.
void DOrientationRowVoxelizer::
voxelizeOneTri(const ClippedTriangle* clipped, int jj, int kk)
{
    if (clipped->vertices().size() < 3)
        return;
    
    Rect3d clippedBounds = clipped->bounds();
    
    double xmin = clippedBounds.p1[rowAxis()];
    double xmax = clippedBounds.p2[rowAxis()];
    
    int minVoxel = voxel(xmin);
    int maxVoxel = voxel(xmax);
    
    int upDirection = clipped->triangle()->upAxis();
    
//    cerr << "===== voxelizeOneTri()\n";
    
    for (int cell = minVoxel; cell <= maxVoxel; cell++)
    {
        map<unsigned int, Matrix3<Vector3d> > sensitivities;
        
//            cerr << "\n=========\n";
//            cerr << "Handle " << *clipped << " cell " << cell << "\n";
        
        ClippedTriangle polyInVoxel = *clipped;
        SensitiveEdge edgeAbove, edgeBelow;
        
        if (fabs(clipped->triangle()->unitNormal()[rowAxis()]) < 0.999999)
        {
            Plane3d topPlane(Vector3d::unit(rowAxis()), -voxelMin(cell+1)),
                bottomPlane(Vector3d::unit(rowAxis()), -voxelMin(cell));
            
//            cerr << "topPlane " << topPlane << "\n";
//            cerr << "bottomPlane " << bottomPlane << "\n";
            
            if (xmax >= voxelMin(cell+1) && xmin <= voxelMin(cell+1))
            {
                edgeAbove = SensitiveEdge(*clipped->triangle(), topPlane);
                polyInVoxel = polyInVoxel.partBelow(rowAxis(), voxelMin(cell+1),
                    &edgeAbove);
            }
            else
                polyInVoxel = polyInVoxel.partBelow(rowAxis(), voxelMin(cell+1));
            
            if (xmax >= voxelMin(cell) && xmin <= voxelMin(cell))
            {
                edgeBelow = SensitiveEdge(*clipped->triangle(), bottomPlane);
                polyInVoxel = polyInVoxel.partAbove(rowAxis(), voxelMin(cell),
                    &edgeBelow);
            }
            else
                polyInVoxel = polyInVoxel.partAbove(rowAxis(), voxelMin(cell));
//            cerr << "Above: " << topPlane << " edge " << edgeAbove.line()
//                << "\n";
//            cerr << "Below: " << bottomPlane << " edge " << edgeBelow.line() << "\n";
        }
        
        if (polyInVoxel.vertices().size() > 2)
        {
            sAdd(sensitivities, polyInVoxel.DNNT(upDirection));
        }
        
        fill(cell, jj, kk, sensitivities);
    }
//    cerr << "===== end\n";
}


//void DOrientationRowVoxelizer::
//handleSortedRecords(int jj, int kk, const vector<DOrientationRecord> & records)
//{
//    if (records.empty())
//        return;
//    
//    int cell;
//    int currentRecord = 0;
//    
//    cerr << "===== handleSortedRecords()\n";
//    while (currentRecord < records.size())
//    {
//        cell = records[currentRecord].voxel;
//        
//        map<unsigned int, Matrix3<Vector3d> > sensitivities;
//        
//        while( currentRecord < records.size() &&
//            records.at(currentRecord).voxel == cell )
//        {
//            const ClippedTriangle* clipped = records[currentRecord].polygon;
//            int upDirection = records[currentRecord].axis;
//            
////            cerr << "\n=========\n";
////            cerr << "Handle " << *clipped << " cell " << cell << "\n";
//            
//            ClippedTriangle polyInVoxel;
//            SensitiveEdge edgeAbove, edgeBelow;
//            
//            if (fabs(clipped->triangle()->unitNormal()[rowAxis()]) < 0.999999)
//            {
//            
//                Plane3d topPlane(Vector3d::unit(rowAxis()), -voxelMin(cell+1)),
//                    bottomPlane(Vector3d::unit(rowAxis()), -voxelMin(cell));
//                
//                cerr << "topPlane " << topPlane << "\n";
//                cerr << "bottomPlane " << bottomPlane << "\n";
//                
//                edgeAbove = SensitiveEdge(
//                    Line(clipped->triangle()->triangle(), topPlane),
//                    clipped->triangle()->controlVertices());
//                edgeBelow = SensitiveEdge(
//                    Line(clipped->triangle()->triangle(), bottomPlane),
//                    clipped->triangle()->controlVertices());
//                
//                cerr << "Above: " << topPlane << " edge " << edgeAbove.line()
//                    << "\n";
//                cerr << "Below: " << bottomPlane << " edge " << edgeBelow.line() << "\n";
//                
//                polyInVoxel = clipped->partAbove(rowAxis(), voxelMin(cell),
//                    &edgeBelow);
//                polyInVoxel = polyInVoxel.partBelow(rowAxis(), voxelMin(cell+1),
//                    &edgeAbove);
//            }
//            else
//                polyInVoxel = *clipped;
//            
//            if (polyInVoxel.vertices().size() > 2)
//            {
//                sAdd(sensitivities, polyInVoxel.DNNT(upDirection));
//                
////                if (jj == 0 && kk == 0 && cell == 0)
////                {
////                    cerr << "Poly " << polyInVoxel << " area "
////                        << polyInVoxel.arealNormal() << " -> "
////                        << polyInVoxel.area2d(upDirection) << ":\n";
////                    
////                    map<unsigned int, Matrix3<Vector3d> > dnnt =
////                        polyInVoxel.DNNT(upDirection);
////                    for (int vv = 0; vv < 16; vv++)
////                    if (dnnt.count(vv))
////                    {
////                        Matrix3d dOri;
////                        for (int mm = 0; mm < 9; mm++)
////                            dOri[mm] = dnnt[vv][mm][0];
////                        if (vectorNorm(dOri) > 1e-4)
////                            cerr << "\tv" << vv << ":\n" << dOri << "\n";
////                    }
////                }
//                
//            }
//            
//            currentRecord++;
//        }
//        
//        fill(cell, jj, kk, sensitivities);
//    }
//    cerr << "===== end\n";
//}




