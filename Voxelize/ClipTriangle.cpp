/*
 *  ClipTriangle.cpp
 *  Voxelize
 *
 *  Created by Paul C Hansen on 6/28/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 */

#include "ClipTriangle.h"

#include <limits>

using namespace std;

ClippedTriangle::
ClippedTriangle() :
    mTriangle(0L),
    mVertices(),
    mEdges()
{
}


ClippedTriangle::
ClippedTriangle(const SimpleMesh::Triangle & tri) :
    mTriangle(&tri),
    mVertices(3),
    mEdges(3)
{
    for (int edge = 0; edge < 3; edge++)
    {
        mVertices[edge] = tri.triangle()[edge];
        mEdges[edge] = &tri.edgeSensitivity(edge);
    }
    //checkForInconsistentEdges();
}

ClippedTriangle ClippedTriangle::
partBelow(int axis, double val,
    const SimpleMesh::SensitiveEdge* clippedSensitivity) const
{
    ClippedTriangle clipped;
    clipped.mTriangle = mTriangle;
    
// This assertion is erroneous now.  I have to make sure I don't USE any of
// these bogus lines later on; that's where I carry the risk now.
//    if (clippedSensitivity != 0L)
//        assert(sumSquares(clippedSensitivity->line().v()) > 0);
    
    for (int nn = 0; nn < mVertices.size(); nn++)
    {
        const Vector3d & pCurrent = mVertices[nn];
        const Vector3d & pNext = mVertices[(nn+1)%mVertices.size()];
        
        if (pCurrent[axis] < val) // currently inside (below)
        {
            clipped.mVertices.push_back(pCurrent);
            clipped.mEdges.push_back(mEdges[nn]);
            
            if (pNext[axis] > val) // inside-to-totally-outside
            {
                Vector3d intersection = pCurrent + (pNext - pCurrent)*
                    (val - pCurrent[axis])/(pNext[axis] - pCurrent[axis]);
                clipped.mVertices.push_back(intersection);
                clipped.mEdges.push_back(clippedSensitivity);
            }
        }
        else if (pCurrent[axis] == val) // currently outside, on boundary
        {
            clipped.mVertices.push_back(pCurrent);
            if (pNext[axis] < val) // boundary-to-inside
                clipped.mEdges.push_back(mEdges[nn]);
            else
            {
                clipped.mEdges.push_back(clippedSensitivity); // edge outside
            }
        }
        else if (pNext[axis] < val) // outside-to-inside, crossing boundary
        {
            Vector3d intersection = pCurrent + (pNext - pCurrent)*
                (val - pCurrent[axis])/(pNext[axis] - pCurrent[axis]);
            
            clipped.mVertices.push_back(intersection);
            clipped.mEdges.push_back(mEdges[nn]);
        }
        else
            ; // continue
            
        //clipped.checkForInconsistentEdges();
    }
    
    if (clipped.vertices().size() < 2)
    {
        clipped.mVertices.clear();
        clipped.mEdges.clear();
    }
    
    //clipped.checkForInconsistentEdges();
    return clipped;
}

ClippedTriangle ClippedTriangle::
partAbove(int axis, double val,
    const SimpleMesh::SensitiveEdge* clippedSensitivity) const
{
    ClippedTriangle clipped;
    clipped.mTriangle = mTriangle;
    
// This assertion is erroneous now.  I have to make sure I don't USE any of
// these bogus lines later on; that's where I carry the risk now.
//    if (clippedSensitivity != 0L)
//        assert(sumSquares(clippedSensitivity->line().v()) > 0);
    
    for (int nn = 0; nn < mVertices.size(); nn++)
    {
        const Vector3d & pCurrent = mVertices[nn];
        const Vector3d & pNext = mVertices[(nn+1)%mVertices.size()];
        
//        if (normInf(pCurrent - pNext) < 1e-6)
//            int dbg = 1;
        
        if (pCurrent[axis] > val) // currently inside, in interior
        {
            clipped.mVertices.push_back(pCurrent);
            clipped.mEdges.push_back(mEdges[nn]);
            
            if (pNext[axis] < val) // inside-to-outside
            {
                Vector3d intersection = pCurrent + (pNext - pCurrent)*
                    (val - pCurrent[axis])/(pNext[axis] - pCurrent[axis]);
                clipped.mVertices.push_back(intersection);
                clipped.mEdges.push_back(clippedSensitivity);
            }
        }
        else if (pCurrent[axis] == val) // inside, on boundary
        {
            clipped.mVertices.push_back(pCurrent);
            
            if (pNext[axis] >= val) // inside-inside on boundary
                clipped.mEdges.push_back(mEdges[nn]);
            else // inside-outside from boundary
            {
                clipped.mEdges.push_back(clippedSensitivity);
            }
        }
        else if (pNext[axis] > val) // outside-to-inside, crossing boundary
        {
            Vector3d intersection = pCurrent + (pNext - pCurrent)*
                (val - pCurrent[axis])/(pNext[axis] - pCurrent[axis]);
            
            clipped.mVertices.push_back(intersection);
            clipped.mEdges.push_back(mEdges[nn]);
        }
        else
            ; // continue
    }
    
    if (clipped.vertices().size() < 2)
    {
        clipped.mVertices.clear();
        clipped.mEdges.clear();
    }
    
//    clipped.checkForInconsistentEdges();
    return clipped;
}

ClippedTriangle ClippedTriangle::partInside(const Rect2d & xyBounds) const
{
    ClippedTriangle clipt;
    
    const int X = 0, Y = 1;
    
    clipt = partBelow(X, xyBounds.p2[0]);
    clipt = clipt.partBelow(Y, xyBounds.p2[1]);
    clipt = clipt.partAbove(X, xyBounds.p1[0]);
    clipt = clipt.partAbove(Y, xyBounds.p1[1]);
    
    return clipt;
}

Rect3d ClippedTriangle::
bounds() const
{
    Rect3d bbox(numeric_limits<double>::max(),
        numeric_limits<double>::max(),
        numeric_limits<double>::max(),
        -numeric_limits<double>::max(),
        -numeric_limits<double>::max(),
        -numeric_limits<double>::max());
    
    for (int vv = 0; vv < mVertices.size(); vv++)
    {
        bbox.p1 = vec_min(bbox.p1, mVertices[vv]);
        bbox.p2 = vec_max(bbox.p2, mVertices[vv]);
    }
    
    return bbox;
}

static Vector2d
s2d(const Vector3d & v, int projectionAxis)
{
    return Vector2d(v[(projectionAxis+1)%3], v[(projectionAxis+2)%3]);
}

double ClippedTriangle::
area2d(int axis) const
{
    double totalArea = 0.0;
    
    for (int vert = 0; vert+2 < vertices().size(); vert++)
    {
        totalArea += 0.5*cross(
            s2d(vertices()[vert+1] - vertices()[0], axis),
            s2d(vertices()[vert+2] - vertices()[0], axis) );
    }
    
    return totalArea;
}

double ClippedTriangle::
volume(int axis, double refHeight) const
{
    double totalVolume = 0.0;
    
    for (int vert = 0; vert+2 < vertices().size(); vert++)
    {
        double area = 0.5*cross(
            s2d(vertices()[vert+1] - vertices()[0], axis),
            s2d(vertices()[vert+2] - vertices()[0], axis) );
        
        Vector3d centroid = 0.3333333333333333333333333333333333333333333333 * (
            vertices()[0] + vertices()[vert+1] + vertices()[vert+2] );
        
        totalVolume += area*(centroid[axis] - refHeight);
    }
    
    return totalVolume;
}

Vector3d ClippedTriangle::
arealNormal() const
{
    Vector3d totalNormal(0.0, 0.0, 0.0);
    
    for (int vert = 0; vert+2 < vertices().size(); vert++)
    {
        totalNormal += 0.5*cross(
            vertices()[vert+1] - vertices()[0],
            vertices()[vert+2] - vertices()[0] );
    }
    
    return totalNormal;
}

Matrix3d ClippedTriangle::
nnTA() const
{
    return outerProduct(arealNormal(), unit(arealNormal()));
}

map<unsigned int, Vector3d> ClippedTriangle::
DV() const
{
    map<unsigned int, Vector3d> outJacobian;
    
    const Vector3d & v0 = triangle()->controlVertices()[0].point();
    const Vector3d & v1 = triangle()->controlVertices()[1].point();
    const Vector3d & v2 = triangle()->controlVertices()[2].point();
    int cv0 = triangle()->controlVertices()[0].id();
    int cv1 = triangle()->controlVertices()[1].id();
    int cv2 = triangle()->controlVertices()[2].id();
    
    Vector3d normal = triangle()->unitNormal();
    
    Matrix3d M = inverse(Matrix3d::withColumns(
        v1-v0, v2-v0, normal));
    Vector3d Mrow2 = M.row(2);
    
    for (int vert = 0; vert+2 < vertices().size(); vert++)
    {
        double area = 0.5 * dot(normal, cross(
            vertices()[vert+1] - vertices()[0],
            vertices()[vert+2] - vertices()[0]));
        Vector3d centroid = (1.0/3)*
            (vertices()[0] + vertices()[vert+1] + vertices()[vert+2] );
        
        Vector3d abc = M*(centroid - v0);
        
        outJacobian[cv0] -= area*(abc[0] + abc[1] - 1)*Mrow2;
        outJacobian[cv1] += area*abc[0]*Mrow2;
        outJacobian[cv2] += area*abc[1]*Mrow2;
    }
    
    return outJacobian;
}

map<unsigned int, Matrix3<Vector3d> > ClippedTriangle::
DNNT(int axis) const
{
    map<unsigned int, Matrix3<Vector3d> > outJacobian;//, jacArea, jacBdy;
    
//    cerr << *this << ": \n";
    intNNTDh(outJacobian);
    intDnnT(outJacobian);
    
//    intNNTDh(jacBdy);
//    intDnnT(jacArea);
//    
//    cerr << "Area integral:\n";
//    printJacobian(jacArea);
//    cerr << "Boundary integral:\n";
//    printJacobian(jacBdy);
//    //printJacobian(outJacobian);
    
    return outJacobian;
}

// This algorithm is written out in my notes for 8/10/11.
void ClippedTriangle::
intDz(map<unsigned int, Vector3d> & outIntegrals) const
{
    const int Z = 2;
    const double THRESHOLD = 1e-9;
    if (fabs(arealNormal()[Z]) < THRESHOLD)
        return;
    
    const Vector3d & v0 = triangle()->controlVertices()[0].point();
    const Vector3d & v1 = triangle()->controlVertices()[1].point();
    const Vector3d & v2 = triangle()->controlVertices()[2].point();
    int cv0 = triangle()->controlVertices()[0].id();
    int cv1 = triangle()->controlVertices()[1].id();
    int cv2 = triangle()->controlVertices()[2].id();
    
    Matrix3d M = inverse(Matrix3d::withColumns(
        v1-v0, v2-v0, -Vector3d::unit(Z)));
    Vector3d MrowZ = M.row(Z); // turns out we need this!  :-D
    
    for (int vert = 0; vert+2 < vertices().size(); vert++)
    {
        double area = 0.5*cross(
            s2d(vertices()[vert+1] - vertices()[0], Z),
            s2d(vertices()[vert+2] - vertices()[0], Z) );
        Vector3d centroid = 0.3333333333333333333333333333333333333333333333 * (
            vertices()[0] + vertices()[vert+1] + vertices()[vert+2] );
        
        Vector3d b = Vector3d(centroid[0], centroid[1], 0) - v0;
        Vector3d abc = M*b;
        
//        if (isnan(area) || isnan(abc[0]) || isnan(abc[1]) || isnan(abc[2]))
//            int ohNo = 'YES';
        
//        if (isnan(normInf(M)))
//            int uhOh = 'TRUE';
        outIntegrals[cv0] += area*(abc[0]+abc[1]-1)*MrowZ;
        outIntegrals[cv1] -= area*abc[0]*MrowZ;
        outIntegrals[cv2] -= area*abc[1]*MrowZ;
        
//        outIntegrals[cv0][Z] += area*(1 - abc[0] - abc[1]);
//        outIntegrals[cv1][Z] += area*abc[0];
//        outIntegrals[cv2][Z] += area*abc[1];
    }
}

// This algorithm is written out in my notes for 8/10/11
void ClippedTriangle::
intZDh(map<unsigned int, Vector3d> & outIntegrals, int upAxis,
    double refHeight) const
{
//    const int X = (upAxis+1)%3, Y = (upAxis+2)%3,
    const int Z = upAxis;
    //const int Z = 0, X = 1, Y = 2;
//    const int X = 0, Y = 1, Z = 2;
    for (int ee = 0; ee < edges().size(); ee++)
    if (edges()[ee] != 0L)
    {
        const SimpleMesh::SensitiveEdge* edge = edges()[ee];
        Vector3d q0 = vertices()[ee];
        Vector3d q1 = vertices()[(ee+1)%vertices().size()];
        double z0 = q0[Z] - refHeight; // generally all we want is the height
        double z1 = q1[Z] - refHeight;
        
        if (normInf(q0 - q1) < 1e-10)
            continue;
        
        // Project most of the geometry to the XY plane.
        Vector2d p0 = s2d(q0, Z);
        Vector2d p1 = s2d(q1, Z);
        
        // If the line goes straight up there will be no contribution AND the
        // linear algebra will barf.
        if (sumSquares(p1-p0) < 1e-9)
            continue;
        
        Vector2d x = s2d(edge->line().x(), Z);
        Vector2d v = s2d(edge->line().v(), Z);
        Vector2d l(-v[1], v[0]); // a direction perpendicular to the edge
        
        // Use some linear algebra to express p0 and p1 in terms of a distance
        // along the edge.
        Matrix2d M = inverse(Matrix2d::withColumns(-v, l));
        Vector2d ts0 = M*(x - p0);
        Vector2d ts1 = M*(x - p1);
        double & t0 = ts0[0]; // we only care about the first number anyway.
        double & t1 = ts1[0];
        
        assert(edge->controlVertices().size() ==
            edge->line().numControlVertices());
        
        for (int vv = 0; vv < edge->controlVertices().size(); vv++)
        for (int xyz = 0; xyz < 3; xyz++)
        {
            Vector2d Dx = s2d(edge->line().Dx(vv).column(xyz), Z);
            Vector2d Dv = s2d(edge->line().Dv(vv).column(xyz), Z);
            
            Vector2d Dts0 = M*(Dx + t0*Dv);
            Vector2d Dts1 = M*(Dx + t1*Dv);
            double & Ds0 = Dts0[1]; // now we only care about the second number
            double & Ds1 = Dts1[1];
            
            Vector2d Dp0 = l*Ds0;
            Vector2d Dp1 = l*Ds1;
            
            double integral = -1.0/6.0 * (
                (2*z0 + z1)*cross(p1-p0, Dp0) + (z0 + 2*z1)*cross(p1-p0, Dp1) );
            
//            if (isnan(integral))
//                int noNoNoNo = 4*'NO';
            
            int controlVertexId = edge->controlVertices()[vv].id();
            outIntegrals[controlVertexId][xyz] += integral;
        }
    }
}


void ClippedTriangle::
intNNTDh(std::map<unsigned int, Matrix3<Vector3d> > & outIntegrals)
    const
{
    Matrix3d nnTdA = triangle()->nnTdA(); // nnT/nz
    const int Z = triangle()->upAxis();
    double zSign = triangle()->unitNormal()[Z] > 0 ? 1.0 : -1.0;
    
//    cerr << "\nintNNTDh()\n";
//    cerr << "Tri " << *this << "\n";
//    cerr << "Up is " << char('x'+Z) << "\n";
//    cerr << "nnTdA = \n" << nnTdA << "\n";
    
    for (int ee = 0; ee < edges().size(); ee++)
    if (edges()[ee] != 0L)
//    if (edges()[ee]->line().numControlVertices() != 0)
//    if (edges()[ee]->line().v() != Vector3d(0,0,0))
    {
        const SimpleMesh::SensitiveEdge* edge = edges()[ee];
        Vector3d q0 = vertices()[ee];
        Vector3d q1 = vertices()[(ee+1)%vertices().size()];
        
//        cerr << "quickLine(" << q0 << ", " << q1 << ", 'ro-');\n";
//        cerr << "\t" << edges()[ee]->line().numControlVertices() << " CVs\n"
//            << "\t velocity " << edges()[ee]->line().v() << "\n";
            
        if (normInf(q0 - q1) < 1e-10)
            continue;
        
        if (edges()[ee]->line().numControlVertices() == 0)
            continue;
        
        // This should never occur when q0 and q1 are at some distance from
        // each other, I'm guessing.
        assert(sumSquares(edge->line().v()) > 0);
        
//        cerr << "edge-line dist = "
//            << pointLineDistance(q0, q1, edge->line().x()) << "\n";
        
//        assert(pointLineDistance(q0, q1, edge->line().x()) < 1e-6);
        
        
        // Project most of the geometry to the XY plane.
        Vector2d p0 = s2d(q0, Z);
        Vector2d p1 = s2d(q1, Z);
        
        // If the line goes straight up there will be no contribution AND the
        // linear algebra will barf.
        if (sumSquares(p1-p0) < 1e-9)
            continue;
        
        Vector2d x = s2d(edge->line().x(), Z);
        Vector2d v = s2d(edge->line().v(), Z);
        Vector2d l(-v[1], v[0]); // a direction perpendicular to the edge
        
        // Use some linear algebra to express p0 and p1 in terms of a distance
        // along the edge.
        Matrix2d M = inverse(Matrix2d::withColumns(-v, l));
        Vector2d ts0 = M*(x - p0);
        Vector2d ts1 = M*(x - p1);
        double & t0 = ts0[0]; // we only care about the first number anyway.
        double & t1 = ts1[0];
        
        assert(edge->controlVertices().size() ==
            edge->line().numControlVertices());
        
//        cerr << "\n\nEdge #" << ee << " " << q0 << " to " << q1 << " along "
//            << x << "+t" << v << " perp " << l << ", ts = " << ts0 << ", "
//            << ts1 << "\n";
//        cerr << "Line " << edge->line().x() << " along " << edge->line().v()
//            << "\n";
//        cerr << "M: \n" << M << "\n";
//        
//        cerr << "plot3([" << q0[0] << " " << q1[0] << "], ["
//            << q0[1] << " " << q1[1] << "], ["
//            << q0[2] << " " << q1[2] << "], 'ro-');\n";
        
        for (int vv = 0; vv < edge->controlVertices().size(); vv++)
        {
            Matrix3<Vector3d> saved;
            int controlVertexId = edge->controlVertices()[vv].id();
            for (int xyz = 0; xyz < 3; xyz++)
            {
    //            cerr << "\t" << vv << "." << xyz;
                
                Vector2d Dx = s2d(edge->line().Dx(vv).column(xyz), Z);
                Vector2d Dv = s2d(edge->line().Dv(vv).column(xyz), Z);
                
    //            cerr << " Dx " << Dx << " Dv " << Dv << " ";
                
                Vector2d Dts0 = M*(Dx + t0*Dv);
                Vector2d Dts1 = M*(Dx + t1*Dv);
                double & Ds0 = Dts0[1]; // now we only care about the second number
                double & Ds1 = Dts1[1];
                
                Vector2d Dp0 = l*Ds0;
                Vector2d Dp1 = l*Ds1;
                
    //            Vector2d Dp0_b = Dx + t0*Dv + v*Dts0[0];
                
                double areaTerm = 0.5*cross(p1-p0, Dp0+Dp1)*zSign;
    //            cerr << " Dp0 " << Dp0 << " Dp1 " << Dp1 << " area term "
    //                << areaTerm << "\n";
                
                for (int ii = 0; ii < 3; ii++)
                for (int jj = 0; jj < 3; jj++)
                {
                    outIntegrals[controlVertexId](ii,jj)[xyz] -=
                        nnTdA(ii,jj)*areaTerm;
                    saved(ii,jj)[xyz] -= nnTdA(ii,jj)*areaTerm;
                }
            }
//            cerr << "\t" << controlVertexId << ": \n";
//            cerr << saved << "\n";
        }
    }
//    cerr << "end intNNTDh()\n";
}

static Matrix3<Vector3d> sTimes(double d, const Matrix3<Vector3d> & m)
{
    return Matrix3<Vector3d>( d*m[0], d*m[1], d*m[2],
        d*m[3], d*m[4], d*m[5],
        d*m[6], d*m[7], d*m[8] );
}

void ClippedTriangle::
intDnnT(std::map<unsigned int, Matrix3<Vector3d> > & outIntegrals)
    const
{
    double area = fabs(area2d(triangle()->upAxis()));
//    cerr << "Control tri:\n"
//        << "quickPatch(" << triangle()->controlTriangle() << ", 'y');\n";
//    cerr << "Self: quickPatch("
//        << *this << ", 'g');\n";
//    cerr << "AREA " << area << "\n";
    
    for (int vv = 0; vv < 3; vv++)
    {
        int controlVertexId = triangle()->controlVertices()[vv].id();
        outIntegrals[controlVertexId] += sTimes(area, triangle()->Dorientation(vv));
        
        //cerr << "\tDNNT for CV " << controlVertexId << ":\n";
//        cerr << "dnnt.v" << controlVertexId << " =\n";
//        cerr << triangle()->Dorientation(vv) << ";\n";
    }
}

SimpleMesh::SensitiveEdge & ClippedTriangle::
storedEdge(int ee)
{
    assert(ee >= 0 && ee <= 5);
    
    return mEdgeStorage[ee];
}


void ClippedTriangle::
checkForInconsistentEdges() const
{
    for (int ee = 0; ee < edges().size(); ee++)
    if (edges()[ee] != 0L)
    {
        const SimpleMesh::SensitiveEdge* edge = edges()[ee];
        Vector3d q0 = vertices()[ee];
        Vector3d q1 = vertices()[(ee+1)%vertices().size()];
        
        if (q1 != q0 && sumSquares(edge->line().v()) == 0)
        {
            cerr << "Got a wacky one.\n";
            cerr << *this << "\n";
        }
    }
}

ostream & operator<<(ostream & str, const ClippedTriangle & ct)
{
    str << "[";
    
    for (int vv = 0; vv+1 < ct.vertices().size(); vv++)
    {
        str << ct.vertices()[vv] << "; ";
    }
    if (ct.vertices().size() > 0)
        str << ct.vertices()[ct.vertices().size()-1];
    str << "]";
    
    return str;
}

void printJacobian(const map<unsigned int, Vector3d> & dFdv)
{
    map<unsigned int, Vector3d>::const_iterator itr;
    for (itr = dFdv.begin(); itr != dFdv.end(); itr++)
    {
        cerr << "jac.v" << itr->first << " = " << itr->second << ";\n";
    }
}
void printJacobian(const map<unsigned int, Matrix3<Vector3d> > & dFdv)
{
    map<unsigned int, Matrix3<Vector3d> >::const_iterator itr;
    for (itr = dFdv.begin(); itr != dFdv.end(); itr++)
    {
        cerr << "jac.v" << itr->first << " = " << itr->second << ";\n";
    }
}
