/*
 *  SimpleMesh.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 9/8/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "SimpleMesh.h"
#include "Search.h"
#include <stdexcept>

using namespace std;

namespace SimpleMesh
{

//int ControlVertex::sNextId = 0;


std::ostream & operator<<(std::ostream & str, const ControlVertex & cv)
{
    str << cv.id() << ": " << cv.point() << " " << cv.freeDirections();
    return str;
}

Line::
Line() :
    m_x(),
    m_v(),
    mNumControlVertices(0)
{
}

Line::
Line(const Vector3d & vert0, const Vector3d & vert1) :
    m_x(vert0),
    m_v(vert1 - vert0),
    mNumControlVertices(2)
{
    mDx[0] = Matrix3d::eye();
    mDx[1] = Matrix3d::zero();
    
    mDv[0] = -Matrix3d::eye();
    mDv[1] = Matrix3d::eye();
    
    assert(sumSquares(m_v) != 0);
}

Line::
Line(const Triangle3d & tri0, const Triangle3d & tri1) :
    mNumControlVertices(6)
{
    Vector3d n0 = unit(tri0.normal());
    Vector3d n1 = unit(tri1.normal());
    const double THRESHOLD = 1e-8;
    
    if (fabs(dot(n0, n1)) > (1.0 - THRESHOLD))
        throw(std::logic_error("Input triangles may be coplanar"));
    
    // Pick two edges of tri0 to intersect with tri1.  They have to not be
    // coplanar with tri1.
    int edge0, edge1, ee;
    
    for (ee = 0; ee < 3; ee++)
    {
        Vector3d edge = tri0.edgeDirection(ee);
        if (fabs(dot(edge, n1)) > 1e-8*norm(edge))
            break;
    }
    edge0 = ee;
    
    for (ee = (edge0+1)%3; ee != edge0; ee = (ee+1)%3)
    {
        Vector3d edge = tri0.edgeDirection(ee);
        if (fabs(dot(edge, n1)) > 1e-8*norm(edge))
            break;
    }
    edge1 = ee;
    
    assert(edge0 != edge1);
    
    // --- Intersection of tri0 and tri1: find two points on the line.
    Vector3d p0, p1;
    p0 = intersection(tri0, tri1, edge0);
    p1 = intersection(tri0, tri1, edge1);
    
    m_x = p0;
    m_v = cross(n0, n1); 
    //m_v = p1 - p0;
    
    assert(sumSquares(m_v) != 0);
    
    Triangle3d zedTri(Vector3d(0,0,0), Vector3d(0,0,0), Vector3d(0,0,0));
    
    // Sensitivity of the point w.r.t. triangle 0 and triangle 1.
    for (int triVert = 0; triVert < 3; triVert++)
    {
        for (int direction = 0; direction < 3; direction++)
        {
            Triangle3d Dtri0(zedTri), Dtri1(zedTri);
            Dtri0[triVert][direction] = 1.0;
            Dtri1[triVert][direction] = 1.0;
            
            Vector3d dPt_dVert0 = dIntersection(tri0, tri1, Dtri0, zedTri, edge0);
            Vector3d dPt_dVert1 = dIntersection(tri0, tri1, zedTri, Dtri1, edge0);
            //col1 = dIntersection(tri0, tri1, Dtri0, zedTri, edge1);
            
            for (int row = 0; row < 3; row++)
            {
                mDx[triVert](row,direction) = dPt_dVert0[row];
                mDx[triVert+3](row,direction) = dPt_dVert1[row];
            }
        }
    }
    
    // Sensitivity of the line direction w.r.t. triangle 0 and triangle 1
    
    for (int triVert = 0; triVert < 3; triVert++)
    {
        Matrix3d dn0dv = tri0.dudv(triVert);
        Matrix3d dn1dv = tri1.dudv(triVert);
        
        for (int direction = 0; direction < 3; direction++)
        {
            // Change in tri[triVert] along direction:
            
            Vector3d dndt0 = cross(dn0dv.column(direction), n1);
            Vector3d dndt1 = cross(n0, dn1dv.column(direction));
            
            for (int row = 0; row < 3; row++)
            {
                mDv[triVert](row,direction) = dndt0[row];
                mDv[triVert+3](row,direction) = dndt1[row];
            }
        }
    }
    
}

Line::
Line(const Triangle3d & tri, const Plane3d & plane) :
    mNumControlVertices(3)
{
    const double THRESHOLD = 1e-10;
    
    Matrix3d dxdp0(Matrix3d::zero());
    Matrix3d dxdp1(Matrix3d::zero());
    int vNum0, vNum1, vNum2; // indices in [0 2] of tri vertices for our edge
    
    double lowestScore = 0;
    
    bool foundPoint = false;
    for (int edge = 0; edge < 3; edge++)
    {
        if (fabs(dot(tri.edgeDirection(edge), plane.normal())) < THRESHOLD)
        {
            //cerr << "tri " << tri.edgeDirection(edge) << " dot "
            //    << dot(tri.edgeDirection(edge), plane.normal()) << "\n";
            continue;
        }
        
        // Define a line along the triangle's edge: x(t) = t*p0 + (1-t)*p1.
        // Find t for intersection with the plane.  Take some derivatives.
        
        const Vector3d & p0 = tri[edge];
        const Vector3d & p1 = tri[(edge+1)%3];
        
        double t = (dot(plane.normal(), p1) + plane.constant()) /
            dot(plane.normal(), p1 - p0);
        
        // Pick an edge if the edge-plane intersection is substantially
        // better than the previous intersection.
        
        double currentScore = fabs(t-0.5);
        
        if (foundPoint && currentScore > (lowestScore - 0.1))
            continue;
        lowestScore = currentScore;
        
        m_x = p0*t + p1*(1.0-t);
        
//        cerr << "p0 " << p0 << " p1 " << p1 << " t " << t << " x " << m_x << "\n";
        
        dxdp0 = t*(Matrix3d::eye() +
            outerProduct(p0-p1, plane.normal()) / dot(plane.normal(), p1 - p0));
        
        dxdp1 = (1-t)*(Matrix3d::eye() +
            outerProduct(p0-p1, plane.normal()) / dot(plane.normal(), p1 - p0));
        
        vNum0 = edge;
        vNum1 = (edge+1)%3;
        vNum2 = (edge+2)%3;
        
//        cerr << "dxdp" << vNum0 << " = " << dxdp0 << "\n";
//        cerr << "dxdp" << vNum1 << " = " << dxdp1 << "\n";
        
        mDx[vNum0] = dxdp0;
        mDx[vNum1] = dxdp1;
        mDx[vNum2] = Matrix3d::zero();
        
        foundPoint = true;
        //cerr << "Found point on edge " << edge << "\n";
    }
    
    if (foundPoint == false)
    {
        cerr << "Triangle is " << tri << "\n";
        cerr << "Plane is " << plane << "\n";
        throw(std::logic_error("Failed to find an intersection"));
    }
    
    // Choose the direction as cross(nTri, nPlane).
    // Its sensitivity is DnTri x nPlane.
    
    m_v = cross(unit(tri.normal()), plane.normal());
    
    // Sensitivity of the line direction w.r.t. triangle vertices.
    
    for (int triVert = 0; triVert < 3; triVert++)
    {
        Matrix3d dudv = tri.dudv(triVert);
        
        for (int direction = 0; direction < 3; direction++)
        {
            // dvdvert is the gradient component dv/dvert_direction.
            
            Vector3d dvdvert = cross(dudv.column(direction), plane.normal());
            
            for (int row = 0; row < 3; row++)
            {
                mDv[triVert](row,direction) = dvdvert[row];
            }
        }
    }
    
    
//    cerr << "Dx =\n" << mDx[0] << "\n" << mDx[1] << "\n" << mDx[2] << "\n";
    /*
//    cerr << "x0 = " << x[0] << "\n";
//    cerr << "p0 = " << whichP0[0] << " p1 = " << whichP1[0] << "\n";
//    cerr << "dx0dp0 = " << dxdp0[0] << "\n";
//    cerr << "dx0dp1 = " << dxdp1[0] << "\n";
//    cerr << "x1 = " << x[1] << "\n";
//    cerr << "p0 = " << whichP0[1] << " p1 = " << whichP1[1] << "\n";
//    cerr << "dx1dp0 = " << dxdp0[1] << "\n";
//    cerr << "dx1dp1 = " << dxdp1[1] << "\n";
    
    m_v = x[1] - x[0];
    
    // The first edge (edge 0) determines x and Dx.
    mDx[whichP0[0]] = dxdp0[0];
    mDx[whichP1[0]] = dxdp1[0];
    
    // The second edge (edge 1) and the first edge (edge 0) determine v and Dv
    mDv[whichP0[0]] -= dxdp0[0];
    mDv[whichP1[0]] -= dxdp1[0];
    mDv[whichP0[1]] += dxdp0[1];
    mDv[whichP1[1]] += dxdp1[1];
    */
    
    if (sumSquares(m_v) == 0)
    {
        cerr << "Warning: Line has zero velocity.  This probably means the triangle"
            " intersects the plane at a point.\n";
    }
}

// Given an edge of tri0, find its supporting line and intersect it with
// the supporting plane of tri1.
// 
// To do so, express the line as (e.g.) line(t) = tri0[0] + t*(tri0[1]-tri0[0]).
// This must equal tri1[0] + a*(tri1[1]-tri1[0]) + b*(tri1[2]-tri1[0]).
// Solve for a, b and t simultaneously.  Below, "coordinates" = [a b -t].  (Note
// the negative sign on -t.)
//
// The calculus to obtain dIntersection follows naturally.
//
Vector3d Line::
intersection(const Triangle3d & tri0, const Triangle3d & tri1, int whichLine)
{
    int vert0 = whichLine, vert1 = (whichLine+1)%3;
    
    Matrix3d intersectionBasis = Matrix3d::withColumns(
        tri1[1]-tri1[0], tri1[2]-tri1[0], tri0[vert1]-tri0[vert0]);
    
    Vector3d coordinates = inverse(intersectionBasis)*(tri0[vert0] - tri1[0]);
    
    Vector3d point = tri0[vert0] - coordinates[2]*(tri0[vert1]-tri0[vert0]);
    
    return point;
}

Vector3d Line::
dIntersection(const Triangle3d & tri0, const Triangle3d & tri1,
    const Triangle3d & dTri0, const Triangle3d & dTri1, int whichLine)
{
    int vert0 = whichLine, vert1 = (whichLine+1)%3;
    
    Matrix3d intersectionBasis = Matrix3d::withColumns(
        tri1[1]-tri1[0], tri1[2]-tri1[0], tri0[vert1]-tri0[vert0]);
    Matrix3d inverseBasis = inverse(intersectionBasis);
    
    Matrix3d dBasis = Matrix3d::withColumns(
        dTri1[1]-dTri1[0], dTri1[2]-dTri1[0], dTri0[vert1]-dTri0[vert0]);
    
    Vector3d coordinates = inverseBasis*(tri0[vert0] - tri1[0]);
    Vector3d dCoordinates = inverseBasis*(dTri0[vert0] - dTri1[0]
        - dBasis*coordinates);
    
    Vector3d dPoint = dTri0[vert0] - dCoordinates[2]*(tri0[vert1] - tri0[vert0])
        - coordinates[2] * (dTri0[vert1] - dTri0[vert0]);
    
    return dPoint;
}

ostream & operator<<(ostream & str, const Line & l)
{
    str << l.x() << " along " << l.v();
    return str;
}

// On initial TP()
SensitiveEdge::
SensitiveEdge(const ControlVertex & v0, const ControlVertex & v1) :
    mLine(v0.point(), v1.point()),
    mControlVertices(2)
{
    mControlVertices[0] = v0;
    mControlVertices[1] = v1;
//    if (sumSquares(line().v()) == 0)
//        int what = true;
//    assert(sumSquares(line().v()) > 0);
}


// Edge along intersection of a triangle with a plane
SensitiveEdge::SensitiveEdge(const Triangle & tri, const Plane3d & plane):
    mLine(tri.triangle(), plane),
    mControlVertices(3)
{
    for (int vv = 0; vv < 3; vv++)
    {
        mControlVertices[vv] = tri.controlVertices()[vv];
    }
}

SensitiveEdge::SensitiveEdge(const Triangle & tri0, const Triangle & tri1):
    mLine(tri0.triangle(), tri1.triangle()),
    mControlVertices(6)
{
    for (int vv = 0; vv < 3; vv++)
    {
        mControlVertices[vv] = tri0.controlVertices()[vv];
        mControlVertices[vv+3] = tri1.controlVertices()[vv];
    }
}

bool SensitiveEdge::
isSensitive() const
{
    for (int vv = 0; vv < mControlVertices.size(); vv++)
    if (mControlVertices[vv].freeDirections() != Vector3b(0,0,0))
        return true;
    return false;
}

Triangle::
Triangle() :
    mTriangle(),
    mUnitNormal(0,0,0),
    mNeighbors(-1,-1,-1),
    mControlVertices()
{
}

Triangle::
Triangle(const Triangle3d & t) :
    mTriangle(t),
    mUnitNormal(unit(t.normal())),
    mNeighbors(-1,-1,-1),
    mControlVertices()
{
}


Triangle::
Triangle(Vector3i idxVert, const std::vector<ControlVertex> & cvs) :
    mTriangle(cvs[idxVert[0]].point(), cvs[idxVert[1]].point(), cvs[idxVert[2]].point()),
    mUnitNormal(unit(mTriangle.normal())),
    mNeighbors(-1,-1,-1),
    mControlVertices(cvs[idxVert[0]], cvs[idxVert[1]], cvs[idxVert[2]])
{
}

Matrix3d Triangle::
nnT() const
{
    return outerProduct(unitNormal(), unitNormal());
}

Matrix3d Triangle::
nnTdA() const
{
    return outerProduct(unitNormal(), unitNormal()) /
        fabs(unitNormal()[upAxis()]);
}


Triangle3d Triangle::
controlTriangle() const
{
    return Triangle3d(mControlVertices[0].point(),
        mControlVertices[1].point(),
        mControlVertices[2].point());
}


void Triangle::
edgeControlVertices(const Vector3<ControlVertex> & verts)
{
    for (int edge = 0; edge < 3; edge++)
    {
        mEdges[edge] = SensitiveEdge(verts[edge], verts[(edge+1)%3]);
        
//        mEdges[edge].controlVertices().resize(2);
//        mEdges[edge].controlVertices()[0] = verts[edge];
//        mEdges[edge].controlVertices()[1] = verts[(edge+1)%3];
    }
}

void Triangle::
edgeControlVertices(const ControlVertex & v0, const ControlVertex & v1,
    const ControlVertex & v2)
{
    edgeControlVertices(Vector3<ControlVertex>(v0, v1, v2));
}

void Triangle::
setEdgeControlVertices(int edge, const Triangle & neighborTri)
{
    // Nowwwwww then.  Every edge has two triangles as neighbors.
    // Every triangle has a control triangle.
    //
    // If a triangle is coplanar with its neighbor, then there is no
    // reason to use the derivative of the edge because it won't change
    // any surface areas.
    //
    // If a triangle is not coplanar with its neighbor, then we can
    // just calculate its sensitivity...
    
    Vector3d nv = triangle().normal();
    Vector3d nv2 = neighborTri.triangle().normal();
    
    // If coplanar:
    if ( fabs(dot(unit(nv), unit(nv2))) > (1.0 - 1e-8) )
    {
        mEdges[edge] = SensitiveEdge();
    }
    else
    {
        mEdges[edge] = SensitiveEdge(*this, neighborTri);
    }
}

void Triangle::
setEdgeControlVertices(int edge, const ControlVertex & cv0, const ControlVertex & cv1)
{
    mEdges[edge] = SensitiveEdge(cv0, cv1);
}

void Triangle::
cacheSensitivity()
{   
    // Obtain the Jacobian of the unit normal vector.
    
    Triangle3d controlTri(mControlVertices[0].point(),
        mControlVertices[1].point(),
        mControlVertices[2].point());
    
//    Vector3d u = unitNormal();
    Vector3d u = unit(controlTri.normal());
    double nz = fabs(u[upAxis()]);    
//    cerr << "Control tri " << controlTri << "\n";
//    cerr << "normal " << u << " axis " << upAxis() << "\n";
//    cerr << "My normal " << unitNormal() << "\n";
//    cerr << "Control tri normal " << u << "\n";
    assert(nz != 0);
    
    
//    Matrix3d nnT = outerProduct(u, u);
    for (int vert = 0; vert < 3; vert++)
    {
        Matrix3<Vector3d> DOrientation; // full sensitivity of nnT/nz
        
        Matrix3d dudv = controlTri.dudv(vert);
        Vector3d dnzdv = dudv.row(upAxis()) * (u[upAxis()] < 0.0 ? -1.0 : 1.0);
        
        for (int kk = 0; kk < 3; kk++)
        for (int jj = 0; jj < 3; jj++)
        for (int ii = 0; ii < 3; ii++)
        {
            DOrientation(ii,jj)[kk] = (1.0/nz/nz) * (
                nz*(dudv(ii,kk)*u[jj] + u[ii]*dudv(jj,kk)) -
                u[ii]*u[jj]*dnzdv[kk]);
//            cerr << Vector3i(ii,jj,kk) << ": " << DOrientation(ii,jj)[kk]
//                << "\n";
        }
        
//        cerr << "\tdnnT/dv" << vert << " =\n" << DOrientation << "\n";
        mDOrientation[vert] = DOrientation;
    }
    
    // Obtain the Jacobians of each edge that's not between coplanar triangles.
    
    for (int edge = 0; edge < 3; edge++)
    {
        // I allow edges w/o control vertices between coplanar triangles!
//        if (edgeControlVertices(edge).size() == 0)
//            throw(std::logic_error("Edge lacks control vertices"));
        
        if (edgeControlVertices(edge).size() == 2)
        {
            mEdges[edge].line() = Line(
                edgeControlVertices(edge)[0].point(),
                edgeControlVertices(edge)[1].point());
        }
        else if (edgeControlVertices(edge).size() == 6)
        {
            Triangle3d controlTri0(
                edgeControlVertices(edge)[0].point(),
                edgeControlVertices(edge)[1].point(),
                edgeControlVertices(edge)[2].point());
            Triangle3d controlTri1(
                edgeControlVertices(edge)[3].point(),
                edgeControlVertices(edge)[4].point(),
                edgeControlVertices(edge)[5].point());
            
            Vector3d n0 = unit(controlTri0.normal());
            Vector3d n1 = unit(controlTri1.normal());
            const double THRESHOLD = 1e-8;
            
            if (fabs(dot(n0, n1)) <= (1.0 - THRESHOLD)) // if not coplanar
            {
                mEdges[edge].line() = Line(controlTri0, controlTri1);
            }
            else // suppress coplanar triangle control verts
            {
                mEdges[edge].controlVertices().clear();
            }
        }
        else if (edgeControlVertices(edge).size() != 0)
        {
            throw(std::logic_error("Edge has wrong number of "
                "control vertices"));
        }
    }
}

bool Triangle::
isSensitive() const
{
    for (int vv = 0; vv < 3; vv++)
    if (mControlVertices[vv].freeDirections() != Vector3b(0,0,0))
        return true;
    
    for (int ee = 0; ee < 3; ee++)
    if (mEdges[ee].isSensitive())
        return true;
    
    return false;
}

std::ostream & operator<<(std::ostream & str, const Triangle & tri)
{
    str << "CVs: " << tri.controlVertices()[0].id() << ", "
        << tri.controlVertices()[1].id() << ", "
        << tri.controlVertices()[2].id() << "; ";
    str << "tri: " << tri.triangle() << "; ";
    
    for (int ee = 0; ee < 3; ee++)
    {
        const std::vector<ControlVertex> & cvs = tri.edgeControlVertices(ee);
        str << "e" << ee << ":";
        for (int vv = 0; vv < cvs.size(); vv++)
        {
            str << cvs[vv].id();
            if (vv < cvs.size()-1)
                str << " ";
        }
        if (ee < 2)
        {
            str << ", ";
        }
    }
    return str;
}

namespace Test
{

Triangle oneTri(Vector3d p0, Vector3d p1, Vector3d p2)
{
    Vector3b threeTrues(true, true, true);
    Triangle3d tri(p0, p1, p2);

    Vector3<ControlVertex> controlVerts;
    controlVerts[0] = ControlVertex(p0, threeTrues, 0);
    controlVerts[1] = ControlVertex(p1, threeTrues, 1);
    controlVerts[2] = ControlVertex(p2, threeTrues, 2);
    
    Triangle outTri(Triangle3d(p0, p1, p2));
    outTri.controlVertices(controlVerts);
    outTri.edgeControlVertices(outTri.controlVertices());
    outTri.cacheSensitivity();
    
    return outTri;
}

Triangle oneFlippedTri(Vector3d p0, Vector3d p1, Vector3d p2)
{
    Vector3b threeTrues(true, true, true);

    Vector3<ControlVertex> controlVerts;
    controlVerts[0] = ControlVertex(p0, threeTrues, 0);
    controlVerts[1] = ControlVertex(p1, threeTrues, 1);
    controlVerts[2] = ControlVertex(p2, threeTrues, 2);
    
    Triangle outTri(Triangle3d(p0, p2, p1));
    outTri.controlVertices(controlVerts);
    outTri.edgeControlVertices(outTri.controlVertices());
    outTri.cacheSensitivity();
    
    return outTri;
}

Triangle3d perturbTri(int vert, int dir)
{
    Triangle3d t(Vector3d(0,0,0), Vector3d(0,0,0), Vector3d(0,0,0));
    t[vert][dir] = 1.0;
    return t;
}




void oneBoxMesh(std::vector<SimpleMesh::ControlVertex> & controlVertices,
    std::vector<SimpleMesh::Triangle> & meshTris,
    Rect3d r, Matrix3d tx, Vector3d translation)
{
    controlVertices.resize(8);
    meshTris.resize(12);
//    mesh.controlVertexIds.resize(8);
//    mesh.freeDirections.resize(8, Vector3b(true,true,true));
//    mesh.faces.resize(12);
    
    double x0 = r.p1[0], x1 = r.p2[0],
        y0 = r.p1[1], y1 = r.p2[1],
        z0 = r.p1[2], z1 = r.p2[2];
    
    controlVertices[0] = ControlVertex(tx*Vector3d(x0, y0, z0) + translation, 0);
    controlVertices[1] = ControlVertex(tx*Vector3d(x1, y0, z0) + translation, 1);
    controlVertices[2] = ControlVertex(tx*Vector3d(x0, y1, z0) + translation, 2);
    controlVertices[3] = ControlVertex(tx*Vector3d(x1, y1, z0) + translation, 3);
    controlVertices[4] = ControlVertex(tx*Vector3d(x0, y0, z1) + translation, 4);
    controlVertices[5] = ControlVertex(tx*Vector3d(x1, y0, z1) + translation, 5);
    controlVertices[6] = ControlVertex(tx*Vector3d(x0, y1, z1) + translation, 6);
    controlVertices[7] = ControlVertex(tx*Vector3d(x1, y1, z1) + translation, 7);
    for (int vv = 0; vv < controlVertices.size(); vv++)
    {
        controlVertices[vv].freeDirections(Vector3b(true, true, true));
    }
    
    
    meshTris[0] = Triangle(Vector3i(0, 4, 6), controlVertices); // x-
    meshTris[1] = Triangle(Vector3i(0, 6, 2), controlVertices);
    meshTris[2] = Triangle(Vector3i(1, 3, 7), controlVertices); // x+
    meshTris[3] = Triangle(Vector3i(1, 7, 5), controlVertices);
    meshTris[4] = Triangle(Vector3i(0, 1, 5), controlVertices); // y-
    meshTris[5] = Triangle(Vector3i(0, 5, 4), controlVertices);
    meshTris[6] = Triangle(Vector3i(2, 6, 7), controlVertices); // y+
    meshTris[7] = Triangle(Vector3i(2, 7, 3), controlVertices);
    meshTris[8] = Triangle(Vector3i(0, 2, 3), controlVertices); // z-
    meshTris[9] = Triangle(Vector3i(0, 3, 1), controlVertices);
    meshTris[10] = Triangle(Vector3i(4, 5, 7), controlVertices); // z+
    meshTris[11] = Triangle(Vector3i(4, 7, 6), controlVertices);
    
    initMeshSensitivity(meshTris);
}

vector<Vector3d> subtractPositions(const vector<ControlVertex> & v1,
    const vector<ControlVertex> & v2)
{
    assert(v1.size() == v2.size());
    vector<Vector3d> vOut(v1.size());
    
    for (int nn = 0; nn < v1.size(); nn++)
    {
        vOut[nn] = v1[nn].point() - v2[nn].point();
        
        assert(v1[nn].id() == v2[nn].id());
        //assert(v1[nn].id() == nn);
    }
    
    return vOut;
}

}; // namespace Test


}; // namespace SimpleMesh
