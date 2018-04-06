/*
 *  SimpleMesh.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 9/8/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _SIMPLEMESH_
#define _SIMPLEMESH_

#include "geometry.h"
#include "geometry2.h"

/*
    SimpleMesh is an inefficient but simple representation of a triangulated
    mesh that arises from Boolean operations on polyhedra.  Its main classes are
    Triangle and ControlVertex.  ControlVertex is a representation of a vertex
    of one of the input polyhedra and contains the coordinates and free
    directions (directions of possible motion) of each vertex.  Triangle is a
    triangle arising from triangulation of the mesh after Boolean operations.
    Each triangle knows three ControlVertices which define its plane of support.
    (There may be more than three such vertices, and the three chosen are
    somewhat arbitrary.)
    
    A polygonal face of the mesh will change shape when the input polyhedra
    change.  The edges of these polygons may coincide with edges of input
    polyhedra, in which case they have two control vertices.  On the other hand,
    the edges may arise from one polyhedron cleaving another, in which case the
    edges lie along the intersection of two polygons' planes of support.  Such
    edges have as many as six control vertices (three each for two planes).
    A triangulated polygon adds internal edges as well.  Such edges also move
    when the input polyhedra move, but their sensitivity is not important for
    calculating the sensitivity of volumes and normal vectors in voxels, so they
    are not given any sensitivity information in SimpleMesh.  Internal edges may
    be recognized by their neighboring triangles, which are coplanar; ONLY
    internal edges neighbor coplanar triangles.
 
    This is all madness because it should have been a face-vertex data structure
    in the first place.  Should be replaced really with something more normal
    and less rigid and wasteful.
*/

namespace SimpleMesh
{

class Triangle;

class ControlVertex
{
public:
    ControlVertex() :
        mId(-1),
        mPoint(),
        mFreeDirections()
    {}
    
    ControlVertex(const Vector3d & point) :
        mId(-1),
        mPoint(point),
        mFreeDirections()
    {}
    
    ControlVertex(const Vector3d & point, int id) :
        mId(id),
        mPoint(point),
        mFreeDirections()
    {}
    
    ControlVertex(const Vector3d & point, const Vector3b & freeDirections) : // unused
        mId(-1),
        mPoint(point),
        mFreeDirections(freeDirections)
    {}
    
    // This function is used when creating a TrackedPolyhedron.
    // This is the only actually-useful constructor I think
    ControlVertex(const Vector3d & point, const Vector3b & freeDirections, // NefBuilder::controlVertices()
        int id) :
        mId(id),
        mPoint(point),
        mFreeDirections(freeDirections)
    {}
    
    int id() const { return mId; }
    const Vector3d & point() const { return mPoint; }
    const Vector3b & freeDirections() const { return mFreeDirections; }
    
//    Vector3d & point() { return mPoint; }
//    Vector3b & freeDirections() { return mFreeDirections; }
    
    void id(int newId) { mId = newId; }
    void point(const Vector3d & p) { mPoint = p; }
    void freeDirections(const Vector3b & f) { mFreeDirections = f; }
    
private:
//    static int sNextId;
    int mId;
    Vector3d mPoint;
    Vector3b mFreeDirections;
};
std::ostream & operator<<(std::ostream & str, const ControlVertex & cv);

// A line supporting two points or supporting the intersection of two triangles.
// Include as well its sensitivity to any perturbations of these vertices.
class Line
{
public:
    Line();
    
    // Line with two control vertices (internally, indexed 0 and 1).
    // The line is given by the point vert0 and the direction (vert1-vert0).
    // The 2 Jacobians are indexed as
    //  Dx[0] = dx/d(vert0)
    //  Dx[1] = dx/d(vert1)
    // and similarly for Dv.
    Line(const Vector3d & vert0, const Vector3d & vert1);
    
    // Line supporting the intersection of tri0 and tri1.
    // A point x and direction v along the line are determined.
    // The 6 Jacobians are indexed as
    //   Dx[0] = dx/d(control vert 0) of triangle 0
    //   Dx[1] = dx/d(control vert 1) of triangle 0
    //   Dx[2] = dx/d(control vert 2) of triangle 0
    //   Dx[3] = dx/d(control vert 0) of triangle 1
    //   Dx[4] = dx/d(control vert 1) of triangle 1
    //   Dx[5] = dx/d(control vert 2) of triangle 1
    // and similarly for Dv.
    Line(const Triangle3d & tri0, const Triangle3d & tri1);
    
    // Line supporting the intersection of tri with plane.
    // A point x and direction v along the line are determined.
    // The 3 Jacobians are indexed as
    //  Dx[i] = dx/d(control vert i)
    // and similarly for Dv.
    Line(const Triangle3d & tri, const Plane3d & plane);
    
    const Vector3d & x() const { return m_x; } // a point on the line
    const Vector3d & v() const { return m_v; } // the direction of the line
    
    const Matrix3d & Dx(unsigned int controlVertex) const
        { return mDx[controlVertex]; }
    const Matrix3d & Dv(unsigned int controlVertex) const
        { return mDv[controlVertex]; }
    
    int numControlVertices() const { return mNumControlVertices; }
    
    // Expose these helper functions for unit testing.
    
    // Given an edge of tri0, find its supporting line and intersect it with
    // the supporting plane of tri1.
    static Vector3d intersection(const Triangle3d & tri0,
        const Triangle3d & tri1, int whichLine);
    
    // Sensitivity of result of intersection() to perturbation of tri0 or tri1.
    static Vector3d dIntersection(const Triangle3d & tri0,
        const Triangle3d & tri1, const Triangle3d & dTri0,
        const Triangle3d & dTri1, int whichLine);
private:
    
    Vector3d m_x;
    Vector3d m_v;
    
    int mNumControlVertices; // zero to six
    
    Matrix3d mDx[6];
    Matrix3d mDv[6];
};
std::ostream & operator<<(std::ostream & str, const Line & l);

/**
 * Representation of an edge of a mesh and its derivatives as the mesh changes.
 * It's just a container for a Line (which does calculus) and a vector of
 * ControlVertex records (which cache a lot of things).
 */
class SensitiveEdge
{
public:
    SensitiveEdge() : mLine(), mControlVertices() {}
   
    // Edge between two vertices.
    // Called when constructing a TrackedPolyhedron.
    SensitiveEdge(const ControlVertex & v0, const ControlVertex & v1);
    
    // Edge along intersection of a triangle with a plane.
    // Called during voxelization.
    SensitiveEdge(const Triangle & tri, const Plane3d & plane);
    
    SensitiveEdge(const Triangle & tri0, const Triangle & tri1);
    
    const Line & line() const { return mLine; }
    const std::vector<ControlVertex> & controlVertices() const
        { return mControlVertices; }
    
    Line & line() { return mLine; }
    std::vector<ControlVertex> & controlVertices() { return mControlVertices; }
    
    bool isSensitive() const;
private:
    Line mLine;
    std::vector<ControlVertex> mControlVertices;
};

class Triangle
{
public:
    Triangle();
    Triangle(const Triangle3d & t);
    Triangle(Vector3i idxVert, const std::vector<ControlVertex> & cvs);
    
    Vector3d unitNormal() const { return mUnitNormal; }
    Matrix3d nnT() const;
    Matrix3d nnTdA() const;
    
    Matrix3<Vector3d> Dorientation(unsigned int controlVertex) const
        { return mDOrientation[controlVertex]; }
    int upAxis() const { return dominantDirection(unitNormal()); }
    
    const SensitiveEdge & edgeSensitivity(int edge) const
        { return mEdges[edge]; }
    
    const Triangle3d & triangle() const { return mTriangle; } //Assembly:writeForMatlab, ClippedTriangle(), handleSortedRecords, voxelizeOneTri(), boxHits(), clipToRow(), assignControlVertices(), findNeighborTriangles_fast/_slow(), inheritControlVertices, tryInheritingControlVertices(), TriangleInBox(), orientationJacobians(), orientations(), linesInBounds()
    Triangle3d controlTriangle() const;
    const Vector3i & neighbors() const { return mNeighbors; }

    const Vector3<ControlVertex> & controlVertices() const
        { return mControlVertices; }
    const std::vector<ControlVertex> & edgeControlVertices(int ee) const
        { return mEdges[ee].controlVertices(); }
    
    void triangle(const Triangle3d & t)
        { mTriangle = t; mUnitNormal = unit(t.normal()); }
    void neighbors(const Vector3i & neighbors) { mNeighbors = neighbors; }
    void controlVertices(const Vector3<ControlVertex> & verts)
        { mControlVertices = verts; }
    void controlVertices(const ControlVertex & v0, const ControlVertex & v1,
        const ControlVertex & v2)
    {
        mControlVertices[0] = v0;
        mControlVertices[1] = v1;
        mControlVertices[2] = v2;
    }
    
    void edgeControlVertices(int edge, const std::vector<ControlVertex> & cvs)
        { mEdges[edge].controlVertices() = cvs; }
    void edgeControlVertices(const Vector3<ControlVertex> & verts); // pass control tri
    void edgeControlVertices(const ControlVertex & v0, const ControlVertex & v1,
        const ControlVertex & v2);
    
    void setEdgeControlVertices(int edge, const Triangle & neighborTri);
    void setEdgeControlVertices(int edge, const ControlVertex & cv0, const ControlVertex & cv1);
    
    void neighbor(int which, int neighbor) { mNeighbors[which] = neighbor; }
    
    void cacheSensitivity();
    bool isSensitive() const; // quick way to check whether it has any Jacobians
    
private:
    Triangle3d mTriangle;
    Vector3d mUnitNormal;
    Vector3i mNeighbors;
    
    // Control vertex information for this tri and its edges
    Vector3<ControlVertex> mControlVertices;
    
    // Cached sensitivity of outerProduct(n,n)/nz, where the z-axis is upAxis().
    Matrix3<Vector3d> mDOrientation[3];
    
    SensitiveEdge mEdges[3];
};
std::ostream & operator<<(std::ostream & str, const Triangle & tri);


namespace Test
{

// Create a triangle that is its own control triangle.
Triangle oneTri(Vector3d p0, Vector3d p1, Vector3d p2);

// Just have the triangle facing oppositely to its control tri.  Important
// case!
Triangle oneFlippedTri(Vector3d p0, Vector3d p1, Vector3d p2);
Triangle3d perturbTri(int vert, int dir);

void oneBoxMesh(std::vector<SimpleMesh::ControlVertex> & controlVertices,
    std::vector<SimpleMesh::Triangle> & meshTris,
    Rect3d r, Matrix3d tx = Matrix3d::eye(), Vector3d translation = Vector3d(0,0,0));
    
    
std::vector<Vector3d> subtractPositions(
    const std::vector<SimpleMesh::ControlVertex> & v1,
    const std::vector<SimpleMesh::ControlVertex> & v2);

}; // namespace Test


}; // namespace SimpleMesh

//#include "SimpleMesh-inl.h"

#endif
