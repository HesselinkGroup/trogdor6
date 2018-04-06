/*
 *  Search.cpp
 *  SimpleMesh
 *
 *  Created by Paul C Hansen on 6/27/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "Search.h"
#include <iostream>
#include <map>

using namespace std;

namespace SimpleMesh
{

AllControlVerticesIn::
AllControlVerticesIn(int cv0, int cv1, int cv2, int cv3,
    int cv4, int cv5, int cv6, int cv7) :
    mControlVertices()
{
    mControlVertices.insert(cv0);
    mControlVertices.insert(cv1);
    if (cv2 != -1) mControlVertices.insert(cv2);
    if (cv3 != -1) mControlVertices.insert(cv3);
    if (cv4 != -1) mControlVertices.insert(cv4);
    if (cv5 != -1) mControlVertices.insert(cv5);
    if (cv6 != -1) mControlVertices.insert(cv6);
    if (cv7 != -1) mControlVertices.insert(cv7);
}

bool AllControlVerticesIn::
operator()(const std::vector<ControlVertex> & lineCVs) const
{
    if (mControlVertices.size() > 2 && lineCVs.size() < 6)
        return false;
    if (lineCVs.size() == 0)
        return false;
    
    for (int nn = 0; nn < lineCVs.size(); nn++)
    {
        if (mControlVertices.count(lineCVs[nn].id()) == 0)
            return false;
    }
    return true;
}

bool AllControlVerticesIn::
operator()(const Triangle & triangle) const
{
    for (int vv = 0; vv < 3; vv++)
    if (mControlVertices.count(triangle.controlVertices()[vv].id()) == 0)
        return false;
    
    return true;
}





std::vector<SensitiveEdge>
linesWithControlVertices(const Triangle & triangle,
    int cv0, int cv1, int cv2, int cv3,
    int cv4, int cv5, int cv6, int cv7)
{
    vector<SensitiveEdge> out;
    AllControlVerticesIn predicate(cv0, cv1, cv2, cv3, cv4, cv5, cv6, cv7);
    
    for (int edge = 0; edge < 3; edge++)
    if (predicate(triangle.edgeControlVertices(edge)))
        out.push_back(triangle.edgeSensitivity(edge));
    
    return out;
}

std::vector<SensitiveEdge>
linesWithControlVertices(const std::vector<Triangle> & triangles,
    int cv0, int cv1, int cv2, int cv3,
    int cv4, int cv5, int cv6, int cv7)
{
    vector<SensitiveEdge> out;
    AllControlVerticesIn predicate(cv0, cv1, cv2, cv3, cv4, cv5, cv6, cv7);
    
    for (int tt = 0; tt < triangles.size(); tt++)
    for (int edge = 0; edge < 3; edge++)
    if (predicate(triangles[tt].edgeControlVertices(edge)))
    {
        out.push_back(triangles[tt].edgeSensitivity(edge));
    }
    
    return out;
}



std::vector<Triangle>
trianglesWithControlVertices(const std::vector<Triangle> & triangles,
    int cv0, int cv1, int cv2, int cv3)
{
    vector<Triangle> out;
    AllControlVerticesIn predicate(cv0, cv1, cv2, cv3);
    
    for (int tt = 0; tt < triangles.size(); tt++)
    if (predicate(triangles[tt]))
        out.push_back(triangles[tt]);
    
    return out;
}


std::vector<SensitiveEdge>
linesInBounds(const std::vector<Triangle> & triangles, Rect3d bounds)
{
    vector<SensitiveEdge> out;
    
    for (int tt = 0; tt < triangles.size(); tt++)
    for (int edge = 0; edge < 3; edge++)
    {
        if (bounds.encloses(triangles[tt].triangle()[edge]) &&
            bounds.encloses(triangles[tt].triangle()[(edge+1)%3]))
        {
            out.push_back(triangles[tt].edgeSensitivity(edge));
        }
    }
    
    return out;
}

std::vector<Triangle>
trianglesInBounds(const std::vector<Triangle> & triangles, Rect3d bounds)
{
    vector<Triangle> out;
    
    for (int tt = 0; tt < triangles.size(); tt++)
    if (bounds.encloses(triangles[tt].triangle()))
        out.push_back(triangles[tt]);
    
    return out;
}




//typedef Vector2<Vector3d> UEdge;
//static UEdge sUnorientedEdge(const SimpleMesh::Triangle & tri, int edge)
//{
//    const Vector3d & p0 = tri.triangle()[edge];
//    const Vector3d & p1 = tri.triangle()[(edge+1)%3];
//
//    if (p0 <= p1)
//    {
//        return UEdge(p0, p1);
//    }
//    else
//    {
//        return UEdge(p1, p0);
//    }
//}

typedef Vector2<Vector3d> Edge;
static Edge fwEdge(const SimpleMesh::Triangle & tri, int edge)
{
    return Edge(tri.triangle()[edge], tri.triangle()[(edge+1)%3]);
}

typedef Vector2<Vector3d> Edge;
static Edge bkEdge(const SimpleMesh::Triangle & tri, int edge)
{
    return Edge(tri.triangle()[(edge+1)%3], tri.triangle()[edge]);
}

void determineNeighbors(std::vector<Triangle> & triangles)
{
    std::map<Edge, int> edge2tri;
    
    // 1. Build edge2tri map
    for (int tt = 0; tt < triangles.size(); tt++)
    {
        const SimpleMesh::Triangle & tri = triangles[tt];
        for (int ee = 0; ee < 3; ee++)
        {
            edge2tri[fwEdge(tri,ee)] = tt;
        }
    }
    
    // 2. Find neighbors
    std::map<Edge,int>::const_iterator itr;
    
    for (int tt = 0; tt < triangles.size(); tt++)
    {
        SimpleMesh::Triangle & tri = triangles[tt];
        for (int ee = 0; ee < 3; ee++)
        {
            itr = edge2tri.find(bkEdge(tri,ee));
            if (itr == edge2tri.end())
            {
                throw std::runtime_error("Triangle has no neighbor");
            }
            else
            {
                tri.neighbor(ee, itr->second);
            }
        }
    }
}

static bool sEdgesAreCollinear(const Vector3d & v0, const Vector3d & v1,
    const Vector3d & w0, const Vector3d & w1,
    double dist)
{
    double distSquared = dist*dist;
    if (pointLineDistanceSquared(v0, w0, w1) < distSquared &&
        pointLineDistanceSquared(v1, w0, w1) < distSquared)
    {
        return true;
    }
    return false;
}

void determineEdgeControlVertices(std::vector<Triangle> & triangles)
{
    for (int tt = 0; tt < triangles.size(); tt++)
    {
        Triangle & tri = triangles.at(tt);
        
        for (int idxEdge = 0; idxEdge < 3; idxEdge++)
        {
            // If the edge is collinear with the line of two of its tri's CVs,
            // then the edge just runs between the two CVs.  Otherwise it has
            // all six CVs of its two adjacent triangles.
            
            bool isCollinearWithCVs = false;
            for (int idxCV = 0; idxCV < 3; idxCV++)
            {
                int idxVert = idxEdge;
                int idxVertNext = (idxEdge+1)%3;
                int idxCVNext = (idxCV+1)%3;
                
                const Vector3d & cv0 = tri.controlTriangle()[idxCV];
                const Vector3d & cv1 = tri.controlTriangle()[idxCVNext];
                const Vector3d & v0 = tri.triangle()[idxVert];
                const Vector3d & v1 = tri.triangle()[idxVertNext];
                
                if (sEdgesAreCollinear(v0, v1, cv0, cv1, 1e-8))
                {
                    tri.setEdgeControlVertices(idxEdge, tri.controlVertices()[idxCV], tri.controlVertices()[idxCVNext]);
                    isCollinearWithCVs = true;
                    break;
                }
            }
            
            if (!isCollinearWithCVs)
            {
                Triangle & neighbor = triangles.at(tri.neighbors()[idxEdge]);
                tri.setEdgeControlVertices(idxEdge, neighbor);
            }
        }
    }
}


void initMeshSensitivity(std::vector<Triangle> & triangles)
{
    determineNeighbors(triangles);
    determineEdgeControlVertices(triangles);
    for (int tt = 0; tt < triangles.size(); tt++)
    {
        triangles[tt].cacheSensitivity();
    }
}

void initMeshSensitivity(std::vector<std::vector<Triangle> > & meshes)
{
    for (int mm = 0; mm < meshes.size(); mm++)
    {
        initMeshSensitivity(meshes.at(mm));
    }
}









} // namespace SimpleMesh

