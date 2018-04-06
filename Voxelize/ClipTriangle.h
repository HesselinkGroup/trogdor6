/*
 *  ClipTriangle.h
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 6/28/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 *
 */

#ifndef CLIPTRIANGLE_H
#define CLIPTRIANGLE_H

#include "SimpleMesh/SimpleMesh.h"
#include "geometry.h"
#include <map>

// Representation of a convex polygon arising from clipping a triangle along
// an axis-oriented plane.
//
// Boundary cases:
//
// Isolated vertices (one vertex, one edge) are EXCLUDED
// Isolated edges (two vertices, one edge) are INCLUDED because they're needed
// for taking the derivative of the triangle's area.
//
class ClippedTriangle
{
public:
    ClippedTriangle();
    ClippedTriangle(const SimpleMesh::Triangle & tri);
    
    const std::vector<Vector3d> & vertices() const { return mVertices; }
    const std::vector<const SimpleMesh::SensitiveEdge*> & edges() const
        { return mEdges; }
    
    const SimpleMesh::Triangle* triangle() const { return mTriangle; }
    
    // Points < val are below, points >= val are above.
    // Polygons of less than three vertices will be culled.
    ClippedTriangle partBelow(int axis, double val,
        const SimpleMesh::SensitiveEdge* clippedSensitivity = 0L) const;
    ClippedTriangle partAbove(int axis, double val,
        const SimpleMesh::SensitiveEdge* clippedSensitivity = 0L) const;
    
    ClippedTriangle partInside(const Rect2d & xyBounds) const;

    
    Rect3d bounds() const;
    
    // Return the area perpendicular to the given axis
    double area2d(int axis) const;
    double volume(int axis, double refHeight = 0.0) const;
    
    // Return a normal vector with magnitude equal to total area.
    Vector3d arealNormal() const;
    Matrix3d nnTA() const;
    
    // Return some sensitivity information
    std::map<unsigned int, Matrix3<Vector3d> > DNNT(int axis) const;
    std::map<unsigned int, Vector3d> DV() const;
    
    // carry out interesting integrals on the ClippedTriangle.
    void intDz(std::map<unsigned int, Vector3d> & outIntegrals) const;
    void intZDh(std::map<unsigned int, Vector3d> & outIntegrals, int upAxis,
        double refHeight = 0.0) const;
    
    // These are the two parts of the sensitivity of NNT.
    // intNNTDh is the integral around the triangle contour of NNT times the
    //  perturbation of the boundary.
    // intDnnT is the integral within the triangle of the perturbation of nnT.
    void intNNTDh(std::map<unsigned int, Matrix3<Vector3d> > & outIntegrals) const;
    void intDnnT(std::map<unsigned int, Matrix3<Vector3d> > & outIntegrals) const;
    
    SimpleMesh::SensitiveEdge & storedEdge(int ee);
    
    void checkForInconsistentEdges() const;
    
private:
    const SimpleMesh::Triangle* mTriangle;
    
    std::vector<Vector3d> mVertices;
    std::vector<const SimpleMesh::SensitiveEdge*> mEdges;
    
    SimpleMesh::SensitiveEdge mEdgeStorage[6];
};

std::ostream & operator<<(std::ostream & str, const ClippedTriangle & ct);

void printJacobian(const std::map<unsigned int, Vector3d> & dFdv);
void printJacobian(const std::map<unsigned int, Matrix3<Vector3d> > & dFdv);















#endif
