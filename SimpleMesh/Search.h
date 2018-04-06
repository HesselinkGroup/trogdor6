/*
 *  Search.h
 *  SimpleMesh
 *
 *  Created by Paul C Hansen on 6/27/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef SEARCH_H
#define SEARCH_H

#include "SimpleMesh.h"
#include <vector>
#include <set>

namespace SimpleMesh
{

class AllControlVerticesIn
{
public:
    AllControlVerticesIn() {}
    AllControlVerticesIn(int cv0, int cv1,
        int cv2 = -1, int cv3 = -1, int cv4 = -1,
        int cv5 = -1, int cv6 = -1, int cv7 = -1);
    
    bool operator()(const std::vector<ControlVertex> & line) const;
    bool operator()(const Triangle & triangle) const;
    
private:
    std::set<int> mControlVertices;
};

// Find lines or triangles whose control vertices are in the given set of
// control vertices.

std::vector<SensitiveEdge>
linesWithControlVertices(const Triangle & triangle,
    int cv0, int cv1, int cv2 = -1, int cv3 = -1,
    int cv4 = -1, int cv5 = -1, int cv6 = -1, int cv7 = -1);

std::vector<SensitiveEdge>
linesWithControlVertices(const std::vector<Triangle> & triangles,
    int cv0, int cv1, int cv2 = -1, int cv3 = -1,
    int cv4 = -1, int cv5 = -1, int cv6 = -1, int cv7 = -1);

std::vector<Triangle>
trianglesWithControlVertices(const std::vector<Triangle> & triangles,
    int cv0, int cv1, int cv2 = -1, int cv3 = -1);

// Find lines or triangles within the given bounding box.

std::vector<SensitiveEdge>
linesInBounds(const std::vector<Triangle> & triangles, Rect3d bounds);

std::vector<Triangle>
trianglesInBounds(const std::vector<Triangle> & triangles, Rect3d bounds);

// Set neighbor indices for a mesh.
void determineNeighbors(std::vector<Triangle> & triangles);
void determineEdgeControlVertices(std::vector<Triangle> & triangles);


// Determine neighbors, determine edge CVs, cache sensitivities
void initMeshSensitivity(std::vector<Triangle> & triangles);

void initMeshSensitivity(std::vector<std::vector<Triangle> > & meshes);

} // namespace SimpleMesh

#endif
