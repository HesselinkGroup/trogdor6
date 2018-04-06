//
//  XML.hpp
//  Trogdor6
//
//  Created by Paul Hansen on 3/4/18.
//  Copyright 2018 Stanford University.
//
// *  This file is covered by the MIT license.  See LICENSE.txt.

#ifndef XML_hpp
#define XML_hpp

#include "SimpleMesh.h"
#include <vector>
#include <map>

class TiXmlElement;

namespace SimpleMesh
{

void loadControlVertices(std::vector<SimpleMesh::ControlVertex> & controlVertices,
    const TiXmlElement* verticesXML);

void loadControlVertex(SimpleMesh::ControlVertex & controlVertex,
    const TiXmlElement* vertexXML,
    int newIndex);

void loadMeshes(std::vector<std::vector<SimpleMesh::Triangle> > & allMeshes,
    const TiXmlElement* meshXML,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices);

void loadMesh(std::vector<SimpleMesh::Triangle> & meshFaces,
    const TiXmlElement* solidXML,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices);

void loadFace(SimpleMesh::Triangle & meshTri,
    const TiXmlElement* faceXML,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices);

void loadEdgeControlVertices(std::vector<std::vector<int> > & edgeControlVertices,
    const TiXmlElement* faceXML,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices);
    
    
    
    
    
void createMeshesXML(const std::vector<std::vector<SimpleMesh::Triangle> > & meshes,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices,
    TiXmlElement & outMeshesXML);

void createControlVerticesXML(const std::vector<SimpleMesh::ControlVertex> & controlVertices,
    TiXmlElement & outVerticesXML);
    
    
void createMeshXML(const std::vector<SimpleMesh::Triangle> & meshFaces,
    const std::map<Vector3d, int> & vertexIds,
    TiXmlElement & outMeshXML);

}; // namespace SimpleMesh
#endif /* XML_hpp */
