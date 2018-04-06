//
//  XML.cpp
//  Trogdor6
//
//  Created by Paul Hansen on 3/4/18.
//  Copyright 2018 Stanford University.
//
// *  This file is covered by the MIT license.  See LICENSE.txt.

#include "XML.h"

#include <string>
#include <sstream>
#include <stdexcept>

#include "XMLExtras.h"
#include "tinyxml.h"

namespace SimpleMesh
{


void loadControlVertices(std::vector<SimpleMesh::ControlVertex> & controlVertices,
    const TiXmlElement* verticesXML)
{
    
    const TiXmlElement* vertexXML = verticesXML->FirstChildElement("Vertex");
    
    while (vertexXML)
    {
        SimpleMesh::ControlVertex controlVertex;
        
        SimpleMesh::loadControlVertex(controlVertex, vertexXML, controlVertices.size());
        controlVertices.push_back(controlVertex);
        
        vertexXML = vertexXML->NextSiblingElement("Vertex");
    }
}

void loadControlVertex(SimpleMesh::ControlVertex & controlVertex, const TiXmlElement* vertexXML, int newIndex)
{
    Vector3d point;
    Vector3b freeDirections;
    
    sGetMandatoryAttribute(vertexXML, "position", point);
    sGetOptionalAttribute(vertexXML, "freeDirections", freeDirections, Vector3b(0,0,0));
    
    controlVertex = SimpleMesh::ControlVertex(point, freeDirections, newIndex);
    
}


void loadMeshes(std::vector<std::vector<SimpleMesh::Triangle> > & allMeshes,
    const TiXmlElement* meshXML,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices)
{
    const TiXmlElement* solidXML = meshXML->FirstChildElement("Solid");
    
    while (solidXML)
    {
        std::vector<SimpleMesh::Triangle> meshFaces;
        SimpleMesh::loadMesh(meshFaces, solidXML, controlVertices);
        allMeshes.push_back(meshFaces);
        
        solidXML = solidXML->NextSiblingElement("Solid");
    }
}

void loadMesh(std::vector<SimpleMesh::Triangle> & meshFaces,
    const TiXmlElement* solidXML,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices)
{
    const TiXmlElement* faceXML = solidXML->FirstChildElement("Face");
    while (faceXML)
    {
        SimpleMesh::Triangle meshTri;
        loadFace(meshTri, faceXML, controlVertices);
        meshFaces.push_back(meshTri);
        
        faceXML = faceXML->NextSiblingElement("Face");
    }
}

void loadFace(SimpleMesh::Triangle & meshTri,
    const TiXmlElement* faceXML,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices)
{
    Vector3i idxFaceVertices;
    Vector3i idxFaceControlVertices;
    sGetMandatoryAttribute(faceXML, "vertices", idxFaceVertices);
    sGetOptionalAttribute(faceXML, "controlVertices", idxFaceControlVertices, idxFaceVertices);
    
    std::vector<std::vector<int> > edgeControlVertices;
    loadEdgeControlVertices(edgeControlVertices, faceXML, controlVertices);
    
    Triangle3d tri(controlVertices.at(idxFaceVertices[0]).point(),
        controlVertices.at(idxFaceVertices[1]).point(),
        controlVertices.at(idxFaceVertices[2]).point());
    meshTri.triangle(tri);
    
    meshTri.controlVertices(controlVertices.at(idxFaceControlVertices[0]),
        controlVertices.at(idxFaceControlVertices[1]),
        controlVertices.at(idxFaceControlVertices[2]));
    
    if (edgeControlVertices.empty())
    {
        meshTri.edgeControlVertices(controlVertices.at(idxFaceControlVertices[0]),
            controlVertices.at(idxFaceControlVertices[1]),
            controlVertices.at(idxFaceControlVertices[2]));
    }
    else
    {
        for (int ee = 0; ee < 3; ee++)
        {
            std::vector<SimpleMesh::ControlVertex> tmpControlVertices;
            for (int vv = 0; vv < edgeControlVertices[ee].size(); vv++)
            {
                tmpControlVertices.push_back(controlVertices.at(edgeControlVertices[ee][vv]));
            }
            meshTri.edgeControlVertices(ee, tmpControlVertices);
        }
    }
    meshTri.cacheSensitivity();
}

void loadEdgeControlVertices(std::vector<std::vector<int> > & edgeControlVertices,
    const TiXmlElement* faceXML,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices)
{
    const TiXmlElement* edgeXML = faceXML->FirstChildElement("Edge");
    while (edgeXML)
    {
        std::string cvString = edgeXML->Attribute("controlVertices");
        std::istringstream str(cvString);
        int idxControlVert;
        std::vector<int> edgeCVs;
        while (str >> idxControlVert)
        {
            edgeCVs.push_back(idxControlVert);
        }
        
        if (edgeCVs.size() != 2 && edgeCVs.size() != 6)
        {
            std::ostringstream err;
            err << "Edge must have 2 or 6 control vertices (possible location: row "
                << edgeXML->Row() << ")";
            throw std::runtime_error(err.str().c_str());
        }
        edgeControlVertices.push_back(edgeCVs);
        edgeXML = edgeXML->NextSiblingElement("Edge");
    }
    if (edgeControlVertices.size() != 0 && edgeControlVertices.size() != 3)
    {
        std::ostringstream err;
        err << "Face must have 0 or 3 Edge elements (possible location: row "
            << faceXML->Row() << ")";
        throw std::runtime_error(err.str().c_str());
    }
}



void createMeshesXML(const std::vector<std::vector<SimpleMesh::Triangle> > & meshes,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices,
    TiXmlElement & outMeshesXML)
{
    std::vector<SimpleMesh::ControlVertex> allVerts(controlVertices);
    // Build vertex ID map.
    // This is necessary because I haven't got a true face-vertex data structure.
    // First put in the control vertices.
    std::map<Vector3d, int> vertexIds;
    for (int vv = 0; vv < controlVertices.size(); vv++)
    {
        vertexIds[controlVertices[vv].point()] = controlVertices[vv].id();
    }
    
    // Then put in the other (non-control) vertices.
    int nextVertId = vertexIds.size();
    for (int mm = 0; mm < meshes.size(); mm++)
    {
        const std::vector<SimpleMesh::Triangle> & tris = meshes.at(mm);
        for (int tt = 0; tt < tris.size(); tt++)
        {
            const Triangle3d & tri = tris[tt].triangle();
            for (int ii = 0; ii < 3; ii++)
            {
                if (vertexIds.count(tri[ii]) == 0)
                {
                    vertexIds[tri[ii]] = nextVertId;
                    allVerts.push_back(SimpleMesh::ControlVertex(tri[ii], Vector3b(false, false, false), nextVertId));
                    nextVertId++;
                }
            }
        }
    }
    
    TiXmlElement verticesXML("Vertices");
    createControlVerticesXML(allVerts, verticesXML);
    outMeshesXML.InsertEndChild(verticesXML);
    
    for (int mm = 0; mm < meshes.size(); mm++)
    {
        const std::vector<SimpleMesh::Triangle> & meshTriangles = meshes.at(mm);
        TiXmlElement meshXML("Solid");
        meshXML.SetAttribute("id", std::to_string(mm).c_str());
        
        createMeshXML(meshTriangles, vertexIds, meshXML);
        outMeshesXML.InsertEndChild(meshXML);
    }
}

void createControlVerticesXML(const std::vector<SimpleMesh::ControlVertex> & controlVertices,
    TiXmlElement & outVerticesXML)
{
    for (int vv = 0; vv < controlVertices.size(); vv++)
    {
        const SimpleMesh::ControlVertex & cv = controlVertices.at(vv);
        
        TiXmlElement vertexXML("Vertex");
        
        std::ostringstream pointStream;
        pointStream.precision(20);
        pointStream << cv.point()[0] << " " << cv.point()[1] << " " << cv.point()[2];
        vertexXML.SetAttribute("position", pointStream.str().c_str());
        
        std::ostringstream freeDirStream;
        freeDirStream.precision(20);
        freeDirStream << cv.freeDirections()[0] << " " << cv.freeDirections()[1] << " " << cv.freeDirections()[2];
        vertexXML.SetAttribute("freeDirections", freeDirStream.str().c_str());

        outVerticesXML.InsertEndChild(vertexXML);
    }
}

void createMeshXML(const std::vector<SimpleMesh::Triangle> & meshFaces,
    const std::map<Vector3d, int> & vertexIds,
    TiXmlElement & outMeshXML)
{
    for (int ff = 0; ff < meshFaces.size(); ff++)
    {
        const SimpleMesh::Triangle & tri = meshFaces.at(ff);
        TiXmlElement faceXML("Face");
        
        std::ostringstream verticesStream;
        verticesStream << vertexIds.find(tri.triangle()[0])->second << " "
            << vertexIds.find(tri.triangle()[1])->second << " "
            << vertexIds.find(tri.triangle()[2])->second;
        faceXML.SetAttribute("vertices", verticesStream.str().c_str());
        
        std::ostringstream cvStream;
        cvStream << tri.controlVertices()[0].id() << " "
            << tri.controlVertices()[1].id() << " "
            << tri.controlVertices()[2].id();
        faceXML.SetAttribute("controlVertices", cvStream.str().c_str());
        
        outMeshXML.InsertEndChild(faceXML);
    }
}
    
}; // namespace SimpleMesh
