//
//  main.cpp
//  utility
//
//  Created by Paul Hansen on 3/4/18.
//  Copyright 2018 Stanford University.
// *  This file is covered by the MIT license.  See LICENSE.txt.
//

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "tinyxml.h"
#include "XMLExtras.h"
#include "XML.h"
#include "SimpleMesh.h"
#include "Search.h"

using namespace std;

int main(int argc, const char** argv)
{
    string errorMessage, versionString;
    
    string filename("mesh.xml");
    
    TiXmlBase::SetCondenseWhiteSpace(0); // for ascii material layers with \n
    
    TiXmlDocument* thisDoc = new TiXmlDocument(filename.c_str());
    if (!(thisDoc->LoadFile()))
    {
        ostringstream err;
        err << "Could not load parameter file " << filename << endl;
        err << "Reason: " << thisDoc->ErrorDesc() << endl;
        err << "(possible location: row " << thisDoc->ErrorRow() << " column "
             << thisDoc->ErrorCol() << " .)" << endl;
        throw(std::logic_error(err.str()));
    }
    
    TiXmlElement* meshXML = thisDoc->RootElement();
    
    // ======= Load vertices
    
    std::vector<SimpleMesh::ControlVertex> controlVertices;
    TiXmlElement* verticesXML = meshXML->FirstChildElement("Vertices");
    SimpleMesh::loadControlVertices(controlVertices, verticesXML);
    
    // ======= Load solids
    
    std::vector<std::vector<SimpleMesh::Triangle> > allMeshes;
    SimpleMesh::loadMeshes(allMeshes, meshXML, controlVertices);
    
    for (int mm = 0; mm < allMeshes.size(); mm++)
    {
        SimpleMesh::determineNeighbors(allMeshes[mm]);
        SimpleMesh::determineEdgeControlVertices(allMeshes[mm]);
//        const std::vector<SimpleMesh::Triangle> & aMesh = allMeshes.at(mm);
//        for (int tt = 0; tt < aMesh.size(); tt++)
//        {
//            std::cout << aMesh[tt].neighbors() << "\n";
//        }
    }
    
    // ======= Write meshes back to XML
    
    TiXmlElement writeMe("Mesh");
    SimpleMesh::createMeshesXML(allMeshes, controlVertices, writeMe);
    
    TiXmlDocument writeDoc("out.xml");
    writeDoc.InsertEndChild(writeMe);
    writeDoc.SaveFile();
    
    return 0;
}



