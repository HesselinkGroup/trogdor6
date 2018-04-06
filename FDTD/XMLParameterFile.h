/*
 *  XMLParameterFile.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 1/31/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _XMLPARAMETERFILE_
#define _XMLPARAMETERFILE_

#include "SimulationDescription.h"
#include "tinyxml.h"

#include <string>
#include <vector>
#include <stdexcept>

class XMLParameterFile
{
	friend class SimulationDescription;
public:
	XMLParameterFile(const std::string & filename) throw(std::logic_error); 
private:
	void load(SimulationDescription & sim) const throw(std::logic_error);
	
	std::vector<GridDescPtr> loadGrids(const TiXmlElement* parent,
        const SimulationDescription & sim) const;
	std::set<std::string> collectGridNames(const TiXmlElement* parent) const;
	std::set<std::string> collectMaterialNames(const TiXmlElement* parent)
		const;
	std::vector<MaterialDescPtr> loadMaterials(const TiXmlElement* parent)
		const;
	std::vector<OutputDescPtr> loadOutputs(const TiXmlElement* parent) const;
    OutputDescPtr loadAFieldOutput(const TiXmlElement* elem) const;
    
    std::vector<GridReportDescPtr> loadGridReports(const TiXmlElement* parent)
        const;
    GridReportDescPtr loadAGridReport(const TiXmlElement* elem) const;
    
    std::vector<SourceDescPtr> loadSources(const TiXmlElement* parent) const;
    SourceDescPtr loadAFieldSource(const TiXmlElement* elem) const;
    
    std::vector<CurrentSourceDescPtr> loadCurrentSources(
        const TiXmlElement* parent) const;
    CurrentSourceDescPtr loadACurrentSource(const TiXmlElement* elem) const;
    
	std::vector<HuygensSurfaceDescPtr> loadTFSFSources(
		const TiXmlElement* parent, const std::set<std::string> & allGridNames)
		const;
	std::vector<HuygensSurfaceDescPtr> loadCustomSources(
		const TiXmlElement* parent) const;
	std::vector<HuygensSurfaceDescPtr> loadLinks(
		const TiXmlElement* parent, const std::set<std::string> & allGridNames)
		const;
	AssemblyDescPtr loadAssembly(const TiXmlElement* parent,
		const std::set<std::string> & allGridNames,
		std::map<std::string, MaterialDescPtr> materials) const;
	
    InstructionPtr loadANewMesh(const TiXmlElement* elem,
        std::map<std::string, MaterialDescPtr> materials) const;
    InstructionPtr loadAMesh(const TiXmlElement* elem,
        std::map<std::string, MaterialDescPtr> materials,
        int & inOutControlVertexId) const;
	InstructionPtr loadACopyFrom(const TiXmlElement* elem,
		const std::set<std::string> & allGridNames) const;
    InstructionPtr loadAnExtrude(const TiXmlElement* elem) const;
    InstructionPtr loadABackground(const TiXmlElement* elem,
        std::map<std::string, MaterialDescPtr> materials) const;
    InstructionPtr loadAPerm_ity(const TiXmlElement* elem,
        int numMeshes,
        const std::map<std::string, MaterialDescPtr> & materials) const;
	
	Duration loadADuration(const TiXmlElement* elem) const;
    Region loadARegion(const TiXmlElement* elem) const;
    Vector3i loadAnOmitSide(const TiXmlElement* elem) const;
    
	Pointer<TiXmlDocument> mDocument;
};



#endif
