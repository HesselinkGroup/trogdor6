/*
 *  XMLParameterFile.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 1/31/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "XMLParameterFile.h"

#include <sstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include <set>
#include <climits>

#include "SimpleMesh/XML.h"
#include "SimpleMesh/SimpleMesh.h"
#include "YeeUtilities.h"
#include "StreamFromString.h"
#include "XMLExtras.h"

using namespace std;
using namespace YeeUtilities;

static Vector3i sUnitVectorFromString(const string & axisString)
	throw(std::logic_error);

#pragma mark *** XMLParameterFile Implementation ***

XMLParameterFile::
XMLParameterFile(const string & filename) throw(std::logic_error)
{
	string errorMessage, versionString;
	
    TiXmlBase::SetCondenseWhiteSpace(0); // for ascii material layers with \n
	
    Pointer<TiXmlDocument> thisDoc(new TiXmlDocument(filename.c_str()));
	if (!(thisDoc->LoadFile()))
    {
		ostringstream err;
        err << "Could not load parameter file " << filename << endl;
        err << "Reason: " << thisDoc->ErrorDesc() << endl;
        err << "(possible location: row " << thisDoc->ErrorRow() << " column "
             << thisDoc->ErrorCol() << " .)" << endl;
        throw(std::logic_error(err.str()));
    }
    if (!sTryGetAttribute(thisDoc->RootElement(), "version", versionString))
    {
        ostringstream err;
        err << "Could not load parameter file " << filename << endl;
        err << "Reason: no version string." << endl;
        throw(std::logic_error(err.str()));
    }
    else
        mDocument = thisDoc;
}

// It is the job of XMLParameterFile to validate the XML and to verify that 
// all the necessary elements and attributes are present and properly-formatted.
// It is the job of the SimulationDescription etc. classes to verify that the
// information in the XML file makes sense; this is handled in the constructors.

void XMLParameterFile::
load(SimulationDescription & sim) const throw(std::logic_error)
{
	const TiXmlElement* simulationXML = mDocument->RootElement();
	string versionString;
    ostringstream err;
    
	if (!simulationXML)
		throw std::logic_error("XML file has no root element.");
    
    if (!sTryGetAttribute(simulationXML, "version", versionString))
        throw std::logic_error("Simulation is unversioned!");
    if (versionString != "6.0")
    {
        err << "Cannot read version " << versionString << " file.";
        throw(std::logic_error(sErr(err.str(), simulationXML)));
    }
	
	Precision::Float dx, dy, dz, dt;
	int numTimesteps;
		
	sGetMandatoryAttribute(simulationXML, "dx", dx);
	sGetMandatoryAttribute(simulationXML, "dy", dy);
	sGetMandatoryAttribute(simulationXML, "dz", dz);
	sGetMandatoryAttribute(simulationXML, "dt", dt);
	sGetMandatoryAttribute(simulationXML, "numT", numTimesteps);
 
    std::string electromagneticMode;
    sGetOptionalAttribute(simulationXML, "electromagneticMode", electromagneticMode, std::string("3d"));
    if (electromagneticMode == "3d")
    {
        sim.setElectromagneticMode(k3d);
    }
    else if (electromagneticMode == "te2d")
    {
        sim.setElectromagneticMode(kTE2d);
    }
    else if (electromagneticMode == "tm2d")
    {
        sim.setElectromagneticMode(kTM2d);
    }
    else if (electromagneticMode == "2d")
    {
        sim.setElectromagneticMode(k2d);
    }
    else if (electromagneticMode == "1d")
    {
        sim.setElectromagneticMode(k1d);
    }
    else
    {
        throw std::runtime_error("Invalid electromagnetic mode.");
    }
    
	
	// Only the calls that may throw exceptions are in the try-catch block.
	try {
		sim.setDiscretization(Precision::Vec3(dx,dy,dz), dt);
		sim.setNumTimesteps(numTimesteps);
	} catch (std::logic_error & e) {
		throw(std::logic_error(sErr(e.what(), simulationXML)));
	}
	// (i.e. these three calls don't throw my own exceptions.)
    // I am really way too loose about exceptions in general.
	sim.setMaterials(loadMaterials(simulationXML));
	sim.setGrids(loadGrids(simulationXML, sim));
	sim.setGridPointers();
}

vector<GridDescPtr> XMLParameterFile::
loadGrids(const TiXmlElement* parent, const SimulationDescription & sim) const
{
	assert(parent);
	vector<GridDescPtr> gridVector;
	const TiXmlElement* elem = parent->FirstChildElement("Grid");
	
	set<string> allGridNames = collectGridNames(parent);
	set<string> allMaterialNames = collectMaterialNames(parent);
	
    while (elem) // FOR EACH GRID: Load data
	{
		GridDescPtr gridDesc;
        Rect3i yeeCells;
        Rect3i halfCells;
        Vector3i numYeeCells;
        Vector3i numHalfCells;
		Rect3i nonPMLYeeCells;
        Rect3i nonPMLHalfCells;
        Precision::Vec3 origin;
        Map<string, string> pmlParams;
		string name;
		
		sGetMandatoryAttribute(elem, "name", name);
        
        sGetOptionalAttribute(elem, "origin", origin, Precision::Vec3(0,0,0));
        
        sGetMandatoryAttribute(elem, "nonPML", nonPMLYeeCells);
        nonPMLHalfCells = yeeToHalf(nonPMLYeeCells);
        
        sGetMandatoryAttribute(elem, "yeeCells", yeeCells);
        
        halfCells = yeeToHalf(yeeCells);
        
        const TiXmlElement* pmlParamXML = elem->FirstChildElement("PML");
        while (pmlParamXML != 0L)
        {
            if (false == pmlParams.empty())
                throw(std::logic_error(sErr(
                    "More than one PML element detected", elem)));
            pmlParams = sGetAttributes(pmlParamXML);
            
            pmlParamXML = pmlParamXML->NextSiblingElement("PML");
        }
        
		// We need to wrap the constructor errors from GridDescription() to
		// put XML row/column information in.
		try {
			gridDesc = GridDescPtr(new GridDescription(name,
                sim.electromagneticMode(),
				yeeCells,
				nonPMLHalfCells, sim.dxyz(),
                sim.dt(),
                sim.numTimesteps()));
		} catch (std::logic_error & e) {
			throw(std::logic_error(sErr(e.what(), elem)));
		}
        
        gridDesc->setOrigin(origin);
		
		// However, these function calls all take care of their own XML file
		// row/column information so we don't need to re-throw their
		// exceptions.
		gridDesc->setOutputs(loadOutputs(elem));
        gridDesc->setGridReports(loadGridReports(elem));
		gridDesc->setSources(loadSources(elem));
        gridDesc->setCurrentSources(loadCurrentSources(elem));
		gridDesc->setAssembly(loadAssembly(elem, allGridNames, sim.materials()));
        try {
            gridDesc->setPMLParams(pmlParams);
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
		
		// half the huygens surfaces come from links, and half from TFSF
		// sources.
		vector<HuygensSurfaceDescPtr> huygensSurfaces =
			loadTFSFSources(elem, allGridNames);
		vector<HuygensSurfaceDescPtr> huygensLinks =
			loadLinks(elem, allGridNames);
		vector<HuygensSurfaceDescPtr> customSources =
			loadCustomSources(elem);
		huygensSurfaces.insert(huygensSurfaces.begin(), huygensLinks.begin(),
			huygensLinks.end());
		huygensSurfaces.insert(huygensSurfaces.begin(), customSources.begin(),
			customSources.end());
		
		gridDesc->setHuygensSurfaces(huygensSurfaces);
		
		gridVector.push_back(gridDesc);
        elem = elem->NextSiblingElement("Grid");
	}
	
	return gridVector;
}

set<string> XMLParameterFile::
collectGridNames(const TiXmlElement* parent) const
{
	assert(parent);
	set<string> names;
	
	const TiXmlElement* gridElem = parent->FirstChildElement("Grid");
	
	while (gridElem)
	{
		string gridName;
		sGetMandatoryAttribute(gridElem, "name", gridName);
		
		if (names.count(gridName) > 0)
		{
			throw(std::logic_error(sErr("Multiple grids with the same name!",
				gridElem)));
		}
		
		names.insert(gridName);
		gridElem = gridElem->NextSiblingElement("Grid");
	}
	return names;
}

set<string> XMLParameterFile::
collectMaterialNames(const TiXmlElement* parent) const
{
	assert(parent);
	set<string> names;
	
    const TiXmlElement* matElem = parent->FirstChildElement("Material");
    
    while (matElem)
    {
        string matName;
        sGetMandatoryAttribute(matElem, "name", matName);
        names.insert(matName);
        matElem = matElem->NextSiblingElement("Material");
    }
    
	return names;
}

vector<MaterialDescPtr> XMLParameterFile::
loadMaterials(const TiXmlElement* parent) const
{
	assert(parent);
	vector<MaterialDescPtr> materials;
    
    const TiXmlElement* elem = parent->FirstChildElement("Material");
    while (elem)
    {
        string name;
        string numeratorString, denominatorString;
        vector<double> numerator;
        vector<double> denominator;
        double coefficient;
        
        sGetMandatoryAttribute(elem, "name", name);
        sGetMandatoryAttribute(elem, "zNumerator", numeratorString);
        sGetMandatoryAttribute(elem, "zDenominator", denominatorString);
        
        // FIXME: No error checking on the coefficients
        istringstream str(numeratorString);
        while (str >> coefficient)
            numerator.push_back(coefficient);
        istringstream str2(denominatorString);
        while (str2 >> coefficient)
            denominator.push_back(coefficient);
        
        Precision::Polynomial numPoly, denomPoly;
        for (int nn = 0; nn < numerator.size(); nn++)
            numPoly.coefficient(nn, numerator[nn]);
        for (int nn = 0; nn < denominator.size(); nn++)
            denomPoly.coefficient(nn, denominator[nn]);
        
        // The material ID will be the same as its index.
        try {
            MaterialDescPtr material(new MaterialDescription(
                materials.size(), name,
                Precision::RationalFunction(numPoly, denomPoly)));
            materials.push_back(material);
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
        
        elem = elem->NextSiblingElement("Material");
    }
    
    for (int nn = 0; nn < materials.size(); nn++)
        assert(materials[nn]->id() == nn);
	
	return materials;	
}


vector<OutputDescPtr> XMLParameterFile::
loadOutputs(const TiXmlElement* parent) const
{
	assert(parent);
	vector<OutputDescPtr> outputs;
	
	const TiXmlElement* elem = parent->FirstChildElement("FieldOutput");
    while (elem)
	{
        outputs.push_back(loadAFieldOutput(elem));
        elem = elem->NextSiblingElement("FieldOutput");
    }
	return outputs;
}

OutputDescPtr XMLParameterFile::
loadAFieldOutput(const TiXmlElement* elem) const
{
    string fields;
    string file;
    Precision::Vec3 interpolationPoint;
    bool isInterpolated;
    sGetMandatoryAttribute(elem, "fields", fields);
    sGetMandatoryAttribute(elem, "file", file);
    isInterpolated = sTryGetAttribute(elem, "interpolate", interpolationPoint);
    
    vector<Region> regions;
    vector<Duration> durations;
    
    const TiXmlElement* child = elem->FirstChildElement("Duration");
    while (child != 0L)
    {
        durations.push_back(loadADuration(child));
        child = child->NextSiblingElement("Duration");
    }
    child = elem->FirstChildElement("Region");
    while (child != 0L)
    {
        regions.push_back(loadARegion(child));
        child = child->NextSiblingElement("Region");
    }
    
    if (regions.size() == 0)
        regions.push_back(Region());
    if (durations.size() == 0)
        durations.push_back(Duration());
    
    if (isInterpolated)
    {
        throw(std::logic_error(sErr("Interpolation not supported yet.", elem)));
//        try {
//            OutputDescPtr f(new OutputDescription(
//                fields, file, interpolationPoint, regions, durations));
//            return f;
//        } catch (std::logic_error & e) {
//            throw(std::logic_error(sErr(e.what(), elem)));
//        }
    }
    else
    {
        try {
            OutputDescPtr f(new OutputDescription(
                fields, file, regions, durations));
            return f;
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    
    LOG << "Shouldn't have gotten here.\n";
    return OutputDescPtr(0L);
}


vector<GridReportDescPtr> XMLParameterFile::
loadGridReports(const TiXmlElement* parent) const
{
	assert(parent);
	vector<GridReportDescPtr> reports;
	
	const TiXmlElement* elem = parent->FirstChildElement("GridReport");
    while (elem)
	{
        reports.push_back(loadAGridReport(elem));
        elem = elem->NextSiblingElement("GridReport");
    }
	return reports;
}

GridReportDescPtr XMLParameterFile::
loadAGridReport(const TiXmlElement* elem) const
{
    assert(elem);
    const TiXmlElement* child;
    string fields;
    string fileName;
    vector<Region> regions;
    
    sGetMandatoryAttribute(elem, "file", fileName);
    
    child = elem->FirstChildElement("Region");
    while (child != 0L)
    {
        regions.push_back(loadARegion(child));
        child = child->NextSiblingElement("Region");
    }
    
    GridReportDescPtr reportPtr;
    if (sTryGetAttribute(elem, "fields", fields))
    {
        cerr << "Warning: loading GridReport with specified fields, but for "
            "now the fields will be ignored and every octant will be present in"
            " the output file.\n";
        try {
            reportPtr = GridReportDescPtr(new GridReportDescription(
                fields, fileName, regions));
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    else
    {
        try {
            reportPtr = GridReportDescPtr(new GridReportDescription(
                fileName, regions));
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    return reportPtr;
}

vector<SourceDescPtr> XMLParameterFile::
loadSources(const TiXmlElement* parent) const
{
    assert(parent);
    vector<SourceDescPtr> sources;
    
    const TiXmlElement* child = parent->FirstChildElement("HardSource");
    while (child != 0L)
    {
        sources.push_back(loadAFieldSource(child));
        child = child->NextSiblingElement("HardSource");
    }
    
    child = parent->FirstChildElement("AdditiveSource");
    while (child != 0L)
    {
        sources.push_back(loadAFieldSource(child));
        child = child->NextSiblingElement("AdditiveSource");
    }
    return sources;
}

SourceDescPtr XMLParameterFile::
loadAFieldSource(const TiXmlElement* elem) const
{
    assert(elem);
    string fields;
    string formula;
    string timeFile;
    string spaceFile;
    string spaceTimeFile;
    Precision::Vec3 polarization;
    SourceFields srcFields;
    vector<Duration> durations;
    vector<Region> regions;
    bool isSoft = (string("AdditiveSource") == elem->Value());
    
    sGetMandatoryAttribute(elem, "fields", fields);
    if (sTryGetAttribute(elem, "polarization", polarization))
    {
        try {
            srcFields = SourceFields(fields, polarization);
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    else
    {
        try {
            srcFields = SourceFields(fields);
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    
    const TiXmlElement* child = elem->FirstChildElement("Duration");
    while (child != 0L)
    {
        durations.push_back(loadADuration(child));
        child = child->NextSiblingElement("Duration");
    }
    if (durations.size() == 0)
        durations.push_back(Duration());
    
    child = elem->FirstChildElement("Region");
    while (child != 0L)
    {
        regions.push_back(loadARegion(child));
        child = child->NextSiblingElement("Region");
    }
    
    SourceDescPtr src;
    if (sTryGetAttribute(elem, "timeFile", timeFile))
    {
        try {
            src = SourceDescPtr(SourceDescription::
                newTimeSource(timeFile, srcFields, isSoft, regions, durations));
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    else if (sTryGetAttribute(elem, "spaceTimeFile", spaceTimeFile))
    {
        try {
            src = SourceDescPtr(SourceDescription::
                newSpaceTimeSource(spaceTimeFile, srcFields, isSoft, regions,
                    durations));
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    else if (sTryGetAttribute(elem, "formula", formula))
    {
        try {
            src = SourceDescPtr(SourceDescription::
                newFormulaSource(formula, srcFields, isSoft, regions, durations));
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    else
        throw(std::logic_error(sErr("HardSource needs timeFile, spaceTimeFile or "
            "formula attribute.", elem)));
    
    return src;
}

vector<CurrentSourceDescPtr> XMLParameterFile::
loadCurrentSources(const TiXmlElement* parent) const
{    assert(parent);

    vector<CurrentSourceDescPtr> sources;
    
    const TiXmlElement* child = parent->FirstChildElement("CurrentSource");
    while (child != 0L)
    {
        sources.push_back(loadACurrentSource(child));
        child = child->NextSiblingElement("CurrentSource");
    }
    
    return sources;
}


CurrentSourceDescPtr XMLParameterFile::
loadACurrentSource(const TiXmlElement* elem) const
{
    assert(elem);
    string fields;
    string formula;
    string timeFile;
    string spaceFile;
    string spaceTimeFile;
    Precision::Vec3 polarization;
    SourceCurrents srcCurrents;
    vector<Duration> durations;
    vector<Region> regions;
    
    sGetMandatoryAttribute(elem, "fields", fields);
    if (sTryGetAttribute(elem, "polarization", polarization))
    {
        try {
            srcCurrents = SourceCurrents(fields, polarization);
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    else
    {
        try {
            srcCurrents = SourceCurrents(fields);
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    
    const TiXmlElement* child = elem->FirstChildElement("Duration");
    while (child != 0L)
    {
        durations.push_back(loadADuration(child));
        child = child->NextSiblingElement("Duration");
    }
    if (durations.size() == 0)
        durations.push_back(Duration());
    
    child = elem->FirstChildElement("Region");
    while (child != 0L)
    {
        regions.push_back(loadARegion(child));
        child = child->NextSiblingElement("Region");
    }
    
    CurrentSourceDescPtr src;
    if (sTryGetAttribute(elem, "timeFile", timeFile))
    {
        try {
            src = CurrentSourceDescPtr(CurrentSourceDescription::
                newTimeSource(timeFile, srcCurrents, regions, durations));
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    else if (sTryGetAttribute(elem, "spaceTimeFile", spaceTimeFile))
    {
        try {
            src = CurrentSourceDescPtr(CurrentSourceDescription::
                newSpaceTimeSource(spaceTimeFile, srcCurrents, regions,
                    durations));
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    else if (sTryGetAttribute(elem, "formula", formula))
    {
        try {
            src = CurrentSourceDescPtr(CurrentSourceDescription::
                newFormulaSource(formula, srcCurrents, regions, durations));
        } catch (std::logic_error & e) {
            throw(std::logic_error(sErr(e.what(), elem)));
        }
    }
    else
        throw(std::logic_error(sErr("CurrentSource needs timeFile, spaceTimeFile or "
            "formula attribute.", elem)));
    
    return src;
}

vector<HuygensSurfaceDescPtr> XMLParameterFile::
loadTFSFSources(const TiXmlElement* parent, const set<string> & allGridNames)
	const
{
	assert(parent);
	vector<HuygensSurfaceDescPtr> tfsfSources;
	
	const TiXmlElement* elem = parent->FirstChildElement("TFSFSource");
    while (elem)
	{
		string field;
		Precision::Vec3 polarization;
        Rect3i yeeCells;
		Vector3i direction;
		string tfsfType;
		string formula;
		string file;
        SourceFields srcFields;
        Duration duration;
        
        const TiXmlElement *durationXML;
        durationXML = elem->FirstChildElement("Duration");
        if (durationXML)
            duration = loadADuration(durationXML);
        
        sGetMandatoryAttribute(elem, "yeeCells", yeeCells);
        sGetMandatoryAttribute(elem, "direction", direction);
		sGetOptionalAttribute(elem, "type", tfsfType, string("TF"));
        
        if (!sTryGetAttribute(elem, "formula", formula) &&
            !sTryGetAttribute(elem, "file", file))
            throw(std::logic_error(sErr("TFSFSource needs either a formula or file"
                " attribute.", elem)));
        
		sGetMandatoryAttribute(elem, "fields", field);
        if (sTryGetAttribute(elem, "polarization", polarization))
            srcFields = SourceFields(field, polarization);
        else
            srcFields = SourceFields(field);
        
        set<Vector3i> omittedSides;
        const TiXmlElement* omitElem = elem->FirstChildElement("OmitSide");
        while (omitElem)
        {
            omittedSides.insert(loadAnOmitSide(omitElem));
            omitElem = omitElem->NextSiblingElement("OmitSide");
        }
        
		try {
            if (file != "")
                tfsfSources.push_back(HuygensSurfaceDescPtr(
                    HuygensSurfaceDescription::newTFSFTimeSource(srcFields,
                        file, direction, yeeCells, omittedSides,
                        tfsfType == "TF")));
            else
                tfsfSources.push_back(HuygensSurfaceDescPtr(
                    HuygensSurfaceDescription::newTFSFFormulaSource(srcFields,
                        formula, direction, yeeCells, omittedSides,
                        tfsfType == "TF")));
		} catch (std::logic_error & e) {
			throw(std::logic_error(sErr(e.what(), elem)));
		}
		
		elem = elem->NextSiblingElement("TFSFSource");
	}
	return tfsfSources;
}

vector<HuygensSurfaceDescPtr> XMLParameterFile::
loadCustomSources(const TiXmlElement* parent) const
{
	assert(parent);
	vector<HuygensSurfaceDescPtr> customSources;
	
	const TiXmlElement* elem = parent->FirstChildElement("CustomTFSFSource");
    while (elem)
	{
        string file;
        Vector3i symmetries(0,0,0);
        string tfsfType;
        Rect3i yeeCells;
        set<Vector3i> omittedSides;
        const TiXmlElement* durationXML;
        Duration duration;
		
		sGetMandatoryAttribute(elem, "file", file);
        sGetMandatoryAttribute(elem, "symmetries", symmetries);
        sGetOptionalAttribute(elem, "tfsfType", tfsfType, string("TF"));
        sGetMandatoryAttribute(elem, "yeeCells", yeeCells);
        
        durationXML = elem->FirstChildElement("Duration");
        if (durationXML != 0L)
            duration = loadADuration(durationXML);
        
        const TiXmlElement* omitElem = elem->FirstChildElement("OmitSide");
        while (omitElem)
        {
            omittedSides.insert(loadAnOmitSide(omitElem));
            omitElem = omitElem->NextSiblingElement("OmitSide");
        }
        
        try {
            HuygensSurfaceDescPtr source(HuygensSurfaceDescription::
                newCustomTFSFSource(file, symmetries, yeeCells, duration,
                omittedSides, tfsfType=="TF"));
            customSources.push_back(source);
		} catch (std::logic_error & e) {
			throw(std::logic_error(sErr(e.what(), elem)));
		}
		
		elem = elem->NextSiblingElement("CustomTFSFSource");
	}
	return customSources;
}

vector<HuygensSurfaceDescPtr> XMLParameterFile::
loadLinks(const TiXmlElement* parent, const set<string> & allGridNames) const
{
	assert(parent);
	vector<HuygensSurfaceDescPtr> links;
	
	const TiXmlElement* elem = parent->FirstChildElement("Link");
    while (elem)
	{
		string linkTypeStr;
		string sourceGridName;
		Rect3i temp;
        Rect3i fromYeeCells;
        Rect3i toYeeCells;
        set<Vector3i> omittedSides;
        
		sGetOptionalAttribute(elem, "type", linkTypeStr, string("TF"));
		sGetMandatoryAttribute(elem, "sourceGrid", sourceGridName);
        sGetMandatoryAttribute(elem, "fromYeeCells", fromYeeCells);
        sGetMandatoryAttribute(elem, "toYeeCells", toYeeCells);
		
        const TiXmlElement* omitElem = elem->FirstChildElement("OmitSide");
        while (omitElem)
        {
            omittedSides.insert(loadAnOmitSide(omitElem));
            omitElem = omitElem->NextSiblingElement("OmitSide");
        }
		try {	
			HuygensSurfaceDescPtr link(HuygensSurfaceDescription::
                newLink(sourceGridName, fromYeeCells, toYeeCells,
                    omittedSides, linkTypeStr == "TF"));
			links.push_back(link);
		} catch (std::logic_error & e) {
			throw(std::logic_error(sErr(e.what(), elem)));
		}
		
		elem = elem->NextSiblingElement("Link");
	}
	
	return links;
}

AssemblyDescPtr XMLParameterFile::
loadAssembly(const TiXmlElement* parent, const set<string> & allGridNames,
	map<string, MaterialDescPtr> materials) const
{
	assert(parent);
	AssemblyDescPtr assembly;
	vector<InstructionPtr> recipe;
	
	const TiXmlElement* assemblyXML = parent->FirstChildElement("Assembly");
	
	if (assemblyXML->NextSiblingElement("Assembly") != 0L)
		throw(std::logic_error(sErr("Only one Assembly element allowed per Grid",
			assemblyXML->NextSiblingElement("Assembly"))));
	
    // ======= Load vertices
    
    std::vector<SimpleMesh::ControlVertex> controlVertices;
    SimpleMesh::loadControlVertices(controlVertices, assemblyXML->FirstChildElement("Vertices"));
    
    // ======= Load solids
    
    std::vector<std::vector<SimpleMesh::Triangle> > allMeshes;
    SimpleMesh::loadMeshes(allMeshes, assemblyXML, controlVertices);
    
    // ======= Load permittivities and permeabilities
	
    const TiXmlElement* instructionXML = assemblyXML->FirstChildElement();
    
    int nextControlVertexId = 0; // running count
    bool hasBackground = false;
	while (instructionXML)
	{
		string instructionType = instructionXML->Value();
		
        if (instructionType == "Mesh")
        {
            // the control vertex id will be bumped by loadAMesh
            recipe.push_back(loadAMesh(instructionXML, materials, nextControlVertexId));
        }
        else if (instructionType == "Permittivity" || instructionType == "Permeability")
        {
            recipe.push_back(loadAPerm_ity(instructionXML, allMeshes.size(), materials));
        }
//        else if (instructionType == "NewMesh")
//        {
//            recipe.push_back(loadANewMesh(instruction, materials));
//        }
		else if (instructionType == "CopyFrom")
			recipe.push_back(loadACopyFrom(instructionXML, allGridNames));
        else if (instructionType == "Background")
        {
            if (hasBackground == false)
            {
                recipe.push_back(loadABackground(instructionXML, materials));
                hasBackground = true;
            }
            else
            {
                throw(std::runtime_error(sErr("Only one Background element "
                    "allowed per Grid", instructionXML)));
            }
        }   
		
		instructionXML = instructionXML->NextSiblingElement();
	}
	
	assembly = AssemblyDescPtr(new AssemblyDescription(recipe, controlVertices, allMeshes));
	return assembly;
}


//InstructionPtr XMLParameterFile::
//loadABlock(const TiXmlElement* elem, map<string, MaterialDescPtr> materials) const
//{
//    assert(elem);
//    assert(elem->Value() == string("Block"));
//
//    InstructionPtr blockPtr;
//    Rect3d rect;
//    MaterialDescPtr permittivity, permeability;
//    string material;
//
//    if (sTryGetAttribute(elem, "permittivity", material))
//    {
//        if (materials.count(material) == 0)
//            throw(std::logic_error(sErr("Unknown permittivity for Block", elem)));
//        permittivity = materials[material];
//    }
//    if (sTryGetAttribute(elem, "permeability", material))
//    {
//        if (materials.count(material) == 0)
//            throw(std::logic_error(sErr("Unknown permeability for Block", elem)));
//        permeability = materials[material];
//    }
//
//    sGetMandatoryAttribute(elem, "yeeBounds", rect);
//    try {
//        blockPtr = InstructionPtr(
//            (Block*) new Block(rect, permittivity, permeability));
//    } catch (std::logic_error & e) {
//        throw(std::logic_error(sErr(e.what(), elem)));
//    }
//
//    return blockPtr;
//}


InstructionPtr XMLParameterFile::
loadANewMesh(const TiXmlElement* elem, map<string, MaterialDescPtr> materials) const
{
    assert(elem);
    assert(elem->Value() == string("NewMesh"));
    
    InstructionPtr meshPtr;
    Rect3d rect;
    MaterialDescPtr permittivity, permeability;
    string material;
    
    if (sTryGetAttribute(elem, "permittivity", material))
    {
        if (materials.count(material) == 0)
            throw(std::logic_error(sErr("Unknown permittivity for Block", elem)));
        permittivity = materials[material];
    }
    if (sTryGetAttribute(elem, "permeability", material))
    {
        if (materials.count(material) == 0)
            throw(std::logic_error(sErr("Unknown permeability for Block", elem)));
        permeability = materials[material];
    }
    
    vector<Vector3i> faces;
    vector<Vector3i> controlFaces;
    
    const TiXmlElement* faceXML = elem->FirstChildElement("Face");
    while (faceXML != 0L)
    {
        Vector3i face, controlFace;
        sGetMandatoryAttribute(faceXML, "vertices", face);
        sGetOptionalAttribute(faceXML, "controlVertices", controlFace, face);
        
        faces.push_back(face);
        controlFaces.push_back(controlFace);
        faceXML = faceXML->NextSiblingElement("Face");
    }
    
    try {
        meshPtr = InstructionPtr((Mesh*) new NewMesh(
            faces, controlFaces, permittivity, permeability));
    } catch (std::logic_error & e) {
        throw(std::logic_error(sErr(e.what(), elem)));
    }
    
    return meshPtr;
}

InstructionPtr XMLParameterFile::
loadAMesh(const TiXmlElement* elem, map<string, MaterialDescPtr> materials, int & inOutControlVertexId) const
{
	assert(elem);
	assert(elem->Value() == string("Mesh"));
	
	InstructionPtr meshPtr;
	Rect3d rect;
    MaterialDescPtr permittivity, permeability;
    string material;
	
    if (sTryGetAttribute(elem, "permittivity", material))
    {
        if (materials.count(material) == 0)
            throw(std::logic_error(sErr("Unknown permittivity for Block", elem)));
        permittivity = materials[material];
    }
    if (sTryGetAttribute(elem, "permeability", material))
    {
        if (materials.count(material) == 0)
            throw(std::logic_error(sErr("Unknown permeability for Block", elem)));
        permeability = materials[material];
    }
    
    vector<Vector3d> vertices;
    vector<Vector3b> freeDirections;
    vector<Vector3i> faces;
    vector<int> controlVertexIds;
    
	const TiXmlElement* vertexXML = elem->FirstChildElement("Vertex");
    while (vertexXML != 0L)
    {
        Vector3d position;
        Vector3b freeDir;
        sGetMandatoryAttribute(vertexXML, "position", position);
        sGetOptionalAttribute(vertexXML, "freeDirections", freeDir,
            Vector3b(0,0,0));
        vertices.push_back(position);
        freeDirections.push_back(freeDir);
        controlVertexIds.push_back(inOutControlVertexId++);
        vertexXML = vertexXML->NextSiblingElement("Vertex");
    }
    
	const TiXmlElement* faceXML = elem->FirstChildElement("Face");
    while (faceXML != 0L)
    {
        Vector3i face;
        sGetMandatoryAttribute(faceXML, "vertices", face);
        faces.push_back(face);
        faceXML = faceXML->NextSiblingElement("Face");
    }
	
    try {
        meshPtr = InstructionPtr((Mesh*) new Mesh(
            vertices, controlVertexIds, freeDirections, faces, permittivity, permeability));
    } catch (std::logic_error & e) {
        throw(std::logic_error(sErr(e.what(), elem)));
    }
    	
	return meshPtr;
}

//InstructionPtr XMLParameterFile::
//loadAHeightMap(const TiXmlElement* elem, map<string, MaterialDescPtr> materials)
//    const
//{
//    assert(elem);
//    assert(elem->Value() == string("HeightMap"));
//    
//    InstructionPtr heightMap;
//    
//    Rect3d yeeRect;
//    MaterialDescPtr permittivity, permeability;
//    string material;
//    string imageFileName;
//    string rowDirStr, colDirStr, upDirStr;
//    Vector3i rowDir, colDir, upDir;
//    
//    if (sTryGetAttribute(elem, "permittivity", material))
//    {
//        if (materials.count(material) == 0)
//            throw(std::logic_error(sErr("Unknown permittivity for HeightMap", elem)));
//        permittivity = materials[material];
//    }
//    if (sTryGetAttribute(elem, "permeability", material))
//    {
//        if (materials.count(material) == 0)
//            throw(std::logic_error(sErr("Unknown permeability for HeightMap", elem)));
//        permeability = materials[material];
//    }
//    sGetMandatoryAttribute(elem, "yeeBounds", yeeRect);
//    sGetMandatoryAttribute(elem, "file", imageFileName);
//    sGetMandatoryAttribute(elem, "row", rowDirStr);
//    sGetMandatoryAttribute(elem, "column", colDirStr);
//    sGetMandatoryAttribute(elem, "up", upDirStr);
//    
//    try {
//        rowDir = sUnitVectorFromString(rowDirStr);
//        colDir = sUnitVectorFromString(colDirStr);
//        upDir = sUnitVectorFromString(upDirStr);
//        
//        heightMap = InstructionPtr(
//            new HeightMap(yeeRect, permittivity, permeability,
//                imageFileName, rowDir, colDir, upDir));
//    } catch (std::logic_error & e) {
//        throw(std::logic_error(sErr(e.what(), elem))); // append XML row/col info
//    }
//    
//    return heightMap;
//}

//InstructionPtr XMLParameterFile::
//loadAEllipsoid(const TiXmlElement* elem, map<string, MaterialDescPtr> materials)
//    const
//{
//    assert(elem);
//    assert(elem->Value() == string("Ellipsoid"));
//    
//    InstructionPtr ellipsoid;
//    Rect3d rect;
//    MaterialDescPtr permittivity, permeability;
//    string material;
//    
//    if (sTryGetAttribute(elem, "permittivity", material))
//    {
//        if (materials.count(material) == 0)
//            throw(std::logic_error(sErr("Unknown permittivity for Ellipsoid", elem)));
//        permittivity = materials[material];
//    }
//    if (sTryGetAttribute(elem, "permeability", material))
//    {
//        if (materials.count(material) == 0)
//            throw(std::logic_error(sErr("Unknown permeability for Ellipsoid", elem)));
//        permeability = materials[material];
//    }
//    
//    sGetMandatoryAttribute(elem, "yeeBounds", rect);
//    
//    try {
//        ellipsoid = InstructionPtr(
//            new Ellipsoid(rect, permittivity, permeability));
//    } catch (std::logic_error & e) {
//        throw(std::logic_error(sErr(e.what(), elem)));
//    }
//    
//    return ellipsoid;
//}

InstructionPtr XMLParameterFile::
loadACopyFrom(const TiXmlElement* elem, const set<string> & allGridNames) const
{
	assert(elem);
	assert(elem->Value() == string("CopyFrom"));
	
	InstructionPtr copyFrom;
	Rect3i sourceRect, destRect;
	string sourceGridName;
	
	sGetMandatoryAttribute(elem, "fromYeeCells", sourceRect);
	sGetMandatoryAttribute(elem, "toYeeCells", destRect);
	sGetMandatoryAttribute(elem, "sourceGrid", sourceGridName);
	if (allGridNames.count(sourceGridName) == 0)
		throw(std::logic_error(sErr("Unknown source grid name for CopyFrom", elem)));
	
	try {
		copyFrom = InstructionPtr(
			new CopyFrom(yeeToHalf(sourceRect), yeeToHalf(destRect),
                sourceGridName));
	} catch (std::logic_error & e) {
		throw(std::logic_error(sErr(e.what(), elem)));
	}
	
	return copyFrom;
}

InstructionPtr XMLParameterFile::
loadAnExtrude(const TiXmlElement* elem) const
{
	assert(elem);
	assert(elem->Value() == string("Extrude"));
	
	InstructionPtr extrude;
	Rect3i halfCellFrom, halfCellTo;
	
	sGetMandatoryAttribute(elem, "fromHalfCells", halfCellFrom);
	sGetMandatoryAttribute(elem, "toHalfCells", halfCellTo);
    
	try {
		extrude = InstructionPtr(
			new Extrude(
				halfCellFrom, halfCellTo));
	} catch (std::logic_error & e) {
		throw(std::logic_error(sErr(e.what(), elem)));
	}
	
	return extrude;
}

InstructionPtr XMLParameterFile::
loadABackground(const TiXmlElement* elem,
    map<string, MaterialDescPtr> materials) const
{
    assert(elem);
    assert(elem->Value() == string("Background"));
    
    InstructionPtr background;
    MaterialDescPtr permittivity, permeability;
    string material;
    
    if (sTryGetAttribute(elem, "permittivity", material))
    {
        if (materials.count(material) == 0)
            throw(std::logic_error(sErr("Unknown permittivity for Background", elem)));
        permittivity = materials[material];
    }
    
    if (sTryGetAttribute(elem, "permeability", material))
    {
        if (materials.count(material) == 0)
            throw(std::logic_error(sErr("Unknown permeability for Background", elem)));
        permeability = materials[material];
    }
    
    try {
        background = InstructionPtr(
            new Background(permittivity, permeability));
    } catch (std::logic_error & e) {
        throw(std::logic_error(sErr(e.what(), elem)));
    }
    
    return background;
}

InstructionPtr XMLParameterFile::
loadAPerm_ity(const TiXmlElement* elem,
    int numMeshes,
    const map<string, MaterialDescPtr> & materials) const
{
    assert(elem);
    assert(elem->Value() == string("Permittivity") || elem->Value() == string("Permeability"));
    
    InstructionPtr permity;
    
    std::string materialName;
    sGetMandatoryAttribute(elem, "materialName", materialName);
    
    if (0 == materials.count(materialName))
    {
        throw std::runtime_error("Unknown material");
    }
    
    int solidId;
    sGetMandatoryAttribute(elem, "solidId", solidId);
    
    if (solidId >= numMeshes)
    {
        throw std::runtime_error("solidId out of bounds");
    }
    
    permity = InstructionPtr( new Perm_ity(elem->Value(), materialName, solidId));
    return permity;
}


Duration XMLParameterFile::
loadADuration(const TiXmlElement* elem) const
{
    assert(elem);
    int firstTimestep;
    int lastTimestep;
    int period;
    
    Vector2d tInterval;
    
    sGetOptionalAttribute(elem, "firstTimestep", firstTimestep, 0);
    sGetOptionalAttribute(elem, "lastTimestep", lastTimestep, INT_MAX);
    sGetOptionalAttribute(elem, "period", period, 1);
    
    if (sTryGetAttribute(elem, "duration", tInterval))
        return Duration(firstTimestep, lastTimestep, tInterval, period);
    else
        return Duration(firstTimestep, lastTimestep, period);
}

Region XMLParameterFile::
loadARegion(const TiXmlElement* elem) const
{
    assert(elem);
    Rect3i yeeCells;
    Rect3d bounds;
    Vector3i stride;
    sGetMandatoryAttribute(elem, "yeeCells", yeeCells);
    sGetOptionalAttribute(elem, "stride", stride, Vector3i(1,1,1));
    if (sTryGetAttribute(elem, "bounds", bounds))
    {
        return Region(yeeCells, bounds, stride);
    }
    else
    {
        return Region(yeeCells, stride);
    }
}

Vector3i XMLParameterFile::
loadAnOmitSide(const TiXmlElement* elem) const
{
    Vector3i omitSide;
    Vector3i omitSideDominantComponent;
    istringstream str(elem->GetText());
    //elem->GetText() >> omitSide;
    str >> omitSide;
    if (!str.good())
        throw(std::logic_error(sErr("Can't read OmitSide", elem)));
    omitSideDominantComponent = dominantComponent(omitSide);
    
    if (omitSide != omitSideDominantComponent ||
        dot(omitSide, omitSide) != 1)
        throw(std::logic_error(sErr("OmitSide must specify an axis-oriented unit "
            "vector", elem)));
    
    return omitSide;
}


#pragma mark *** Static Method Implementations ***

static Vector3i sUnitVectorFromString(const string & axisString)
	throw(std::logic_error)
{
    Vector3i out(0,0,0);
    
    if (axisString == "x")
        out = Vector3i(1,0,0);
    else if (axisString == "-x")
        out = Vector3i(-1,0,0);
    else if (axisString == "y")
        out = Vector3i(0,1,0);
    else if (axisString == "-y")
        out = Vector3i(0, -1 ,0);
    else if (axisString == "z")
        out = Vector3i(0,0,1);
    else if (axisString == "-z")
        out = Vector3i(0,0,-1);
    else
        throw(std::logic_error("Unit vector string must be x, -x, y, -y, z, or -z"));
	
    return out;
}

