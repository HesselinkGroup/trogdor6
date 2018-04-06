/*
 *  main.cpp
 *  Trogdor6
 *
 *  Copyright 2018 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstdlib> // for random()
#include <unistd.h> // for sleep()

#include <boost/filesystem.hpp>

#include "Arguments.h"

#include "TimeWrapper.h"
#include "UserPreferences.h"

#include "SimulationDescription.h"
#include "YeeUtilities.h"
#include "XMLParameterFile.h"
#include "SimpleMesh/Search.h"

#include "VoxelizedPartition.h"
#include "GridOperations.h"
#include "ForwardGridOperations.h"
#include "AdjointGridOperations.h"
#include "GridFields.h"

using namespace std;
using namespace YeeUtilities;
namespace po = boost::program_options;

void makeOutputDirectory();

void voxelizeGrids(SimulationDescription & sim,
    map<GridDescPtr, VoxelizedPartitionPtr> & grids);

void
voxelizeGridRecursor(SimulationDescription & sim,
    map<GridDescPtr, VoxelizedPartitionPtr> & grids,
    VoxelizedPartition & vp,
    Vector3i numMPINodes,
    Vector3i thisNode,
    Rect3i partitionYeeCells);

VoxelizedPartitionPtr
makeSourceGrid(
    ElectromagneticMode mode,
    VoxelizedPartition & parentGrid,
	HuygensSurfaceDescription & huygensSurface, string srcGridName,
    Rect3i partitionMPIYee);
    
VoxelizedPartitionPtr
makeAuxGrid(Vector3i collapsible,
    ElectromagneticMode mode,
    VoxelizedPartition & parentGrid,
	HuygensSurfaceDescription & huygensSurface, string auxGridName,
    Rect3i partitionMPIYee);
    
void addAuxOutputs(GridDescPtr grid);

void prepareOperations(
    const map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids,
    map<GridDescPtr, GridOperationsPtr> & operations);

void allocateFields(
    const map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids,
    map<int, Pointer<GridFields> > & gridFields);
void setOperationPointers(
    map<GridDescPtr, GridOperationsPtr> & operations,
    map<int, Pointer<GridFields> > & gridFields);

void allocateOperationData(
    map<GridDescPtr, GridOperationsPtr> & operations);

void runSimulation(map<GridDescPtr, GridOperationsPtr> & operations,
    int numTimesteps);

static void setPreferences(po::variables_map arguments);

int main(int argc, char** argv)
{
#ifdef NDEBUG
    cout << "No assertions.\n";
    assert(0);
#else
    cout << "Debug mode.\n";
#endif
    
	po::variables_map variablesMap = handleArguments(argc, argv);
	if (variablesMap.count("help") || variablesMap.count("version") ||
		variablesMap.count("numerics") )
		return 0;
    
    setPreferences(variablesMap);
    
    makeOutputDirectory();
    
    string paramFileName("params.xml");
    if (UserPreferences::defines("input-file"))
        paramFileName = UserPreferences::value("input-file");
    
    XMLParameterFile xml(paramFileName);
    SimulationDescription sim(xml);
    
    map<GridDescPtr, VoxelizedPartitionPtr> voxelizedGrids;
    map<GridDescPtr, GridOperationsPtr> operations;
    map<int, Pointer<GridFields> > fields;
    
    voxelizeGrids(sim, voxelizedGrids);
    prepareOperations(voxelizedGrids, operations);
    // delete voxelized grids here?
    allocateFields(voxelizedGrids, fields);
    setOperationPointers(operations, fields);
    allocateOperationData(operations);
    
    if (!variablesMap.count("nosim"))
    {
        // open IO files here.  no reason to crash earlier.
                
        runSimulation(operations, sim.numTimesteps());
    }
    
    map<GridDescPtr, GridOperationsPtr>::iterator itr;
    for (itr = operations.begin(); itr != operations.end(); itr++)
    {
        itr->second->printPerformance(sim.numTimesteps());
    }
    
    return 0;
}

void setPreferences(po::variables_map arguments)
{
    if (arguments.count("input-file"))
        UserPreferences::set("input-file",
            arguments["input-file"].as<string>());
    if (arguments.count("adjoint"))
    {
        LOG << "Adjoint mode.\n";
        UserPreferences::set("adjoint");
    }
    if (arguments.count("noAveraging"))
    {
        LOG << "Permittivity averaging off.\n";
        UserPreferences::set("noAveraging");
    }
    if (arguments["materials"].as<string>() == "mixed")
        UserPreferences::set("materials", "mixed");
    else if (arguments["materials"].as<string>() == "harmonic")
        UserPreferences::set("materials", "harmonic");
    else if (arguments["materials"].as<string>() == "arithmetic")
        UserPreferences::set("materials", "arithmetic");
    else
    {
        cerr << "Invalid materials option: "
            << arguments["materials"].as<string>() << "\n"
            << "Allowed values are mixed, harmonic, and arithmetic.\n";
        exit(1);
    }
    
    if (arguments.count("noNormalizedOrientation"))
        UserPreferences::set("noNormalizedOrientation");
    
    if (arguments.count("sensitivity"))
        UserPreferences::set("sensitivity");
    if (arguments.count("savePermittivity"))
        UserPreferences::set("savePermittivity");
    if (arguments.count("saveOrientation"))
        UserPreferences::set("saveOrientation");
    if (arguments.count("saveOrientationSensitivity"))
        UserPreferences::set("saveOrientationSensitivity");
    if (arguments.count("saveFillFactors"))
        UserPreferences::set("saveFillFactors");
    if (arguments.count("saveFillFactorSensitivity"))
        UserPreferences::set("saveFillFactorSensitivity");
    
    if (arguments.count("geometry"))
        UserPreferences::set("geometry");
    if (arguments.count("derivDelta"))
        UserPreferences::set<double>("derivDelta", arguments["derivDelta"].as<double>());
    if (arguments.count("oSmoothing"))
        UserPreferences::set<int>("oSmoothing", arguments["oSmoothing"].as<int>());
    if (arguments.count("fSmoothing"))
        UserPreferences::set<int>("fSmoothing", arguments["fSmoothing"].as<int>());
    if (arguments.count("omitEdgeIntegral"))
        UserPreferences::set("omitEdgeIntegral");
    if (arguments.count("fullTensor"))
        UserPreferences::set("fullTensor");
    if (arguments.count("flipOffDiagonalSign"))
        UserPreferences::set("flipOffDiagonalSign");
    if (arguments.count("excludeCorners"))
        UserPreferences::set("excludeCorners");
        UserPreferences::set("offDiagonalOctant",
        arguments["offDiagonalOctant"].as<string>());
    
    if (arguments.count("noPMLE"))
        UserPreferences::set("noPMLE");
    if (arguments.count("noPMLH"))
        UserPreferences::set("noPMLH");
    if (arguments.count("printRunlines"))
        UserPreferences::set("printRunlines");
    
    if (arguments.count("outputAuxFields"))
        UserPreferences::set("outputAuxFields");
    if (arguments.count("outputAuxCurrents"))
        UserPreferences::set("outputAuxCurrents");
    
    if (arguments.count("outputDirectory"))
        UserPreferences::set("outputDirectory",
        arguments["outputDirectory"].as<string>());
    
    if (arguments.count("booleans"))
        UserPreferences::set("booleans");
}



void makeOutputDirectory()
{
    string outDir = UserPreferences::valueAs<string>("outputDirectory");
    
    boost::filesystem::path outputDirPath(outDir.c_str());
    
    LOG << "Considering output directory " << outDir << "\n";
    LOG << "Boost path is " << outputDirPath << "\n";
    
    if (!boost::filesystem::exists(outputDirPath))
    {
        if (!boost::filesystem::create_directory(outputDirPath))
        {
            throw(std::runtime_error("Could not create output directory!"));
        }
        LOG << "Created output directory.\n";
    }
}



// TODO: move up to caller
static void sWriteForMatlab(std::string fileName, Vector3d origin,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices,
    const std::vector<std::vector<SimpleMesh::Triangle> > & meshes)
{
    ofstream str(fileName.c_str());
    str << "function polyhedra = mesh()\n";
    str << "% Get triangulated permittivity and permeability mesh with "
        "control vertices.\n";
    str << "% polyhedra = mesh() returns a structure with members:\n";
    str << "%   polyhedra.permittivity\n";
    str << "%   polyhedra.permeability\n";
    
//    map<string, TrackedPolyhedron*>::const_iterator itr;
    
    // First, control vertices.
    str << "polyhedra.controlVertices = [...\n";
    for (int nn = 0; nn < controlVertices.size(); nn++)
    {
        // FIXME: At this point, the control vertex IDs are all -1, which is STUPID.
        Vector3d point(controlVertices[nn].point() + origin);
        str << "\t" << point[0] << " "
            << point[1] << " "
            << point[2] << ";...\n";
    }
    str << "];\n";
    
    // First permittivity.
//    itr = polyhedra().begin();
//    for (int nn = 0; nn < polyhedra().size(); nn++, itr++)
    for (int nn = 0; nn < meshes.size(); nn++)
    {
        const std::vector<SimpleMesh::Triangle> & mesh = meshes.at(nn);
        
        str << "polyhedra.permittivity{" << nn+1 << "}.name = '"
            << nn << "';\n";
        
        str << "polyhedra.permittivity{" << nn+1
            << "}.controlVertexIds = [...\n";
//        vector<SimpleMesh::Triangle>::const_iterator itr2;
        for (int tt = 0; tt < mesh.size(); tt++)
        {
            str << "\t" << mesh[tt].controlVertices()[0].id() + 1 << " "
                << mesh[tt].controlVertices()[1].id() + 1 << " "
                << mesh[tt].controlVertices()[2].id() + 1<< ";...\n";
        }
        str << "];\n";
        
        str << "polyhedra.permittivity{" << nn+1
            << "}.vertices = [...\n";
        for (int mm = 0; mm < mesh.size(); mm++)
        {
            for (int tt = 0; tt < 3; tt++)
            {
                str << "\t" << mesh[mm].triangle()[tt][0] + origin[0] << " "
                    << mesh[mm].triangle()[tt][1] + origin[1] << " "
                    << mesh[mm].triangle()[tt][2] + origin[2] << ";...\n";
            }
        }
        str << "];\n";
        
        str << "polyhedra.permittivity{" << nn+1
            << "}.faces = [...\n";
        for (int mm = 0; mm < mesh.size(); mm++)
        {
            str << "\t" << 3*mm+1 << " " << 3*mm+2 << " " << 3*mm+3 << ";...\n";
        }
        str << "];\n";
    }
}


void doNonBoolean(const AssemblyDescription & assemblyDesc,
    const std::map<std::string, Precision::RationalFunction> & padeCoeffs,
    std::vector<Precision::RationalFunction> & outPermittivities,
    std::vector<Precision::RationalFunction> & outPermeabilities,
    std::vector<std::vector<SimpleMesh::Triangle> > & outPermittivityMesh,
    std::vector<std::vector<SimpleMesh::Triangle> > & outPermeabilityMesh)
{
    const vector<InstructionPtr> & instructions(assemblyDesc.instructions());
    const std::vector<std::vector<SimpleMesh::Triangle> > & inputMeshes = assemblyDesc.meshes();
    
    // Accumulate the meshes for all structures here.
    std::map<std::string, std::vector<SimpleMesh::Triangle> > epsMeshes, muMeshes;
    
    for (unsigned int mm = 0; mm < instructions.size(); mm++)
    if (instructions[mm]->type() == kPermityType)
    {
        const Perm_ity & chunk = (const Perm_ity&)*instructions[mm];
        const std::vector<SimpleMesh::Triangle> & chunkTriangles = inputMeshes.at(chunk.getSolidId());
        
        std::map<std::string, std::vector<SimpleMesh::Triangle> >* tmpMeshes;
        
        if (chunk.getType() == "Permittivity")
        {
            tmpMeshes = &epsMeshes;
        }
        else if (chunk.getType() == "Permeability")
        {
            tmpMeshes = &muMeshes;
        }
        else
        {
            throw(std::runtime_error("Unknown chunk type"));
        }
        
        std::vector<SimpleMesh::Triangle> & tris = (*tmpMeshes)[chunk.getMaterialName()];
        tris.insert(tris.end(), chunkTriangles.begin(), chunkTriangles.end());
    }
    
    // Dump the map<material, mesh> structures into vector<material> and vector<mesh>.
    // This really should be like [keys, values] = map.getKeyValuePairs().  :-/
    
    std::map<std::string, std::vector<SimpleMesh::Triangle> >::iterator itr;
    for (itr = epsMeshes.begin(); itr != epsMeshes.end(); itr++)
    {
        SimpleMesh::determineNeighbors(itr->second);
        SimpleMesh::determineEdgeControlVertices(itr->second);
        outPermittivities.push_back(padeCoeffs.find(itr->first)->second);
        outPermittivityMesh.push_back(itr->second);
    }
    for (itr = muMeshes.begin(); itr != muMeshes.end(); itr++)
    {
        SimpleMesh::determineNeighbors(itr->second);
        SimpleMesh::determineEdgeControlVertices(itr->second);
        outPermeabilities.push_back(padeCoeffs.find(itr->first)->second);
        outPermeabilityMesh.push_back(itr->second);
    }
}


void
voxelizeGrids(SimulationDescription & sim,
    map<GridDescPtr, VoxelizedPartitionPtr> & grids)
{
    Rect3i partitionMPIYee(-1000000, -1000000, -1000000, 1000000, 1000000, 1000000);
    Vector3i numNodes(1,1,1);
    Vector3i thisNode(0,0,0);
    
    LOG << "Recursively voxelize grids.\n";
    
    for (int nn = 0; nn < sim.grids().size(); nn++)
    {
        GridDescPtr currentGrid = sim.grids()[nn];
        Rect3i partitionYeeCells(intersection(partitionMPIYee,
            currentGrid->extents().physicalYeeCells()));
        
        std::vector<Precision::RationalFunction> permittivities, permeabilities;
        std::vector<std::vector<SimpleMesh::Triangle> > permittivityMesh, permeabilityMesh;
        
        if (UserPreferences::defines("booleans"))
        {
            throw(std::runtime_error("Boolean operations not supported (use external tool!)"));
//            std::cout << "Using Booleans.\n";
//
//            std::cerr << "How many meshes? " << currentGrid->assembly()->meshes().size() << ".\n";
//
//            doBoolean(*currentGrid->assembly(), sim.materialPadeCoefficients(),
//                permittivities, permeabilities, permittivityMesh, permeabilityMesh);
        }
        else
        {
            std::cout << "Not using Booleans.\n";
            
            doNonBoolean(*currentGrid->assembly(), sim.materialPadeCoefficients(),
                permittivities, permeabilities, permittivityMesh, permeabilityMesh);
        }
        
        if (UserPreferences::defines("geometry"))
        {
            std::string currentGridName("Main");
            sWriteForMatlab(currentGridName + ".m", currentGrid->origin(),
                currentGrid->assembly()->controlVertices(), permittivityMesh);
            // TODO: permeability geometry export as well
        }
        
        VoxelizedPartition* vp = new VoxelizedPartition(currentGrid,
            currentGrid->extents(),
            partitionYeeCells,
            permittivities,
            permeabilities,
            permittivityMesh,
            permeabilityMesh,
            currentGrid->assembly()->controlVertices());
        
        grids[currentGrid] = Pointer<VoxelizedPartition>(vp);
        
        voxelizeGridRecursor(sim, grids, *vp, numNodes, thisNode,
            partitionMPIYee);
    }
}

void
voxelizeGridRecursor(SimulationDescription & sim,
    map<GridDescPtr, VoxelizedPartitionPtr> & grids,
    VoxelizedPartition & vp,
    Vector3i numMPINodes,
    Vector3i thisNode,
    Rect3i partitionMPIYee)
{
    // Next, iterate through all the Huygens surfaces and look for TFSF sources
    // that ought to generate more VoxelizedPartitions.
    
    for (int ss = 0; ss < vp.gridDescription()->huygensSurfaces().size(); ss++)
    if (vp.gridDescription()->huygensSurfaces()[ss]->type() == kTFSFSource)
    {
        HuygensSurfaceDescription & huygSurf =
            *vp.gridDescription()->huygensSurfaces()[ss];
        Vector3i srcSymmetry(huygSurf.symmetries());
        Vector3i gridSymmetry(vp.boundarySymmetry(huygSurf.yeeCells()));
        Vector3i collapsibleDimensions = srcSymmetry*gridSymmetry*
            vp.gridDescription()->extents().physicalYeeCells().size();
        // the elements of collapsibleDimensions are interpreted as true or
        // false.
        
        LOG << "Source: " << srcSymmetry << " grid: " << gridSymmetry << "\n";
        LOG << "Size: " << vp.gridDescription()->extents().physicalYeeCells().size()
            << "\n";
        LOG << "Collapsible dimensions " << collapsibleDimensions << "\n";
        if (norm(collapsibleDimensions) > 0)
        {
            LOG << "Making aux grid.\n";
            VoxelizedPartitionPtr vpAux = makeAuxGrid(
                srcSymmetry*gridSymmetry, sim.electromagneticMode(),
                vp, huygSurf,
                vp.gridDescription()->name() + "_autoaux",
                partitionMPIYee);
            grids[vpAux->gridDescription()] = vpAux;
            // the huygens surface was changed into a link by the function call.
            
            voxelizeGridRecursor(sim, grids, *vpAux, numMPINodes, thisNode,
                partitionMPIYee);
        }
        else
        {
            LOG << "Making source grid.\n";
            if (vp.gridDescription()->numDimensions() > 1)
                throw(std::logic_error("Grid dimension too high"));
            VoxelizedPartitionPtr vpSrc = makeSourceGrid(
                sim.electromagneticMode(),
                vp, huygSurf,
                vp.gridDescription()->name() + "_autosrc",
                partitionMPIYee);
            grids[vpSrc->gridDescription()] = vpSrc;
            voxelizeGridRecursor(sim, grids, *vpSrc, numMPINodes, thisNode,
                partitionMPIYee);
        }
    }
}

VoxelizedPartitionPtr
makeAuxGrid(Vector3i collapsible, ElectromagneticMode mode,
    VoxelizedPartition & parentGrid,
	HuygensSurfaceDescription & huygensSurface, string auxGridName,
    Rect3i partitionMPIYee)
{
    Mat3i collapser(Mat3i::diagonal(!collapsible));
	const set<Vector3i> & omittedSides = huygensSurface.omittedSides();
    
	// What does the aux grid look like?
	// It's the size of the original total field region, collapsed, and all
	// non-1D dimensions get 10 cells of PML.
    
    Rect3i linkSourceYeeCells(collapser*huygensSurface.yeeCells());
    
    Rect3i totalFieldYeeCells(inset(huygensSurface.yeeCells(), -1));
    Rect3i nonPMLYeeCells(inset(totalFieldYeeCells, -1));
    Rect3i physicalYeeCells(inset(nonPMLYeeCells, -10));
    
    totalFieldYeeCells = collapser*totalFieldYeeCells;
    nonPMLYeeCells = collapser*nonPMLYeeCells;
    physicalYeeCells = collapser*physicalYeeCells;
    
    GridDescPtr childGridDesc(new GridDescription(auxGridName,
        mode,
        physicalYeeCells,
        yeeToHalf(nonPMLYeeCells),
        parentGrid.gridDescription()->dxyz(),
        parentGrid.gridDescription()->dt(),
        parentGrid.gridDescription()->numTimesteps()));
    addAuxOutputs(childGridDesc);
    
	HuygensSurfaceDescPtr childSource;
	if (huygensSurface.type() == kTFSFSource)
	{
		//	Omitting the "back side" is harmless in aux grids, and desirable
		//  for 1D aux grids in the terminal recursion.
		set<Vector3i> newSourceOmittedSides(omittedSides);
		newSourceOmittedSides.insert(huygensSurface.direction());
        if (newSourceOmittedSides.size() >= 6)
            throw(std::logic_error("Huygens surface omits all sides"));
        
        if (huygensSurface.formula() != "")
        {
            childSource = HuygensSurfaceDescPtr(HuygensSurfaceDescription::
                newTFSFFormulaSource(huygensSurface.sourceFields(),
                    huygensSurface.formula(),
                    huygensSurface.direction(),
                    totalFieldYeeCells,
                    newSourceOmittedSides,
                    huygensSurface.isTotalField()));
            
            // TODO: huygensSurface.newTFSFFormulaSource(totalFieldYeeCells, newSourceOmittedSides)
        }
        else
        {
            childSource = HuygensSurfaceDescPtr(HuygensSurfaceDescription::
                newTFSFTimeSource(huygensSurface.sourceFields(),
                    huygensSurface.timeFile(),
                    huygensSurface.direction(),
                    totalFieldYeeCells,
                    newSourceOmittedSides,
                    huygensSurface.isTotalField()));
            
            // TODO: huygensSurface.newTFSFTimeSource(totalFieldYeeCells, newSourceOmittedSides)
        }
    }
	else if (huygensSurface.type() == kCustomTFSFSource)
	{
		childSource = HuygensSurfaceDescPtr(HuygensSurfaceDescription::
            newCustomTFSFSource(
			huygensSurface.file(),
			huygensSurface.symmetries(),
            totalFieldYeeCells,
            huygensSurface.duration(),
			huygensSurface.omittedSides(),
            huygensSurface.isTotalField()));
        
        // TODO: huygensSurface.newCustomTFSFSource(totalFieldYeeCells)
	}
	
	childGridDesc->setHuygensSurfaces(
        vector<HuygensSurfaceDescPtr>(1,childSource));
    
    // And last of all, turn the parent source into a link to the child grid.
    huygensSurface.becomeLink(childGridDesc, linkSourceYeeCells);
    
	Rect3i copyYeeCells(huygensSurface.yeeCells());
	for (int xyz = 0; xyz < 3; xyz++)
	if (omittedSides.count(cardinal(xyz*2)))
		copyYeeCells.p1[xyz] = copyYeeCells.p2[xyz];
	else if (omittedSides.count(cardinal(xyz*2+1)))
		copyYeeCells.p2[xyz] = copyYeeCells.p1[xyz];
	else if (collapsible[xyz] != 0)
		copyYeeCells.p1[xyz] = copyYeeCells.p2[xyz]; // either side is ok here!
    
    Rect3i partitionYeeCells(intersection(partitionMPIYee,
        childGridDesc->extents().physicalYeeCells()));
    VoxelizedPartitionPtr vp(new VoxelizedPartition(
        childGridDesc, parentGrid, copyYeeCells, childGridDesc->extents(),
        partitionYeeCells));
	
    return vp;
}

VoxelizedPartitionPtr
makeSourceGrid(
    ElectromagneticMode mode,
    VoxelizedPartition & parentGrid,
	HuygensSurfaceDescription & huygensSurface, string srcGridName,
    Rect3i partitionMPIYee)
{
	assert(huygensSurface.type() == kTFSFSource);
    
	// What does the source grid look like?
	// Two cases:
	//   1.  Back side of source region is omitted.  Then the source is tiny,
	//   just large enough to provide the fields to handle the front boundary.
	//   2.  Back side of source region is not omitted.  Then we clone the
	//   whole space and make a source with omitted back side.
    
    Rect3i totalFieldYeeCells(huygensSurface.yeeCells());
    if (huygensSurface.omittedSides().count(huygensSurface.direction()))
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (huygensSurface.direction()[xyz] > 0)
            totalFieldYeeCells.p2[xyz] = totalFieldYeeCells.p1[xyz];
        else if (huygensSurface.direction()[xyz] < 0)
            totalFieldYeeCells.p1[xyz] = totalFieldYeeCells.p2[xyz];
    }
    
    totalFieldYeeCells.p1 -= vec_abs(huygensSurface.direction());
    totalFieldYeeCells.p2 += vec_abs(huygensSurface.direction());
    
    Rect3i nonPMLYeeCells(totalFieldYeeCells);
    nonPMLYeeCells.p1 -= vec_abs(huygensSurface.direction());
    nonPMLYeeCells.p2 += vec_abs(huygensSurface.direction());
    
    Rect3i physicalYeeCells(nonPMLYeeCells);
    physicalYeeCells.p1 -= 10*vec_abs(huygensSurface.direction()); // PML x10
    physicalYeeCells.p2 += 10*vec_abs(huygensSurface.direction());
    
	GridDescPtr childGridDesc(new GridDescription(srcGridName,
        mode,
		physicalYeeCells,
		yeeToHalf(nonPMLYeeCells),
        parentGrid.gridDescription()->dxyz(),
        parentGrid.gridDescription()->dt(),
        parentGrid.gridDescription()->numTimesteps()));
	addAuxOutputs(childGridDesc);
    
	if (huygensSurface.omittedSides().count(huygensSurface.direction()))
//    if (1)
	{
        LOG << "Using hard source.\n";
        const int USE_HARD_SOURCE = 0;
        
        // Put the source just outside the edge of the total field region.
        // Here's a sneaky way to position it.
        Vector3i sourcePosition = clip(totalFieldYeeCells,
            -100000*huygensSurface.direction()) - huygensSurface.direction();
		
        LOG << "\tPosition " << sourcePosition << "\n";
        
        vector<Region> regions(1, Region(
            Rect3i(sourcePosition, sourcePosition)));
        
		SourceDescPtr sPtr(new SourceDescription(
            huygensSurface.sourceFields(),
			huygensSurface.formula(),
			huygensSurface.timeFile(),
            "", // HuygensSurface has no spaceFile (mask)
			"", // HuygensSurface has no spaceTimeFile (TFSFSource anyway)
            USE_HARD_SOURCE,
            regions,
            vector<Duration>(1,huygensSurface.duration())));
        childGridDesc->setSources(vector<SourceDescPtr>(1, sPtr));
	}
	else
	{
        LOG << "Using soft source.\n";
		// use a soft source, and omit the back side
		HuygensSurfaceDescPtr hPtr(new HuygensSurfaceDescription(
            huygensSurface));
        hPtr->yeeCells(totalFieldYeeCells);
        hPtr->omitSide(huygensSurface.direction());
        childGridDesc->setHuygensSurfaces(vector<HuygensSurfaceDescPtr>(1, hPtr));
	}
    
    // And last of all, turn the parent source into a link to the child grid.
    huygensSurface.becomeLink(childGridDesc, huygensSurface.yeeCells());
    
	Rect3i copyYeeCells(huygensSurface.yeeCells());
	for (int xyz = 0; xyz < 3; xyz++)
	if (huygensSurface.omittedSides().count(cardinal(xyz*2)))
		copyYeeCells.p1[xyz] = copyYeeCells.p2[xyz];
	else if (huygensSurface.omittedSides().count(cardinal(xyz*2+1)))
		copyYeeCells.p2[xyz] = copyYeeCells.p1[xyz];
    
    Rect3i partitionYeeCells(intersection(partitionMPIYee,
        childGridDesc->extents().physicalYeeCells()));
    VoxelizedPartitionPtr vp(new VoxelizedPartition(
        childGridDesc, parentGrid, copyYeeCells,
        childGridDesc->extents(), partitionYeeCells));
	
	return vp;
}

void addAuxOutputs(GridDescPtr grid)
{
    if (UserPreferences::defines("outputAuxFields"))
    {
        Region everywhere(grid->extents().physicalYeeCells());
        ostringstream fname;
        fname << grid->name() << "_eh";
        OutputDescPtr outputDesc(new OutputDescription("ex ey ez hx hy hz",
            fname.str(), everywhere));
        grid->pushOutput(outputDesc);
    }
    if (UserPreferences::defines("outputAuxCurrents"))
    {
        //Region everywhere(grid->extents().physicalYeeCells());
        Region everywhere(halfToYee(grid->extents().nonPMLHalfCells()));
        ostringstream fname;
        fname << grid->name() << "_jm";
        OutputDescPtr outputDesc(new OutputDescription("jx jy jz mx my mz",
            fname.str(), everywhere));
        grid->pushOutput(outputDesc);
    }
}

void prepareOperations(
    const map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids,
    map<GridDescPtr, GridOperationsPtr> & operations)
{
    map<GridDescPtr, VoxelizedPartitionPtr>::const_iterator itr;
    for (itr = voxelizedGrids.begin(); itr != voxelizedGrids.end(); itr++)
    {
        LOG << "Operations for " << itr->first->name() << "\n";
        
        GridOperations* ops;
        
        if (UserPreferences::defines("adjoint"))
        {
            ops = new AdjointGridOperations(itr->first, voxelizedGrids);
        }
        else
        {
            ops = new ForwardGridOperations(itr->first, voxelizedGrids);
        }
        
        operations[itr->first] = Pointer<GridOperations>(ops);
    }
}

void allocateFields(
    const map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids,
    map<int, Pointer<GridFields> > & gridFields)
{
    long bytes = 0, totalBytes = 0;
    map<GridDescPtr, VoxelizedPartitionPtr>::const_iterator itr;
    for (itr = voxelizedGrids.begin(); itr != voxelizedGrids.end(); itr++)
    {
        LOG << "Allocating fields and currents for " << itr->first->name() << "\n";
        gridFields[itr->first->id()] = Pointer<GridFields>(
            new GridFields(itr->second->fieldIndices()));
        bytes = gridFields[itr->first->id()]->bytes();
        LOGMORE << "\t" << bytes/1024.0 << " kB\n";
        totalBytes += bytes;
    }
    
}

void setOperationPointers(
    map<GridDescPtr, GridOperationsPtr> & operations,
    map<int, Pointer<GridFields> > & gridFields)
{
    map<GridDescPtr, GridOperationsPtr>::iterator itr;
    for (itr = operations.begin(); itr != operations.end(); itr++)
    {
        LOG << "Setting pointers for " << itr->first->name() << "\n";
        itr->second->setPointers(gridFields);
        LOGMORE << "\t" << itr->second->bytes()/1024.0 << " kB\n";
    }
}

void allocateOperationData(
    map<GridDescPtr, GridOperationsPtr> & operations)
{
    map<GridDescPtr, GridOperationsPtr>::iterator itr;
    for (itr = operations.begin(); itr != operations.end(); itr++)
    {
        LOG << "Allocating aux data for " << itr->first->name() << "\n";
        itr->second->allocate();
        LOGMORE << "\t" << itr->second->bytes()/1024.0 << " kB\n";
        if (UserPreferences::defines("printRunlines"))
        {
            LOG << "Printing runlines:\n";
            itr->second->printRunlines();
        }
    }
}

void runSimulation(map<GridDescPtr, GridOperationsPtr> & operations,
    int numTimesteps)
{
    map<GridDescPtr, GridOperationsPtr>::iterator itr;
    for (int tt = 0; tt < numTimesteps; tt++)
    {
        cout << "\r                                                          "
            << flush;
        cout << "\rTimestep " << tt << " of " << numTimesteps << flush;
        for (itr = operations.begin(); itr != operations.end(); itr++)
            itr->second->firstHalfStep(tt, itr->first->dt());
        for (itr = operations.begin(); itr != operations.end(); itr++)
            itr->second->secondHalfStep(tt, itr->first->dt());
        for (itr = operations.begin(); itr != operations.end(); itr++)
            itr->second->output(tt);
    }
    cout << endl;
    
    
}

