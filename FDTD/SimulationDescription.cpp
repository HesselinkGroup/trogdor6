/*
 *  SimulationDescription.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 1/31/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "SimulationDescription.h"
#include "XMLParameterFile.h"
//#include "STLOutput.h"

// This is included just in order to validate PML formulae
//#include "calc.hh"

#include <iostream>
#include <sstream>
#include <stdexcept>

// This is included for the function that converts hex colors to RGB.

#include "YeeUtilities.h"
using namespace YeeUtilities;
using namespace std;

namespace P = Precision;

SimulationDescription::
SimulationDescription(const XMLParameterFile & file) :
    mElectromagneticMode(k3d),
	m_dt(0.0f),
	m_dxyz(P::Vec3(0.0f, 0.0f, 0.0f)),
	mNumTimesteps(0)
{
	file.load(*this);
	
	// This error checking is redundant provided that the loading mechanism
	// remembers to check the timestep.  But it's the constructor's job to
	// make sure it returns a *valid* data structure, so let's double-check
	// these silly conditions here.
	if (m_dt <= 0)
		throw(std::logic_error("Nonpositive timestep"));
	if (!vec_gt(m_dxyz, 0.0f))
		throw(std::logic_error("Nonpositive cell dimensions"));
	if (mNumTimesteps < 0)
		throw(std::logic_error("Negative number of iterations"));
}

void SimulationDescription::
setMaterials(const vector<MaterialDescPtr> & materials)
{
    for (int nn = 0; nn < materials.size(); nn++)
        mMaterials[materials[nn]->name()] = materials[nn];
}

void SimulationDescription::
setGridPointers()
{
	unsigned int nGrid; //, nMaterial;
	
	// 1.  Make maps of names to pointers
	Map<string, GridDescPtr> gridMap;
	
	for (nGrid = 0; nGrid < mGrids.size(); nGrid++)
		gridMap[mGrids[nGrid]->name()] = mGrids[nGrid];
	
	// 2.  Point Huygens surfaces to appropriate grids
	for (nGrid = 0; nGrid < mGrids.size(); nGrid++)
		mGrids[nGrid]->setPointers(gridMap);
}

void SimulationDescription::
setDiscretization(Precision::Vec3 dxyz, Precision::Float dt)
{
	if (!vec_gt(dxyz, 0.0f))
		throw(std::logic_error("Nonpositive cell dimensions"));
	if (dt <= 0)
		throw(std::logic_error("Nonpositive timestep"));
	
	m_dxyz = dxyz;
	m_dt = dt;
}

void SimulationDescription::
setNumTimesteps(int numT)
{
	if (numT < 0)
		throw(std::logic_error("Number of timesteps must be nonnegative"));
	mNumTimesteps = numT;
}

std::map<string, P::RationalFunction> SimulationDescription::
materialPadeCoefficients() const
{
    map<string, P::RationalFunction> outMaterials;
    map<string, MaterialDescPtr>::const_iterator itr;
    for (itr = mMaterials.begin(); itr != mMaterials.end(); itr++)
    {
        outMaterials[itr->first] = itr->second->padeCoefficients();
    }
    return outMaterials;
}

#pragma mark *** Grid ***

int GridDescription::sNextGridID = 0;

GridDescription::
GridDescription(string name, ElectromagneticMode mode, Rect3i yeeCells, Rect3i nonPMLHalfCells,
    Precision::Vec3 dxyz, Precision::Float dt, long numTimesteps):
	mName(name),
    mID(sNextGridID++),
    mExtents(yeeCells, nonPMLHalfCells),
    mElectromagneticMode(mode),
    mDxyz(dxyz),
    mDt(dt),
    mNumTimesteps(numTimesteps)
{
	if (!vec_ge(mExtents.nonPMLHalfCells().size(), 0))
		throw(std::logic_error("Non-PML region has some negative dimensions"));
    
    vector<BoundaryType> boundaryType(6);
    
    // Set the boundary types on all six faces.  Currently kPECBoundary
    // and kPMCBoundary cannot be turned on from the XML parameter file.
    // However, if PML is only put on one side of the grid, say -x, then
    // PEC or PMC will be put on the other side automatically.
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (mExtents.physicalYeeCells().num()[xyz] == 1)
        {
            boundaryType[2*xyz] = kTranslationSymmetricBoundary;
            boundaryType[2*xyz+1] = kTranslationSymmetricBoundary;
        }
        else if (mExtents.nonPMLHalfCells().num(xyz) !=
            yeeToHalf(extents().physicalYeeCells()).num(xyz))
        {
            if (mExtents.nonPMLHalfCells().p1[xyz] >
                yeeToHalf(extents().physicalYeeCells()).p1[xyz])
            {
                boundaryType[2*xyz] = kPMLBoundary;
            }
            else
            {
                boundaryType[2*xyz] = kPECBoundary;
            }
            
            if (mExtents.nonPMLHalfCells().p2[xyz] <
                yeeToHalf(extents().physicalYeeCells()).p2[xyz])
            {
                boundaryType[2*xyz+1] = kPMLBoundary;
            }
            else
            {
                boundaryType[2*xyz+1] = kPMCBoundary;
            }
        }
        else
        {
            boundaryType[2*xyz] = kPeriodicBoundary;
            boundaryType[2*xyz+1] = kPeriodicBoundary;
        }
    }
    mExtents.physicalBoundaries(boundaryType);
}

void GridDescription::
setGridReports(const std::vector<GridReportDescPtr> & gridReports)
{
    mGridReports = gridReports;
    
    //LOG << "Adding " << gridReports.size() << " grid reports.\n";
    for (int nn = 0; nn < gridReports.size(); nn++)
    {
        if (mGridReports[nn]->regions().size() == 0)
            mGridReports[nn]->regions().push_back(
                Region(extents().physicalYeeCells()));
    }
}

void GridDescription::
setHuygensSurfaces(const vector<HuygensSurfaceDescPtr> & surfaces)
{
	mHuygensSurfaces = surfaces;
	// We need to make sure that Huygens surfaces in repeating grids will
	// omit any sides that are adjacent due to periodicity.  That means, if the
	// Huygens surface extends the full dimension of the grid, don't make a
	// correction.
	for (unsigned int nn = 0; nn < mHuygensSurfaces.size(); nn++)
	{
		//Rect3i destYeeCells = mHuygensSurfaces[nn]->yeeCells();
        //Vector3i numHalfCells = yeeToHalf(extents().physicalYeeCells()).num();
		
		for (int xyz = 0; xyz < 3; xyz++)
        if (mHuygensSurfaces[nn]->yeeCells().num(xyz) ==
            extents().physicalYeeCells().num(xyz))
        {
            mHuygensSurfaces[nn]->omitSide(2*xyz); // the -x, -y or -z side
            mHuygensSurfaces[nn]->omitSide(2*xyz+1); // the +x, +y or +z side
        }
        
//		{
//			if (destHalfRect.p1[xyz] <= 0 &&
//				destHalfRect.p2[xyz] >= numHalfCells[xyz]-1)
//			{
//				mHuygensSurfaces[nn]->omitSide(2*xyz); // the -x, -y or -z side
//				mHuygensSurfaces[nn]->omitSide(2*xyz+1); // the +x, +y or +z side
//			}
//		}
	}
}


void GridDescription::
setPMLParams(const Map<string, string> & p)
   
{
    string errorString;
    bool validParams;
    
//    validParams = sValidPMLParams(p, errorString);
//
//    if (!validParams)
//        throw(std::logic_error(errorString));
    
    mPMLParams = p;
}

void GridDescription::
setPointers(const Map<string, GridDescPtr> & gridMap)
{
	for (unsigned int nn = 0; nn < mHuygensSurfaces.size(); nn++)
	if (mHuygensSurfaces[nn]->type() == kLink)
		mHuygensSurfaces[nn]->setPointers(gridMap);
	
	if (mAssembly != 0L)
		mAssembly->setPointers(gridMap);
}

int GridDescription::
numDimensions() const
{
	int nDim = 0;
	for (int xyz = 0; xyz < 3; xyz++)
	if (extents().physicalYeeCells().num(xyz) > 1)
		nDim += 1;
	return nDim;
}

Precision::RationalFunction GridDescription::
backgroundPermittivity() const
{
    if (assembly()->backgroundPermittivity() != 0L)
        return assembly()->backgroundPermittivity()->padeCoefficients();
    
    return Precision::RationalFunction(0.0);
}

Precision::RationalFunction GridDescription::
backgroundPermeability() const
{
    if (assembly()->backgroundPermeability() != 0L)
        return assembly()->backgroundPermeability()->padeCoefficients();
    
    return Precision::RationalFunction(0.0);
}

#pragma mark *** Output ***

OutputField::
OutputField() :
    mFieldAsString("no field"),
    mField(kD),
    m_ii(0), m_jj(0),
    mOctant(0)
{
}

OutputField::
OutputField(string fieldName) :
    mFieldAsString(fieldName)
{
    init(fieldName);
}

OutputField::
OutputField(string fieldName, int octant) :
    mFieldAsString(fieldName)
{
    init(fieldName);
    mOctant = octant;
}

void OutputField::
init(string fieldName)
{
    Octant oct;
    
    if (fieldName.length() == 2)
    {
        switch (fieldName[0])
        {
            case 'e': mField = kE; oct = Octant::electric(); mTime = 0.0; break;
            case 'd': mField = kD; oct = Octant::electric(); mTime = 0.0; break;
            case 'h': mField = kH; oct = Octant::magnetic(); mTime = 0.5; break;
            case 'b': mField = kB; oct = Octant::magnetic(); mTime = 0.5; break;
            case 'j': mField = kJ; oct = Octant::electric(); mTime = -0.5; break;
            case 'm': mField = kM; oct = Octant::magnetic(); mTime = 0.0; break;
            default:
                throw(std::logic_error(string("Invalid field ") + fieldName));
                break;
        };
        
        if (fieldName[1] >= 'x' && fieldName[1] <= 'z')
        {
            int xyz = fieldName[1] - 'x';
            mOctant = oct(xyz, xyz);
            m_ii = fieldName[1] - 'x';
            m_jj = m_ii;
        }
        else
            throw(std::logic_error(string("Invalid field ") + fieldName));
    }
    else if (fieldName.length() == 3) // e.g. exy, hzz
    {
        switch (fieldName[0])
        {
            case 'e': mField = kEE; oct = Octant::electric(); mTime = 0.0; break;
            case 'd': mField = kD; oct = Octant::electric(); mTime = 0.0; break;
            case 'h': mField = kHH; oct = Octant::magnetic(); mTime = 0.5; break;
            case 'b': mField = kB; oct = Octant::magnetic(); mTime = 0.5; break;
            default:
                throw(std::logic_error(string("Invalid field ") + fieldName));
                break;
        };
        
        // FIXME: This is really really really kludgey.
        if (mField == kEE || mField == kHH)
        {
            int ii = fieldName[1] - 'x';
            int jj = fieldName[2] - 'x';
            
            if (ii < 0 || ii > 2 || jj < 0 || jj > 2)
                throw(std::logic_error(string("Invalid field ") + fieldName));
            
            mOctant = oct(ii,jj);
            m_ii = ii;
            m_jj = jj;
        }
        else if (mField == kD || mField == kB)
        {
            int ii = fieldName[1] - 'x';
            
            mOctant = 0;
            m_ii = ii;
            m_jj = ii;
        }
    }
    else
        throw(std::logic_error(string("Invalid field ") + fieldName));
}

Vector3d OutputField::
cellOffset() const
{
    return YeeUtilities::halfCellPosition(octant());
}

bool operator==(const OutputField & lhs, const OutputField & rhs)
{
    return (lhs.field() == rhs.field()
        && lhs.ii() == rhs.ii()
        && lhs.jj() == rhs.jj()
        && lhs.octant() == rhs.octant());
}

OutputDescription::
OutputDescription(string fields, string file) :
    mFile(file),
    mRegions(vector<Region>(1,Region())),
    mDurations(vector<Duration>(1,Duration()))
{
    initFields(fields);
}

OutputDescription::
OutputDescription(string fields, string file,
    Region region, Duration duration) :
    mFile(file),
    mRegions(vector<Region>(1,region)),
    mDurations(vector<Duration>(1,duration))
{
    initFields(fields);
}

OutputDescription::
OutputDescription(string fields, string file,
    const std::vector<Region> & regions,
    const std::vector<Duration> & durations) :
    mFile(file),
    mRegions(regions),
    mDurations(durations)
{
    initFields(fields);
}

bool OutputDescription::
includes(const Field & field) const
{
    for (int nn = 0; nn < mOutputFields.size(); nn++)
    {
        if (mOutputFields[nn].field() == field.name() &&
            mOutputFields[nn].ii() == field.xyz() &&
            mOutputFields[nn].jj() == field.xyz())
            return true;
    }
    return false;
}


void OutputDescription::
initFields(string fields)
{
    istringstream istr(fields);
    string field;
    
    while (istr >> field && field != "")
    {
        OutputField outField(field);
        mOutputFields.push_back(outField);
    }
}

#pragma mark *** MaterialOutput ***

MaterialOutputDescription::
MaterialOutputDescription()
{
}

#pragma mark *** GridReport ***

GridReportDescription::
GridReportDescription() :
    mFileName("defaultGridReport"),
    mRegions()
{
    for (int nn = 0; nn < 8; nn++)
        mOctants[nn] = 1;
}

GridReportDescription::
GridReportDescription(string fields, string file, vector<Region> regions) :
    mFileName(file),
    mRegions(regions)
{
    SourceFields tempSrcFields(fields);
    
    for (int xyz = 0; xyz < 3; xyz++)
    {
        mOctants[octantE(xyz)] = tempSrcFields.whichE()[xyz];
        mOctants[octantH(xyz)] = tempSrcFields.whichH()[xyz];
    }
}

GridReportDescription::
GridReportDescription(string file, vector<Region> regions) :
    mFileName(file),
    mRegions(regions)
{
    for (int nn = 0; nn < 8; nn++)
        mOctants[nn] = 1;
}

bool GridReportDescription::
usesOctant(int octant) const
{
    assert(octant >= 0 && octant < 8);
    return mOctants[octant];
}



#pragma mark *** SourceFields ***

SourceFields::
SourceFields() :
    mUsesPolarization(0),
    mPolarization(0.0, 0.0, 0.0),
    mWhichE(0,0,0),
    mWhichH(0,0,0)
{
}

SourceFields::
SourceFields(string fields) :
    mUsesPolarization(0),
    mPolarization(0.0, 0.0, 0.0),
    mWhichE(0,0,0),
    mWhichH(0,0,0)
{
    istringstream istr(fields);
    string field;
    while(istr >> field && field != "")
    {
        if (field == "electric")
            mWhichE = Vector3i(1,1,1);
        else if (field == "magnetic")
            mWhichH = Vector3i(1,1,1);
        else if (field.length() == 2)
        {
            char c1 = field[0];
            char c2 = field[1];
            
            int direction = int(c2 - 'x');
            if (direction < 0 || direction > 2)
                throw(std::logic_error(string("Invalid field ") + field));
            
            if (c1 == 'e')
                mWhichE[direction] = 1;
            else if (c1 == 'h')
                mWhichH[direction] = 1;
            else
                throw(std::logic_error(string("Invalid field ") + field));
        }
        else
            throw(std::logic_error(string("Invalid field ") + field));
    }
}

SourceFields::
SourceFields(string fields, Precision::Vec3 polarization) :
    mUsesPolarization(1),
    mPolarization(polarization),
    mWhichE(0,0,0),
    mWhichH(0,0,0)
{
    if (fields == "electric" || fields == "e")
        mWhichE = Vector3i(1,1,1);
    else if (fields == "magnetic" || fields == "h")
        mWhichH = Vector3i(1,1,1);
    else
        throw(std::logic_error(string("When polarization is provided, the field "
            "must be 'electric' or 'magnetic'; field is ") + fields));
}

#pragma mark *** Source ***

SourceDescription* SourceDescription::
newTimeSource(string timeFile, SourceFields fields, bool isSoft,
    const vector<Region> & regions, const vector<Duration> & durations)
{
    return new SourceDescription(fields, "", timeFile, "", "",
        isSoft, regions, durations);
}

SourceDescription* SourceDescription::
newSpaceTimeSource(string spaceTimeFile, SourceFields fields, bool isSoft,
    const vector<Region> & regions, const vector<Duration> & durations)
{
    return new SourceDescription(fields, "", "", "", spaceTimeFile,
        isSoft, regions, durations);
}

SourceDescription* SourceDescription::
newFormulaSource(string formula, SourceFields fields, bool isSoft,
    const vector<Region> & regions, const vector<Duration> & durations)
{
    return new SourceDescription(fields, formula, "", "", "",
        isSoft, regions, durations);
}

SourceDescription::
SourceDescription(SourceFields fields, string formula, string timeFile,
    string spaceFile, string spaceTimeFile, bool isSoft,
    const vector<Region> & regions,
    const vector<Duration> & durations) :
    mFormula(formula),
    mTimeFile(timeFile),
    mSpaceFile(spaceFile),
    mSpaceTimeFile(spaceTimeFile),
    mFields(fields),
    mRegions(regions),
    mDurations(durations),
    mIsSoft(isSoft)
{
}

#pragma mark *** SourceCurrents ***

SourceCurrents::
SourceCurrents() :
    mUsesPolarization(0),
    mPolarization(0.0, 0.0, 0.0),
    mWhichJ(0,0,0),
    mWhichM(0,0,0),
    mWhichJE(0,0,0),
    mWhichMH(0,0,0)
{
}

SourceCurrents::
SourceCurrents(string fields) :
    mUsesPolarization(0),
    mPolarization(0.0, 0.0, 0.0),
    mWhichJ(0,0,0),
    mWhichM(0,0,0),
    mWhichJE(0,0,0),
    mWhichMH(0,0,0)
{
    istringstream istr(fields);
    string field;
    while(istr >> field && field != "")
    {
        if (field == "electric" || field == "j")
            mWhichJ = Vector3i(1,1,1);
        else if (field == "magnetic" || field == "m")
            mWhichM = Vector3i(1,1,1);
        else if (field == "je")
            mWhichJE = Vector3i(1,1,1);
        else if (field == "mh")
            mWhichMH = Vector3i(1,1,1);
        else if (field.length() == 2) // e.g. "jx", "my"
        {
            char c1 = field[0];
            char c2 = field[1];
            
            int direction = int(c2 - 'x');
            if (direction < 0 || direction > 2)
                throw(std::logic_error(string("Invalid field ") + field));
            
            if (c1 == 'j')
                mWhichJ[direction] = 1;
            else if (c1 == 'm')
                mWhichM[direction] = 1;
            else
                throw(std::logic_error(string("Invalid field ") + field));
        }
        else if (field.length() == 3) // e.g. "jex", "mhy"
        {
            string firstTwo = field.substr(0,2);
            int direction = int(field[2] - 'x');
            
            if (direction < 0 || direction > 2)
                throw(std::logic_error(string("Invalid field ") + field));
            
            if (firstTwo == "je")
                mWhichJE[direction] = 1;
            else if (firstTwo == "mh")
                mWhichMH[direction] = 1;
            else
                throw(std::logic_error(string("Invalid field ") + field));
        }
        else
            throw(std::logic_error(string("Invalid field ") + field));
    }
}

SourceCurrents::
SourceCurrents(string fields, Precision::Vec3 polarization) :
    mUsesPolarization(1),
    mPolarization(polarization),
    mWhichJ(0,0,0),
    mWhichM(0,0,0)
{
    if (fields == "electric")
        mWhichJ = Vector3i(1,1,1);
    else if (fields == "magnetic")
        mWhichM = Vector3i(1,1,1);
    else
        throw(std::logic_error(string("When polarization is provided, the field "
            "must be 'electric' or 'magnetic'; field is ") + fields));
}

#pragma mark *** CurrentSource ***


CurrentSourceDescription* CurrentSourceDescription::
newTimeSource(string timeFile, SourceCurrents currents,
    const vector<Region> & regions, const vector<Duration> & durations)
{
    return new CurrentSourceDescription(currents, "", timeFile, "", "",
        regions, durations);
}

CurrentSourceDescription* CurrentSourceDescription::
newMaskedTimeSource(string timeFile, string spaceFile, SourceCurrents currents,
    const vector<Region> & regions, const vector<Duration> & durations)
{
    return new CurrentSourceDescription(currents, "", timeFile, spaceFile, "",
        regions, durations);
}

CurrentSourceDescription* CurrentSourceDescription::
newSpaceTimeSource(string spaceTimeFile, SourceCurrents currents,
    const vector<Region> & regions, const vector<Duration> & durations)
{
    return new CurrentSourceDescription(currents, "", "", "", spaceTimeFile,
        regions, durations);
}

CurrentSourceDescription* CurrentSourceDescription::
newFormulaSource(string formula, SourceCurrents currents,
    const vector<Region> & regions, const vector<Duration> & durations)
{
    return new CurrentSourceDescription(currents, formula, "", "", "",
        regions, durations);
}

CurrentSourceDescription::
CurrentSourceDescription(SourceCurrents currents, string formula,
    string timeFile, string spaceFile, string spaceTimeFile,
    const vector<Region> & regions, const vector<Duration> & durations)
    :
    mFormula(formula),
    mTimeFile(timeFile),
    mSpaceFile(spaceFile),
    mSpaceTimeFile(spaceTimeFile),
    mCurrents(currents),
    mRegions(regions),
    mDurations(durations)
{
}


#pragma mark *** HuygensSurface ***

HuygensSurfaceDescription::
HuygensSurfaceDescription() :
    mDirection(0,0,0),
    mSymmetries(0,0,0)
{
}

HuygensSurfaceDescription::
HuygensSurfaceDescription(HuygensSurfaceSourceType type) :
    mType(type),
    mDirection(0,0,0),
    mSymmetries(0,0,0)
{
}

void HuygensSurfaceDescription::
setPointers(const Map<string, GridDescPtr> & gridMap)
{
	assert(mType == kLink);
	mSourceGrid = gridMap[mSourceGridName];
    
    if (!mSourceGrid->extents().physicalYeeCells().encloses(mFromYeeCells))
        throw(std::logic_error("Link fromYeeCells must be entirely"
            " contained within sourceGrid's bounds."));
}

void HuygensSurfaceDescription::
omitSide(int sideNum)
{
	Vector3i side = cardinal(sideNum);
	mOmittedSides.insert(side);
}

void HuygensSurfaceDescription::
omitSide(Vector3i direction)
{
    mOmittedSides.insert(direction);
}

void HuygensSurfaceDescription::
becomeLink(GridDescPtr sourceGrid, const Rect3i & sourceYeeCells)
{
    mType = kLink;
    mFromYeeCells = sourceYeeCells;
    mSourceGrid = sourceGrid;
    mSourceGridName = sourceGrid->name();
    
    for (int xyz = 0; xyz < 3; xyz++)
    if (fromYeeCells().num(xyz) != yeeCells().num(xyz) &&
            fromYeeCells().num(xyz) != 1)
    {
        throw(std::logic_error("All dimensions of fromYeeCells must"
            " either be the same as toYeeCells or span the entire"
            " dimension of the source grid, which must be one Yee cell across"
            " in that direction."));
    }
    
//    LOG << "Source half cells " << sourceHalfCells << "\n";
//    LOG << "Dest half cells " << mHalfCells << "\n";
}

bool HuygensSurfaceDescription::
omitsSide(int side) const
{
    return mOmittedSides.count(cardinal(side));
}

Rect3i HuygensSurfaceDescription::
yeeCells(int face, int octant) const
{
    Vector3i direction = cardinal(face);
    Rect3i sideCells = edgeOfRect(yeeToHalf(yeeCells()), face);
    if (face%2 == 0) // low X, low Y, or low Z side
        sideCells.p1 += direction;
    else // high X, high Y or high Z side
        sideCells.p2 += direction;
    
    // sideHalfCells is now a half cell inside the TF region and a half
    // cell in the SF region.  We just need to return the right octant
    // of it.
    Rect3i yeeCells = halfToYee(sideCells, octant);
    
    return yeeCells;
}

HuygensSurfaceDescription* HuygensSurfaceDescription::
newTFSFTimeSource(SourceFields fields, string timeFile, Vector3i direction,
    Rect3i yeeCells, set<Vector3i> omittedSides, bool isTF)
{
	if (!vec_ge(yeeCells.num(), 1))
		throw(std::logic_error("TFSFSource rect has some negative dimensions"));
    if (direction != dominantComponent(direction))
        throw(std::logic_error("direction needs to be an axis-oriented unit vector."));
    
    HuygensSurfaceDescription* hs2 = new HuygensSurfaceDescription(kTFSFSource);
    hs2->mFields = fields;
    hs2->mTimeFile = timeFile;
    hs2->mDirection = direction;
    for (int xyz = 0; xyz < 3; xyz++)
    if (direction[xyz] == 0)
        hs2->mSymmetries[xyz] = 1;
    hs2->mYeeCells = yeeCells;
    hs2->mOmittedSides = omittedSides;
    hs2->mIsTotalField = isTF;
    
//    LOG << "Half cells " << halfCells << "\n";
    
    return hs2;
}

HuygensSurfaceDescription* HuygensSurfaceDescription::
newTFSFFormulaSource(SourceFields fields, string formula, Vector3i direction,
    Rect3i yeeCells, set<Vector3i> omittedSides, bool isTF)
{
	if (!vec_ge(yeeCells.num(), 1))
		throw(std::logic_error("TFSFSource rect has some negative dimensions"));
    if (direction != dominantComponent(direction))
        throw(std::logic_error("direction needs to be an axis-oriented unit vector."));
    LOG << "Warning: not validating formula.\n";
    
    HuygensSurfaceDescription* hs2 = new HuygensSurfaceDescription(kTFSFSource);
    hs2->mFields = fields;
    hs2->mFormula = formula;
    hs2->mDirection = direction;
    for (int xyz = 0; xyz < 3; xyz++)
    if (direction[xyz] == 0)
        hs2->mSymmetries[xyz] = 1;
    hs2->mDuration = Duration();
    hs2->mYeeCells = yeeCells;
    hs2->mOmittedSides = omittedSides;
    hs2->mIsTotalField = isTF;
    
    //LOG << "Half cells " << halfCells << "\n";
    
    return hs2;
}

HuygensSurfaceDescription* HuygensSurfaceDescription::
newCustomTFSFSource(string file, Vector3i symmetries, Rect3i yeeCells,
    Duration duration, set<Vector3i> omittedSides, bool isTF)
{
	if (!vec_ge(yeeCells.num(), 1))
		throw(std::logic_error("TFSFSource rect has some negative dimensions"));
    
    HuygensSurfaceDescription* hs2 =
        new HuygensSurfaceDescription(kCustomTFSFSource);
    hs2->mFile = file;
    hs2->mDuration = duration;
    hs2->mSymmetries = symmetries;
    hs2->mYeeCells = yeeCells;
    hs2->mOmittedSides = omittedSides;
    hs2->mIsTotalField = isTF;
    
    return hs2;
}

HuygensSurfaceDescription* HuygensSurfaceDescription::
newLink(string sourceGrid, Rect3i fromYeeCells, Rect3i toYeeCells,
    set<Vector3i> omittedSides, bool isTF)
{
	if (!vec_ge(fromYeeCells.num(), 1))
		throw(std::logic_error("Link from rect has some negative dimensions"));
	if (!vec_ge(toYeeCells.num(), 1))
		throw(std::logic_error("Link to rect has some negative dimensions"));
    for (int xyz = 0; xyz < 3; xyz++)
    if (fromYeeCells.num(xyz) != toYeeCells.num(xyz) &&
            fromYeeCells.num(xyz) != 1)
        throw(std::logic_error("All dimensions of fromYeeCells must"
            " either be the same as toYeeCells or span the entire"
            " dimension of the source grid, which must be one Yee cell across"
            " in that direction."));
    
    HuygensSurfaceDescription* hs2 = new HuygensSurfaceDescription(kLink);
    hs2->mSourceGridName = sourceGrid;
    hs2->mFromYeeCells = fromYeeCells;
    hs2->mYeeCells = toYeeCells;
    hs2->mOmittedSides = omittedSides;
    hs2->mIsTotalField = isTF;
    
    return hs2;
}

#pragma mark *** Material ***

MaterialDescription::
MaterialDescription(int ID, string name, P::RationalFunction padeCoeffs) :
    mID(ID),
	mName(name),
	mPadeZ(padeCoeffs)
{
}

ostream &
operator<<(ostream & out, const MaterialDescription & mat)
{
	out << mat.name();
	return out;
}

#pragma mark *** Assembly ***

AssemblyDescription::
AssemblyDescription(const vector<InstructionPtr> & recipe,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices,
    const std::vector<std::vector<SimpleMesh::Triangle> > & inMeshes) :
	mInstructions(recipe),
    mControlVertices(controlVertices),
    mMeshes(inMeshes)
{
    int numBackgrounds = 0;
	for (int nn = 0; nn < mInstructions.size(); nn++)
    {
        if (instructions().at(nn)->type() == kBackgroundType)
            numBackgrounds++;
    }
    
    if (numBackgrounds > 1)
        throw(std::runtime_error("Grid has more than one background material."));
}

void AssemblyDescription::
setPointers(const Map<string, GridDescPtr> & gridMap)
{
	for (unsigned int nn = 0; nn < mInstructions.size(); nn++)
	{
		InstructionPtr ii = mInstructions[nn];
		switch (ii->type())
		{
			case kCopyFromType:
				((CopyFrom&)*ii).setPointers(gridMap);
				break;
			case kExtrudeType:
				// don't need to do anything here since there are no pointers
				//((Extrude&)*ii).setPointers(gridMap);
				break;
		}
	}
}

const MaterialDescPtr AssemblyDescription::
backgroundPermittivity() const
{
    for (unsigned int nn = 0; nn < mInstructions.size(); nn++)
    {
        InstructionPtr ii = mInstructions[nn];
        if (ii->type() == kBackgroundType)
            return ((const Background&)*ii).permittivity();
    }
    return MaterialDescPtr(0L);
}

const MaterialDescPtr AssemblyDescription::
backgroundPermeability() const
{
    for (unsigned int nn = 0; nn < mInstructions.size(); nn++)
    {
        InstructionPtr ii = mInstructions[nn];
        if (ii->type() == kBackgroundType)
            return ((const Background&)*ii).permeability();
    }
    return MaterialDescPtr(0L);
}

Instruction::
Instruction(InstructionType inType) :
	mType(inType)
{
}


NewMesh::
NewMesh(const vector<Vector3i> & faces,
    const vector<Vector3i> & controlFaces,
    MaterialDescPtr permittivity,
    MaterialDescPtr permeability) :
    Instruction(kNewMeshType),
    mFaces(faces),
    mControlFaces(controlFaces),
    mPermittivity(permittivity),
    mPermeability(permeability)
{
}

string NewMesh::
permittivityName() const
{
    if (mPermittivity != 0L)
        return mPermittivity->name();
    else
        return "";
}

string NewMesh::
permeabilityName() const
{
    if (mPermeability != 0L)
        return mPermeability->name();
    else
        return "";
}



Mesh::
Mesh(const vector<Vector3d> & vertices,
    const vector<int> & controlVertexIds,
    const vector<Vector3b> & vertexFreeDirections,
    const vector<Vector3i> & faces,
    MaterialDescPtr permittivity,
    MaterialDescPtr permeability) :
    Instruction(kMeshType),
    mVertices(vertices),
    mControlVertexIds(controlVertexIds),
    mFreeDirections(vertexFreeDirections),
    mFaces(faces),
    mPermittivity(permittivity),
    mPermeability(permeability)
{
}

Rect3d Mesh::
yeeBounds() const
{
    if (mVertices.size() == 0)
        return Rect3d();
    
    Rect3d bounds(mVertices[0], mVertices[0]);
    for (int vv = 1; vv < mVertices.size(); vv++)
    {
        bounds.p1 = vec_min(bounds.p1, mVertices[vv]);
        bounds.p2 = vec_max(bounds.p2, mVertices[vv]);
    }
    return bounds;
}

string Mesh::
permittivityName() const
{
    if (mPermittivity != 0L)
        return mPermittivity->name();
    else
        return "";
}

string Mesh::
permeabilityName() const
{
    if (mPermeability != 0L)
        return mPermeability->name();
    else
        return "";
}


CopyFrom::
CopyFrom(Rect3i halfCellSourceRegion, Rect3i halfCellDestRegion,
	string gridName) :
	Instruction(kCopyFromType),
	mSourceRect(halfCellSourceRegion),
	mDestRect(halfCellDestRegion),
	mGridName(gridName)
{
	//cerr << "Warning: minimal validation done for CopyFrom().\n";
	// Easy validation: no inside-out rects
	if (!vec_ge(mSourceRect.size(), 0))
		throw(std::logic_error("Some source rect dimensions are negative"));
	if (!vec_ge(mDestRect.size(), 0))
		throw(std::logic_error("Some dest rect dimensions are negative"));
	
	// All copyFrom dimensions must equal copyTo or be 0.
	
	for (int xyz = 0; xyz < 3; xyz++)
	if (mSourceRect.size(xyz) != 0)
	{
		if (mSourceRect.size(xyz) != mDestRect.size(xyz))
		throw(std::logic_error("Error: copy from region must be same size as "
			"copy to region or have size 0"));
		if (mSourceRect.p1[xyz]%2 != mDestRect.p1[xyz]%2)
		throw(std::logic_error("Error: copy from region must start on same half-cell"
			" octant as copy to region or have size 0 (both should start on"
			" even indices or on odd indices)"));
		// it's sufficient to check p1 and not p2 because we already know
		// that the dimensions are equal.
	}
}

CopyFrom::
CopyFrom(Rect3i halfCellSourceRegion, Rect3i halfCellDestRegion,
	const GridDescPtr grid) :
	Instruction(kCopyFromType),
	mSourceRect(halfCellSourceRegion),
	mDestRect(halfCellDestRegion),
	mGridName(grid->name()),
	mGrid(grid)
{
	// Easy validation: no inside-out rects
	if (!vec_ge(mSourceRect.size(), 0))
		throw(std::logic_error("Some source rect dimensions are negative"));
	if (!vec_ge(mDestRect.size(), 0))
		throw(std::logic_error("Some dest rect dimensions are negative"));
	
	// All copyFrom dimensions must equal copyTo or be 0.
	
	for (int xyz = 0; xyz < 3; xyz++)
	if (mSourceRect.size(xyz) != 0)
	{
		if (mSourceRect.size(xyz) != mDestRect.size(xyz))
		throw(std::logic_error("Error: copy from region must be same size as "
			"copy to region or have size 0"));
		if (mSourceRect.p1[xyz]%2 != mDestRect.p1[xyz]%2)
		throw(std::logic_error("Error: copy from region must start on same half-cell"
			" octant as copy to region or have size 0 (both should start on"
			" even indices or on odd indices)"));
		// it's sufficient to check p1 and not p2 because we already know
		// that the dimensions are equal.
	}
}

void CopyFrom::
setPointers(const Map<string, GridDescPtr> & gridMap)
{
	mGrid = gridMap[mGridName];
}

Extrude::
Extrude(Rect3i halfCellExtrudeFrom, Rect3i halfCellExtrudeTo) :
	Instruction(kExtrudeType),
	mExtrudeFrom(halfCellExtrudeFrom),
	mExtrudeTo(halfCellExtrudeTo)
{
	// Easy validation: no inside-out rects
	if (!vec_ge(mExtrudeFrom.size(), 0))
		throw(std::logic_error("ExtrudeFrom dimensions are negative"));
	if (!vec_ge(mExtrudeTo.size(), 0))
		throw(std::logic_error("ExtrudeTo dimensions are negative"));
	if (!mExtrudeTo.encloses(mExtrudeFrom))
		throw(std::logic_error("ExtrudeTo does not enclose ExtrudeFrom"));
}



Background::
Background(MaterialDescPtr permittivity,
    MaterialDescPtr permeability) :
	Instruction(kBackgroundType),
	mPermittivity(permittivity),
    mPermeability(permeability)
{
}

string Background::
permittivityName() const
{
    if (mPermittivity != 0L)
        return mPermittivity->name();
    else
        return "";
}

string Background::
permeabilityName() const
{
    if (mPermeability != 0L)
        return mPermeability->name();
    else
        return "";
}




