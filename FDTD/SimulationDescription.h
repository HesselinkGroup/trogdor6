/*
 *  SimulationDescription.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 1/31/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _SIMULATIONDESCRIPTION_
#define _SIMULATIONDESCRIPTION_

#include "SimulationDescriptionPredeclarations.h"
#include "SimpleMesh/SimpleMesh.h"
#include "PrecisionGeometry.h"
#include "Pointer.h"
#include "Extents.h"
#include "FieldEnum.h"

#include "Map.h"
#include "PrecisionRationalFunction.h"
#include <climits>
#include <vector>
#include <string>
#include <set>

class XMLParameterFile;


enum ElectromagneticMode
{
    k1d,
    kTE2d,
    kTM2d,
    k2d,
    k3d
};

class SimulationDescription
{
	friend class XMLParameterFile;
public:
	SimulationDescription(const XMLParameterFile & file);
	
	void setGrids(const std::vector<GridDescPtr> & grids) { mGrids = grids; }
	void setMaterials(const std::vector<MaterialDescPtr> & materials);
	void setGridPointers();
	
    void setElectromagneticMode(const ElectromagneticMode & mode) { mElectromagneticMode = mode; }
	void setDiscretization(Precision::Vec3 dxyz, Precision::Float dt);
	void setNumTimesteps(int numT);
	
	const std::vector<GridDescPtr> & grids() const { return mGrids; }
	std::map<std::string, MaterialDescPtr> materials() const { return mMaterials; }
    std::map<std::string, Precision::RationalFunction> materialPadeCoefficients() const;
	
    const std::string & version() const { return mVersionString; }
    
    ElectromagneticMode electromagneticMode() const { return mElectromagneticMode; }
	Precision::Float dt() const { return m_dt; }
	Precision::Vec3 dxyz() const { return m_dxyz; }
	int numTimesteps() const { return mNumTimesteps; }
	
private:
    std::string mVersionString;
	std::vector<GridDescPtr> mGrids;
	std::map<std::string, MaterialDescPtr> mMaterials;
	
    ElectromagneticMode mElectromagneticMode;
	Precision::Float m_dt;
	Precision::Vec3 m_dxyz;
	long mNumTimesteps;
};

class GridDescription
{
public:
	GridDescription(std::string name, ElectromagneticMode mode, Rect3i yeeCells,
        Rect3i nonPMLHalf, Precision::Vec3 mDxyz, Precision::Float mDt,
        long numTimesteps);
	
    int id() const { return mID; }
    
	// Mutators
	void setOutputs(const std::vector<OutputDescPtr> & outputs)
        { mOutputs = outputs; }
	void setGridReports(const std::vector<GridReportDescPtr> & gridReports);
	void setSources(const std::vector<SourceDescPtr> & sources)
        { mSources = sources; }
    void setCurrentSources(const std::vector<CurrentSourceDescPtr> & currents)
        { mCurrentSources = currents; }
	void setHuygensSurfaces(
		const std::vector<HuygensSurfaceDescPtr> & surfaces);
	void setAssembly(AssemblyDescPtr assembly) {
		mAssembly = assembly; }
    void setPMLParams(const Map<std::string, std::string> & p);
    
	void setPointers(const Map<std::string, GridDescPtr> & gridMap);
    
    void setOrigin(Precision::Vec3 origin) { mOrigin = origin; }
	
	// Accessors
	const std::string & name() const { return mName; }
    const PhysicalExtents & extents() const { return mExtents; }
    ElectromagneticMode electromagneticMode() const { return mElectromagneticMode; }
    const Map<std::string, std::string> & pmlParams() const
        { return mPMLParams; }
	int numDimensions() const;
    Precision::Vec3 dxyz() const { return mDxyz; }
    Precision::Float dt() const { return mDt; }
    long numTimesteps() const { return mNumTimesteps; }
    Precision::Vec3 origin() const { return mOrigin; }
	
	const std::vector<OutputDescPtr> & outputs() const { return mOutputs; }
    const std::vector<GridReportDescPtr> & gridReports() const
        { return mGridReports; }
	const std::vector<SourceDescPtr> & sources() const { return mSources; }
    const std::vector<CurrentSourceDescPtr> & currentSources() const
        { return mCurrentSources; }
	const std::vector<HuygensSurfaceDescPtr> & huygensSurfaces() const
		{ return mHuygensSurfaces; }
	std::vector<HuygensSurfaceDescPtr> & huygensSurfaces()
		{ return mHuygensSurfaces; }
    /**
     * Access all the triangulated meshes from the XML file.
     */
	const AssemblyDescPtr assembly() const { return mAssembly; }
    
    Precision::RationalFunction backgroundPermittivity() const;
    Precision::RationalFunction backgroundPermeability() const;
    
    void pushOutput(OutputDescPtr out) { mOutputs.push_back(out); }
	
private:
	std::string mName;
    int mID;
    static int sNextGridID;
	
    PhysicalExtents mExtents;
    
    ElectromagneticMode mElectromagneticMode;
    Precision::Vec3 mDxyz;
    Precision::Float mDt;
    long mNumTimesteps;
    Precision::Vec3 mOrigin;
    
    Map<std::string, std::string> mPMLParams;
	
	std::vector<OutputDescPtr> mOutputs;
    std::vector<GridReportDescPtr> mGridReports;
	std::vector<SourceDescPtr> mSources;
    std::vector<CurrentSourceDescPtr> mCurrentSources;
    std::vector<HuygensSurfaceDescPtr> mHuygensSurfaces;
	AssemblyDescPtr mAssembly;
};

class Duration
{
public:
    Duration() :
        mFirst(0),
        mLast(INT_MAX),
        mPeriod(1),
        mHasInterval(false) {}
    Duration(int firstTimestep, int lastTimestep, int period = 1) :
        mFirst(firstTimestep),
        mLast(lastTimestep),
        mPeriod(period),
        mHasInterval(false) {}
    Duration(int firstTimestep, int lastTimestep, Vector2d interval, 
        int period = 1) :
        mFirst(firstTimestep),
        mLast(lastTimestep),
        mPeriod(period),
        mInterval(interval),
        mHasInterval(true)
    {}
    
    static Duration from(int firstTimestep, int period = 1)
    {
        return Duration(firstTimestep, INT_MAX, period);
    }
    static Duration periodic(int period)
    {
        return Duration(0, INT_MAX, period);
    }
    
    long first() const { return mFirst; }
    long last() const { return mLast; }
    long period() const { return mPeriod; }
    
    bool hasInterval() const { return mHasInterval; }
    Vector2d interval() const
    {
        if (false == mHasInterval)
            throw(std::runtime_error("Duration has undefined interval"));
        return mInterval;
    }
    
    void setFirst(int first) { mFirst = first; }
    void setLast(int last) { mLast = last; }
    void setPeriod(int period) { mPeriod = period; }
    
private:
    long mFirst;
    long mLast;
    long mPeriod;
    Vector2d mInterval;
    bool mHasInterval;
};
//std::ostream & operator<<(std::ostream & str, const Duration & dur);

class Region
{
public:
    Region() :
        mYeeCells(-INT_MAX, -INT_MAX, -INT_MAX, INT_MAX, INT_MAX, INT_MAX),
        mBounds(),
        mHasBounds(false),
        mStride(1,1,1) {}
    Region(const Rect3i & yeeCells) :
        mYeeCells(yeeCells),
        mBounds(),
        mHasBounds(false),
        mStride(1,1,1) {}
    Region(const Rect3i & yeeCells, Vector3i stride) :
        mYeeCells(yeeCells),
        mBounds(),
        mHasBounds(false),
        mStride(stride) {}
    Region(const Rect3i & yeeCells, const Rect3d & bounds, Vector3i stride) :
        mYeeCells(yeeCells),
        mBounds(bounds),
        mHasBounds(true),
        mStride(stride) {}
    
    void setYeeCells(const Rect3i & rect) { mYeeCells = rect; }
    void setBounds(const Rect3d & rect) { mBounds = rect; mHasBounds = true; }
    void setStride(const Vector3i & stride) { mStride = stride; }
    const Rect3i & yeeCells() const { return mYeeCells; }
    const Rect3d & bounds() const
    {
        if (false == hasBounds())
            throw(std::runtime_error("Region has undefined bounds"));
        return mBounds;
    }
    bool hasBounds() const { return mHasBounds; }
    const Vector3i & stride() const { return mStride; }
private:
    Rect3i mYeeCells;
    Rect3d mBounds;
    bool mHasBounds;
    Vector3i mStride;
};
//std::ostream & operator<<(std::ostream & str, const Region & reg);

class OutputField
{
public:
    OutputField();
    OutputField(std::string fieldName);
    OutputField(std::string fieldName, int octant); // interpolate to octant
    
    void init(std::string fieldName);
    
    const std::string & fieldAsString() const { return mFieldAsString; }
    FieldName field() const { return mField; }
    int ii() const { return m_ii; }
    int jj() const { return m_jj; }
    int octant() const { return mOctant; }
    Vector3d cellOffset() const;
    double timeOffset() const { return mTime; }
    
private:
    std::string mFieldAsString;
    FieldName mField;
    int m_ii, m_jj; // field components, vectorial (ii) or tensorial (ii,jj)
    int mOctant;
    double mTime;
};
bool operator==(const OutputField & lhs, const OutputField & rhs);


class OutputDescription
{
public:
    OutputDescription(std::string fields, std::string file);
    OutputDescription(std::string fields, std::string file,
        Region region, Duration duration = Duration());
    OutputDescription(std::string fields, std::string file,
        const std::vector<Region> & regions,
        const std::vector<Duration> & durations);
    
    const std::string & file() const { return mFile; }
    const std::vector<OutputField> & fields() const { return mOutputFields; }
    const std::vector<Region> & regions() const { return mRegions; }
    const std::vector<Duration> & durations() const { return mDurations; }
    bool includes(const Field & field) const;
    
private:
    void initFields(std::string fields);
    
    std::string mFile;
    std::vector<OutputField> mOutputFields;
    std::vector<Region> mRegions;
    std::vector<Duration> mDurations;
};

class MaterialOutputDescription
{
public:
    MaterialOutputDescription();
private:
};

class GridReportDescription
{
public:
    GridReportDescription();
    GridReportDescription(std::string fields, std::string file,
        std::vector<Region> regions);
    GridReportDescription(std::string file, std::vector<Region> regions);
    
    std::string fileName() const { return mFileName; }
    std::vector<Region> regions() const { return mRegions; }
    std::vector<Region> & regions() { return mRegions; }
    bool usesOctant(int octant) const;
private:
    std::string mFileName;
    std::vector<Region> mRegions;
    bool mOctants[8];
};

class SourceFields
{
public:
    SourceFields();
    SourceFields(std::string fields);
    SourceFields(std::string fields, Precision::Vec3 polarization);
    
    bool usesPolarization() const { return mUsesPolarization; }
    Precision::Vec3 polarization() const
        { assert(mUsesPolarization); return mPolarization; }
    Vector3i whichE() const { return mWhichE; }
    Vector3i whichH() const { return mWhichH; }
private:
    bool mUsesPolarization;
    Precision::Vec3 mPolarization;
    Vector3i mWhichE;
    Vector3i mWhichH;
};

class SourceDescription
{
public:
    static SourceDescription* newTimeSource(std::string timeFile,
        SourceFields fields, bool isSoft, const std::vector<Region> & regions,
            const std::vector<Duration> & durations);
    static SourceDescription* newSpaceTimeSource(std::string spaceTimeFile,
        SourceFields fields, bool isSoft, const std::vector<Region> & regions,
        const std::vector<Duration> & durations);
    static SourceDescription* newFormulaSource(std::string formula,
        SourceFields fields, bool isSoft, const std::vector<Region> & regions,
        const std::vector<Duration> & durations);
    
    SourceDescription(SourceFields fields, std::string formula,
        std::string timeFile, std::string spaceFile, std::string spaceTimeFile,
        bool isSoft, const std::vector<Region> & regions,
        const std::vector<Duration> & durations);
    
    const std::string & formula() const { return mFormula; }
    const std::string & timeFile() const { return mTimeFile; }
    const std::string & spaceFile() const { return mSpaceFile; }
    const std::string & spaceTimeFile() const { return mSpaceTimeFile; }
    const SourceFields & sourceFields() const { return mFields; }
    
    bool isHardSource() const { return !mIsSoft; }
    bool isSoftSource() const { return mIsSoft; }
    bool hasMask() const { return (mSpaceFile != ""); }
    bool isSpaceVarying() const { return (mSpaceTimeFile != ""); }
    
    const std::vector<Region> & regions() const { return mRegions; }
    const std::vector<Duration> & durations() const { return mDurations; }
    
private:
    std::string mFormula;
    std::string mTimeFile;
    std::string mSpaceFile;
    std::string mSpaceTimeFile;
    SourceFields mFields;
    std::vector<Region> mRegions;
    std::vector<Duration> mDurations;
    
    bool mIsSoft;
};

class SourceCurrents
{
public:
    SourceCurrents();
    SourceCurrents(std::string fields);
    SourceCurrents(std::string fields, Precision::Vec3 polarization);
    
    bool usesPolarization() const { return mUsesPolarization; }
    Precision::Vec3 polarization() const
        { assert(mUsesPolarization); return mPolarization; }
    Vector3i whichJ() const { return mWhichJ; }
    Vector3i whichM() const { return mWhichM; }
    Vector3i whichJE() const { return mWhichJE; }
    Vector3i whichMH() const { return mWhichMH; }
private:
    bool mUsesPolarization;
    Precision::Vec3 mPolarization;
    Vector3i mWhichJ;
    Vector3i mWhichM;
    Vector3i mWhichJE;
    Vector3i mWhichMH;
};

class CurrentSourceDescription
{
public:
    static CurrentSourceDescription* newTimeSource(std::string timeFile,
        SourceCurrents currents, const std::vector<Region> & regions,
            const std::vector<Duration> & durations);
    static CurrentSourceDescription* newMaskedTimeSource(std::string timeFile,
        std::string spaceFile, SourceCurrents currents,
        const std::vector<Region> & regions,
        const std::vector<Duration> & durations);
    static CurrentSourceDescription* newSpaceTimeSource(
        std::string spaceTimeFile, SourceCurrents currents, 
        const std::vector<Region> & regions,
        const std::vector<Duration> & durations);
    static CurrentSourceDescription* newFormulaSource(std::string formula,
        SourceCurrents fields, const std::vector<Region> & regions,
        const std::vector<Duration> & durations);
    
    bool hasMask() const
        { return (mSpaceFile != ""); }
    bool isSpaceVarying() const
        { return (mSpaceTimeFile != ""); }
    const std::string & formula() const { return mFormula; }
    const std::string & timeFile() const { return mTimeFile; }
    const std::string & spaceFile() const { return mSpaceFile; }
    const std::string & spaceTimeFile() const { return mSpaceTimeFile; }
    const SourceCurrents & sourceCurrents() const { return mCurrents; }
    
    const std::vector<Region> & regions() const { return mRegions; }
    const std::vector<Duration> & durations() const { return mDurations; }
private:
    CurrentSourceDescription(SourceCurrents currents, std::string formula,
        std::string timeFile, std::string spaceFile, std::string spaceTimeFile,
        const std::vector<Region> & regions,
        const std::vector<Duration> & durations);
    
    std::string mFormula;
    std::string mTimeFile;
    std::string mSpaceFile;
    std::string mSpaceTimeFile;
    SourceCurrents mCurrents;
    std::vector<Region> mRegions;
    std::vector<Duration> mDurations;
};

enum HuygensSurfaceSourceType
{
	kLink,
	kTFSFSource,
	kCustomTFSFSource
};

class HuygensSurfaceDescription
{
public:
    HuygensSurfaceDescription();
    
    HuygensSurfaceDescription(HuygensSurfaceSourceType type);
    
    static HuygensSurfaceDescription*
    newTFSFTimeSource(SourceFields fields,
        std::string timeFile, Vector3i direction, Rect3i yeeCells,
        std::set<Vector3i> omittedSides, bool isTF = 1);
        
    static HuygensSurfaceDescription*
    newTFSFFormulaSource(SourceFields fields,
        std::string formula, Vector3i direction, Rect3i yeeCells,
        std::set<Vector3i> omittedSides, bool isTF = 1);
    
    static HuygensSurfaceDescription*
    newCustomTFSFSource(std::string file,
        Vector3i symmetries, Rect3i yeeCells, Duration duration,
        std::set<Vector3i> omittedSides, bool isTF = 1);
    
    static HuygensSurfaceDescription*
    newLink(std::string sourceGrid,
        Rect3i fromYeeCells, Rect3i toYeeCells,
        std::set<Vector3i> omittedSides, bool isTF = 1);
    
	// modifiers
	void setPointers(const Map<std::string, GridDescPtr> & gridMap);
	void omitSide(int nSide);
    void omitSide(Vector3i dir);
    
    void becomeLink(GridDescPtr sourceGrid,
        const Rect3i & sourceYeeCells);
    
    void yeeCells(Rect3i cells) { mYeeCells = cells; }
    
    // Common accessors
    HuygensSurfaceSourceType type() const { return mType; }
    const Rect3i & yeeCells() const { return mYeeCells; }
    const std::set<Vector3i> & omittedSides() const { return mOmittedSides; }
    bool omitsSide(int side) const;
    bool isTotalField() const { return mIsTotalField; }
    bool isScatteredField() const { return !mIsTotalField; }
    /**
     * Each face of the Huygens surface (six faces) has three field components
     * on the TF side and three on the SF side.  Check the results of this
     * function to find out whether the returned region is TF or SF.
     */
    Rect3i yeeCells(int face, int octant) const;
    
    // TFSFSource accessors
    SourceFields sourceFields() const
        { assert(mType == kTFSFSource); return mFields; }
    std::string timeFile() const
        { assert(mType == kTFSFSource); return mTimeFile; }
    std::string formula() const
        { assert(mType == kTFSFSource); return mFormula; }
    Vector3i direction() const
        { assert(mType == kTFSFSource); return mDirection; }
    
    // Custom or ordinary TFSFSource accessor
    Vector3i symmetries() const
        { assert(mType == kTFSFSource || mType == kCustomTFSFSource);
          return mSymmetries; }
    Duration duration() const
        { assert(mType == kCustomTFSFSource || mType == kTFSFSource);
          return mDuration; }
    
    // Custom source accessors
    std::string file() const
        { assert(mType == kCustomTFSFSource); return mFile; }
    
    // Link accessors
    std::string sourceGridName() const
        { assert(mType == kLink); return mSourceGridName; }
    GridDescPtr sourceGrid() const
        { assert(mType == kLink); return mSourceGrid; }
    Rect3i fromYeeCells() const { return mFromYeeCells; }
    
private:
    // Common data
    HuygensSurfaceSourceType mType;
    Rect3i mYeeCells;
    std::set<Vector3i> mOmittedSides;
    bool mIsTotalField;
    
    // Source data
    SourceFields mFields;
    std::string mTimeFile;
    std::string mFormula;
    Vector3i mDirection;
    
     // TFSF or custom TFSF source
    Vector3i mSymmetries;
    Duration mDuration;
    
    // Custom source data
    std::string mFile;
    
    // Link data
    std::string mSourceGridName;
    GridDescPtr mSourceGrid;
    Rect3i mFromYeeCells;
};

class MaterialDescription
{
public:
	MaterialDescription(int ID, std::string name, Precision::RationalFunction padeCoeffs);
	
    int id() const { return mID; }
	std::string name() const { return mName; }
    const Precision::RationalFunction & padeCoefficients() const { return mPadeZ; }
    
	friend std::ostream & operator<<(std::ostream & out,
		const MaterialDescription & mat);
private:
    int mID;
	std::string mName;
    Precision::RationalFunction mPadeZ;
};
std::ostream & operator<<(std::ostream & out, const MaterialDescription & mat);

#pragma mark *** Assembly things ***

class Instruction;
typedef Pointer<Instruction> InstructionPtr;

/**
 * Just a list of InstructionPtrs i.e. triangulated meshes from the XML.
 */
class AssemblyDescription
{
public:
	AssemblyDescription(const std::vector<InstructionPtr> & recipe,
        const std::vector<SimpleMesh::ControlVertex> & controlVertices,
        const std::vector<std::vector<SimpleMesh::Triangle> > & inMeshes);
	
	void setInstructions(const std::vector<InstructionPtr> & instructions);
	void setPointers(const Map<std::string, GridDescPtr> & gridMap);
	
	const std::vector<InstructionPtr> & instructions() const
		{ return mInstructions; }
    const std::vector<SimpleMesh::ControlVertex> & controlVertices() const
        { return mControlVertices; }
    const std::vector<std::vector<SimpleMesh::Triangle> > & meshes() const
        { return mMeshes; }
    
    const MaterialDescPtr backgroundPermittivity() const;
    const MaterialDescPtr backgroundPermeability() const;
	
private:
	std::vector<InstructionPtr> mInstructions;
    std::vector<SimpleMesh::ControlVertex> mControlVertices;
    std::vector<std::vector<SimpleMesh::Triangle> > mMeshes;
};

enum InstructionType
{
    kNewMeshType,
    kMeshType,
	kCopyFromType,
	kExtrudeType,
    kBackgroundType,
    kPermityType
};

// The purpose of the Instruction base class is to provide some by-hand
// runtime type information without having to turn on RTTI, which I think
// is kind of evil.  (Well, C++ is evil.)
//
// A little RTTI will permit outside methods to determine which Instruction
// we are dealing with and handle each differently, while still receiving
// the Instructions in one vector or list.
class Instruction
{
public:
	Instruction(InstructionType inType);
	InstructionType type() const { return mType; }
    virtual ~Instruction() {}
    
protected:
	InstructionType mType;
};
typedef Pointer<Instruction> InstructionPtr;


class NewMesh : public Instruction
{
public:
    NewMesh(const std::vector<Vector3i> & faces,
        const std::vector<Vector3i> & controlFaces,
        MaterialDescPtr permittivity, MaterialDescPtr permeability);
    
    const MaterialDescPtr & permittivity() const { return mPermittivity; }
    const MaterialDescPtr & permeability() const { return mPermeability; }
    std::string permittivityName() const;
    std::string permeabilityName() const;
    
    const std::vector<Vector3i> & faces() const { return mFaces; }
    const std::vector<Vector3i> & controlFaces() const { return mControlFaces; }
private:
    std::vector<Vector3i> mFaces;
    std::vector<Vector3i> mControlFaces;
    MaterialDescPtr mPermittivity;
    MaterialDescPtr mPermeability;
};

class Mesh : public Instruction
{
public:
    Mesh(const std::vector<Vector3d> & vertices,
        const std::vector<int> & controlVertexIds,
        const std::vector<Vector3b> & vertexFreeDirections,
        const std::vector<Vector3i> & faces,
        MaterialDescPtr permittivity, MaterialDescPtr permeability);
    
    Rect3d yeeBounds() const;
    const MaterialDescPtr & permittivity() const { return mPermittivity; }
    const MaterialDescPtr & permeability() const { return mPermeability; }
    std::string permittivityName() const;
    std::string permeabilityName() const;
    
    const std::vector<Vector3d> & vertices() const { return mVertices; }
    const std::vector<int> & controlVertexIds() const { return mControlVertexIds; }
    const std::vector<Vector3b> & freeDirections() const
        { return mFreeDirections; }
    const std::vector<Vector3i> & faces() const { return mFaces; }
private:
    std::vector<Vector3d> mVertices;
    std::vector<int> mControlVertexIds;
    std::vector<Vector3b> mFreeDirections;
    std::vector<Vector3i> mFaces;
    MaterialDescPtr mPermittivity;
    MaterialDescPtr mPermeability;
};



class CopyFrom : public Instruction
{
public:
	CopyFrom(Rect3i halfCellSourceRegion, Rect3i halfCellDestRegion,
		std::string gridName);
	CopyFrom(Rect3i halfCellSourceRegion, Rect3i halfCellDestRegion,
		const GridDescPtr grid);
	
	void setPointers(const Map<std::string, GridDescPtr> & gridMap);
	
	const Rect3i & sourceHalfRect() const { return mSourceRect; }
	const Rect3i & destHalfRect() const { return mDestRect; }
	const std::string & gridName() const { return mGridName; }
	const GridDescPtr & grid() const { return mGrid; }
	
private:
	Rect3i mSourceRect;
	Rect3i mDestRect;
	std::string mGridName;
	GridDescPtr mGrid;
};

class Extrude : public Instruction
{
public:
	Extrude(Rect3i halfCellExtrudeFrom, Rect3i halfCellExtrudeTo);
	
	const Rect3i & extrudeFrom() const { return mExtrudeFrom; }
	const Rect3i & extrudeTo() const { return mExtrudeTo; }
private:
	Rect3i mExtrudeFrom;
	Rect3i mExtrudeTo;
};

class Background : public Instruction
{
public:
    Background(MaterialDescPtr permittivity, MaterialDescPtr permeability);
    
	const MaterialDescPtr & permittivity() const { return mPermittivity; }
	const MaterialDescPtr & permeability() const { return mPermeability; }
    std::string permittivityName() const;
    std::string permeabilityName() const;
private:
	MaterialDescPtr mPermittivity;
	MaterialDescPtr mPermeability;
};

class Perm_ity : public Instruction
{
public:
    Perm_ity(const std::string & type, const std::string & materialName, int solidId) :
        Instruction(kPermityType),
        mPermityType(type),
        mMaterialName(materialName),
        mSolidId(solidId)
    {
        if (mPermityType != "Permittivity" && mPermityType != "Permeability")
        {
            throw std::runtime_error("Type must be Permittivity or Permeability");
        }
    }
    
    int getSolidId() const { return mSolidId; }
    const std::string & getMaterialName() const { return mMaterialName; }
    const std::string & getType() const { return mPermityType; }
    
private:
    std::string mPermityType; // "permittivity" or "permeability"
    std::string mMaterialName;
    int mSolidId;
};




#endif
