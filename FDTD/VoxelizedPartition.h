/*
 *  VoxelizedPartition.h
 *  TROGDOR
 *
 *  Created by Paul Hansen on 2/9/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _VOXELIZEDPARTITION_
#define _VOXELIZEDPARTITION_

#include "Pointer.h"
#include "geometry.h"
#include "SimulationDescription.h"
#include "rle/SupportRegion3.h"
#include "rle/DynamicRLE3.h"
#include "rle/IndexArray3.h"
#include "SimpleMesh/SimpleMesh.h"

#include "Extents.h"
#include "PMLParameters.h"
#include "ConstitutiveTensorField.h"
#include "TensorGridConstants.h"
#include "ElectromagneticFieldIndices.h"
#include "PMLField.h"

#include <vector>
#include <string>

class SimulationDescription;
class HuygensSurfaceDescription;

class FillFactorReport;
class FillFactorSensitivityReport;
class OrientationReport;
class OrientationSensitivityReport;

class VoxelizedPartition
{
public:
    /**
     *  partitionYeeCells   intersection of MPI partition with grid Yee cells
     **/
	VoxelizedPartition(GridDescPtr gridDesc,
        const PhysicalExtents & gridExtents,
        Rect3i partitionYeeCells,
        const std::vector<Precision::RationalFunction> & permittivities,
        const std::vector<Precision::RationalFunction> & permeabilities,
        const std::vector<std::vector<SimpleMesh::Triangle> > & permittivityMesh,
        const std::vector<std::vector<SimpleMesh::Triangle> > & permeabilityMesh,
        const std::vector<SimpleMesh::ControlVertex> & controlVertices);
    
    VoxelizedPartition(GridDescPtr gridDesc,
        const VoxelizedPartition & parentGrid,
        Rect3i copyYeeCells,
        const PhysicalExtents & gridExtents,
        Rect3i partitionYeeCells);
    
//    void copyVoxels();
	
    // Determine whether the current electromagnetic mode (e.g. 3d or kTM_2d)
    // makes use of a given field.
    bool modeUsesField(const Field & f) const;
    
    GridDescPtr gridDescription() const { return mGridDescription; }
    const NodeExtents & extents() const { return mExtents; }
    const ElectromagneticFieldIndices & fieldIndices() const
        { return mFieldIndices; }
    const MaterialIndices & materialIndices() const
        { return mDispersiveMaterialIndices; }
    const PMLIndices & pmlIndices() const
        { return mPMLFieldIndices; }
    
    ConstitutiveTensorField & inversePermittivity()
        { return mInversePermittivity; }
    const ConstitutiveTensorField & inversePermittivity() const
        { return mInversePermittivity; }
    ConstitutiveTensorField & inversePermeability()
        { return mInversePermeability; }
    const ConstitutiveTensorField & inversePermeability() const
        { return mInversePermeability; }
    
    RLE::DynamicRLE3<Precision::RationalFunction> &
    inversePermittivity(int i, int j) { return mInversePermittivity(i,j); }
    RLE::DynamicRLE3<Precision::RationalFunction> &
    inversePermeability(int i, int j) { return mInversePermeability(i,j); }
    const RLE::DynamicRLE3<Precision::RationalFunction> &
    inversePermittivity(int i, int j) const { return mInversePermittivity(i,j); }
    const RLE::DynamicRLE3<Precision::RationalFunction> &
    inversePermeability(int i, int j) const { return mInversePermeability(i,j); }
    
    const PMLField & pmlE() const { return mPMLParametersE; }
    const PMLField & pmlH() const { return mPMLParametersH; }
    
    
    const RLE::DynamicRLE3<PMLParameters> &
    pmlE(int fieldXYZ, int absorbXYZ) const
    {
        return mPMLParametersE.along(absorbXYZ)[fieldXYZ];
    }
    
    const RLE::DynamicRLE3<PMLParameters> &
    pmlH(int fieldXYZ, int absorbXYZ) const
    {
        return mPMLParametersH.along(absorbXYZ)[fieldXYZ];
    }
    
    Vector3b boundarySymmetry(Rect3i yeeCells) const
    {
        std::cerr << "Warning: not actually checking symmetry.\n";
        return Vector3b(1,1,1);
    }
    
    const RLE::SupportRegion3 & sensitiveCells(int ii, int jj) const
    {
        return mSensitiveCells[ii][jj];
    }
    
    Matrix3i projectionMatrix() const
    {
        return Matrix3i::diagonal(Vector3i(extents().dimensions()));
    }
    
    enum GridType
    {
        kUserGrid = 1,
        kAutoGrid = 2
    };
    int gridType() const { return mGridType; }
private:

    void initPermittivity(
        const std::vector<Precision::RationalFunction> & materials,
        const std::vector<std::vector<SimpleMesh::Triangle> > & mesh,
        const std::vector<SimpleMesh::ControlVertex> & controlVertices
        );
    void initPermeability(
        const std::vector<Precision::RationalFunction> & materials,
        const std::vector<std::vector<SimpleMesh::Triangle> > & mesh,
        const std::vector<SimpleMesh::ControlVertex> & controlVertices
        );
    
    void initVoxels(OrientedVoxels & voxels,
        const std::vector<std::vector<SimpleMesh::Triangle> > & polyhedra,
        const std::vector<RLE::DynamicRLE3<double> > & halfCellFillFactors,
        const std::vector<std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > >&
            halfCellFillJacobians,
        int octant,
        FillFactorReport & ffReport, FillFactorSensitivityReport & ffsReport,
        OrientationReport & oReport, OrientationSensitivityReport & osReport);
    
    void addForwardBoundaryOutputs();
    void addAdjointBoundaryOutputs();
    void indexFields();
    
    void calculateCurrentIndices(const GridDescription & desc);
    
    RLE::SupportRegion3 mSensitiveCells[3][3];
    
    GridDescPtr mGridDescription;
    NodeExtents mExtents;
    
    ConstitutiveTensorField mInversePermittivity;
    ConstitutiveTensorField mInversePermeability;
    PMLField mPMLParametersE;
    PMLField mPMLParametersH;
    
    // TODO: Put ConstitutiveTensorField and PML constants here.
    // TODO: Put a support() function in ConstitutiveTensorField.
    ElectromagneticFieldIndices mFieldIndices;
    MaterialIndices mDispersiveMaterialIndices;
    PMLIndices mPMLFieldIndices;
	
	friend std::ostream & operator<< (std::ostream & out,
		const VoxelizedPartition & grid);
    
    int mGridType;
};
typedef Pointer<VoxelizedPartition> VoxelizedPartitionPtr;

std::ostream & operator<< (std::ostream & out, const VoxelizedPartition & grid);

#endif
