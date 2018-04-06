/*
 *  VoxelizedPartition.cpp
 *  TROGDOR
 *
 *  Created by Paul Hansen on 2/9/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "VoxelizedPartition.h"

#include "Log.h"
#include "YeeUtilities.h"
#include "OrientedVoxels.h"
#include "FillFactorReport.h"
#include "OrientationReport.h"
#include "VoxelizePolyhedra.h"
#include "RLEOperations.h"
#include "Support.h"
#include "WriteSomeRLE.h"
#include "UserPreferences.h"
#include "OutputDirectory.h"

#include <algorithm>
#include <sstream>
#include <map>

using namespace std;
using namespace YeeUtilities;
using namespace RLE;

#pragma mark *** VoxelizedPartition ***

VoxelizedPartition::
VoxelizedPartition(GridDescPtr gridDesc,
    const PhysicalExtents & gridExtents,
    Rect3i partitionYeeCells,
    const std::vector<Precision::RationalFunction> & permittivities,
    const std::vector<Precision::RationalFunction> & permeabilities,
    const std::vector<std::vector<SimpleMesh::Triangle> > & permittivityMesh,
    const std::vector<std::vector<SimpleMesh::Triangle> > & permeabilityMesh,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices) :
    mGridDescription(gridDesc),
    mExtents(gridExtents, partitionYeeCells),
    mInversePermittivity(Octant::electric(), mExtents.allocatedYeeCells()),
    mInversePermeability(Octant::magnetic(), mExtents.allocatedYeeCells()),
    mPMLParametersE(Octant::electric()),
    mPMLParametersH(Octant::magnetic()),
    mGridType(kUserGrid)
{
    LOG << "Node:\n";
    LOGMORE << "\tphysical Yee cells = " << mExtents.physicalYeeCells() << "\n";
    LOGMORE << "\tcalc half cells = " << mExtents.calcHalfCells() << "\n";
    LOGMORE << "\tghost half cells = " << mExtents.ghostHalfCells() << "\n";
    LOGMORE << "\tnon-PML half cells = " << mExtents.nonPMLHalfCells() << "\n";
    LOGMORE << "\tconstants half cells = " << mExtents.constantsHalfCells()
        << "\n";
    LOGMORE << "\tallocated Yee cells = " << mExtents.allocatedYeeCells()
        << "\n";
    LOGMORE << "\tboundary conditions: \n";
    for (int nn = 0; nn < 6; nn++)
    {
        LOGMORE << "\t\t" << nn << ": ";
        switch(mExtents.nodeBoundaries().at(nn))
        {
            case kPECBoundary: LOGMORE << "PEC\n"; break;
            case kPMCBoundary: LOGMORE << "PMC\n"; break;
            case kPeriodicBoundary: LOGMORE << "Periodic\n"; break;
            case kMPIBoundary: LOGMORE << "MPI\n"; break;
            case kPMLBoundary: LOGMORE << "PML\n"; break;
            case kTranslationSymmetricBoundary: LOGMORE << "Symmetric\n"; break;
            default: LOGMORE << "UNKNOWN\n"; break;
        }
    }
    
    LOG << "Voxelize...\n";
    
    if (permittivityMesh.size() == 0 &&
        gridDescription()->backgroundPermittivity() == Precision::Float(0.0))
        throw(std::logic_error("No permittivities!"));
    if (permeabilityMesh.size() == 0 &&
        gridDescription()->backgroundPermeability() == Precision::Float(0.0))
        throw(std::logic_error("No permeabilities!"));
    
    // NB: epsilon and mu share all control vertices.
    initPermittivity(permittivities, permittivityMesh, controlVertices);
    initPermeability(permeabilities, permeabilityMesh, controlVertices);
    
    // All I have done here is turn the Assembly map<string,Poly> into two
    // vectors, and get the list of control vertices, and sent these to
    // initPermittivity().
    
    inversePermittivity().extrude(
        intersection(mExtents.constantsHalfCells(), mExtents.nonPMLHalfCells()),
        mExtents.constantsHalfCells());
    inversePermeability().extrude(
        intersection(mExtents.constantsHalfCells(), mExtents.nonPMLHalfCells()),
        mExtents.constantsHalfCells());
    // FIXME: extrude() will lose information on sensitive cells.  Extrude too?
    // Actually, the best idea might be to extrude the fill factors and
    // orientations (extrude inside OrientedVoxels).  This will fix everything
    // else.
        
    mPMLParametersE.init(gridExtents, gridDesc->dxyz(), gridDesc->pmlParams());
    mPMLParametersH.init(gridExtents, gridDesc->dxyz(), gridDesc->pmlParams());
    
    if (gridType() == VoxelizedPartition::kUserGrid)
    {
        LOG << "Add boundary outputs...\n";
        if (UserPreferences::defines("adjoint"))
            addAdjointBoundaryOutputs();
        else
            addForwardBoundaryOutputs();
    }
    
    LOG << "Index...\n";
    indexFields();
    
    LOG << "Memory:\n";
    LOGMORE << "\tEstimated field data (E, H, D, B) (vacuum) "
        << mExtents.allocatedYeeCells().count() * 12.0 * sizeof(Precision::Float)
        / 1024.0 << " kB\n";
    LOGMORE << "\tGrid physical constants: "
        << (mInversePermittivity.bytes() + mInversePermeability.bytes())/1024.0
        << " kB\n";
    LOGMORE << "\tPML constants: "
        << (mPMLParametersE.bytes() + mPMLParametersH.bytes()) / 1024.0
        << " kB\n";
    LOGMORE << "\tMaterial indices "
        << mDispersiveMaterialIndices.bytes()/1024.0 << " kB\n";
    LOGMORE << "\tField indices: " << mFieldIndices.bytes()/1024.0 << " kB\n";
    LOGMORE << "\tPML indices " << mPMLFieldIndices.bytes()/1024.0 << " kB\n";
    LOGMORE << "\tHuygens surfaces: UNIMPLEMENTED\n";
    LOGMORE << "\tBuffered outputs: UNIMPLEMENTED\n";
}

VoxelizedPartition::
VoxelizedPartition(GridDescPtr gridDesc,
    const VoxelizedPartition & parentGrid,
    Rect3i copyYeeCells,
    const PhysicalExtents & gridExtents,
    Rect3i partitionYeeCells) :
    mGridDescription(gridDesc),
    mExtents(gridDesc->extents(), partitionYeeCells),
    mInversePermittivity(Octant::electric(), mExtents.allocatedYeeCells()),
    mInversePermeability(Octant::magnetic(), mExtents.allocatedYeeCells()),
    mPMLParametersE(Octant::electric()),
    mPMLParametersH(Octant::magnetic()),
    mGridType(kAutoGrid)
{
    LOG << "Node:\n";
    LOGMORE << "\tphysical Yee cells = " << mExtents.physicalYeeCells() << "\n";
    LOGMORE << "\tcalc half cells = " << mExtents.calcHalfCells() << "\n";
    LOGMORE << "\tghost half cells = " << mExtents.ghostHalfCells() << "\n";
    LOGMORE << "\tnon-PML half cells = " << mExtents.nonPMLHalfCells() << "\n";
    LOGMORE << "\tconstants half cells = " << mExtents.constantsHalfCells()
        << "\n";
    LOGMORE << "\tallocated Yee cells = " << mExtents.allocatedYeeCells()
        << "\n";
    LOGMORE << "Copy Yee cells = " << copyYeeCells << "\n";
    
    Mat3i collapser(Mat3i::diagonal(Vector3i(gridDesc->extents().dimensions())));
    Rect3i pasteYeeCells(collapser * copyYeeCells);
    LOGMORE << "Paste Yee cells = " << pasteYeeCells << "\n";
    
    // Copy permittivity omitting sensitivities for now.
    inversePermittivity().copy(parentGrid.inversePermittivity(), copyYeeCells, pasteYeeCells);
    inversePermeability().copy(parentGrid.inversePermeability(), copyYeeCells, pasteYeeCells);
    
//    cout << "Permittivity:\n";
//    for (int xyz = 0; xyz < 3; xyz++)
//        cout << inversePermittivity()(xyz,xyz) << "\n";
//    cout << "Permeability:\n";
//    for (int xyz = 0; xyz < 3; xyz++)
//        cout << inversePermeability()(xyz,xyz) << "\n";
    
//    Rect3i r1 = intersection(mExtents.constantsHalfCells(),
//        yeeToHalf(pasteYeeCells));
//    Rect3i r2 = mExtents.constantsHalfCells();
    
    inversePermittivity().extrude(
        intersection(mExtents.constantsHalfCells(), yeeToHalf(pasteYeeCells)),
        mExtents.constantsHalfCells());
    inversePermeability().extrude(
        intersection(mExtents.constantsHalfCells(), yeeToHalf(pasteYeeCells)),
        mExtents.constantsHalfCells());
    
//    cout << "Permittivity:\n";
//    for (int xyz = 0; xyz < 3; xyz++)
//        cout << inversePermittivity()(xyz,xyz) << "\n";
//    cout << "Permeability:\n";
//    for (int xyz = 0; xyz < 3; xyz++)
//        cout << inversePermeability()(xyz,xyz) << "\n";
    
    // FIXME: extrude() will lose information on sensitive cells.  Extrude too?
    // Actually, the best idea might be to extrude the fill factors and
    // orientations (extrude inside OrientedVoxels).  This will fix everything
    // else.
        
    mPMLParametersE.init(gridExtents, gridDesc->dxyz(), gridDesc->pmlParams());
    mPMLParametersH.init(gridExtents, gridDesc->dxyz(), gridDesc->pmlParams());
    
    LOG << "Index...\n";
    mDispersiveMaterialIndices.index(inversePermittivity(),
        inversePermeability());
    mFieldIndices.index(mExtents, inversePermittivity(),
        inversePermeability());
    mPMLFieldIndices.index(mExtents, inversePermittivity(),
        inversePermeability(), pmlE(), pmlH());
    
    // What does this do?
    // It iterates through the grid denizens and determines which cells of the
    // grid need to store J or M (currents).  It should also check all the
    // outputs and determine the extra locations on the grid where I need to
    // store D and B.
    calculateCurrentIndices(*gridDesc);
    
    /*
    mConstitutiveParameters.copyConstitutives(parentGrid,
        copyYeeCells);
    
    mConstitutiveParameters.extrude(yeeToHalf(copyYeeCells),
        mExtents.constantsHalfCells());
    mConstitutiveParameters.setPML(gridDesc->extents(), gridDesc->dxyz(),
        'PML'+'ARGS');
    
    LOG << "Index...\n";
    mDispersiveMaterialIndices.index(mConstitutiveParameters);
    mFieldIndices.index(mExtents, mConstitutiveParameters);
    mPMLFieldIndices.index(mExtents, mConstitutiveParameters);
    
    calculateCurrentIndices(*gridDesc);
    */
}

// Determine whether the current electromagnetic mode (e.g. 3d or kTM_2d)
// makes use of a given field.
bool VoxelizedPartition::modeUsesField(const Field & f) const
{
    ElectromagneticMode mode = gridDescription()->electromagneticMode();
    
    switch (mode)
    {
        case k1d:
            throw(std::runtime_error("1d mode not supported"));
        case k2d:
            return true;
        case kTE2d:
            if (f.isElectric())
            {
                return gridDescription()->extents().transverseAxis() == f.xyz();
            }
            else
            {
                return gridDescription()->extents().transverseAxis() != f.xyz();
            }
            break;
        case kTM2d:
            if (f.isElectric())
            {
                return gridDescription()->extents().transverseAxis() != f.xyz();
            }
            else
            {
                return gridDescription()->extents().transverseAxis() == f.xyz();
            }
            break;
        case k3d:
            return true;
            break;
    };
}

//void VoxelizedPartition::
//createHuygensSurfaces(const GridDescPtr & gridDescription,
//    const Map<GridDescPtr, VoxelizedPartitionPtr> & grids)
//{   
//    const vector<HuygensSurfaceDescPtr> & surfaces = 
//        gridDescription->huygensSurfaces();
//    
//    for (unsigned int nn = 0; nn < surfaces.size(); nn++)
//    {
//        ostringstream huygensSurfaceName;
//        huygensSurfaceName << gridDescription->name() << " HS " << nn;
//        mHuygensSurfaces.push_back(
//            HuygensSurfaceFactory::newHuygensSurface(
//                huygensSurfaceName.str(), *this, grids,
//                surfaces[nn]));
//    }
//    
//    paintFromHuygensSurfaces(*gridDescription);
//}
//
//void VoxelizedPartition::
//writeDataRequests() const
//{
//    // 1.  Current sources
//    // 2.  Custom TFSF sources (HuygensSurfaces)
//    
//    //LOG << "Writing data requests." << endl;
//    
//    for (int nn = 0; nn < mSetupCurrentSources.size(); nn++)
//    if (mSetupCurrentSources[nn]->description()->hasMask() ||
//        mSetupCurrentSources[nn]->description()->isSpaceVarying())
//    {
//        // if this current source will require scheduled data from a file
//        
//        ostringstream str;
//        str << "currentreq_" << nn << ".m";
//        IODescriptionFile::write(
//            str.str(),
//            mSetupCurrentSources[nn]->description(),
//            *this,
//            mSetupCurrentSources[nn]->getRectsJ(),
//            mSetupCurrentSources[nn]->getRectsK());
//    }
//    
//    for (int nn = 0; nn < mHuygensSurfaces.size(); nn++)
//    if (mHuygensSurfaces[nn]->description()->type() == kCustomTFSFSource)
//    {
//        ostringstream str;
//        str << "tfsfreq_" << nn << ".m";
//        IODescriptionFile::write(
//            str.str(),
//            *mHuygensSurfaces[nn],
//            *this,
//            mHuygensRegionSymmetries[nn]);
//    }
//}
//


// Return true for zero matrix or for any matrix with a 1.0 on the diagonal.
bool sIsAxisAlignedSurface(const Matrix3d & matrix)
{
    const double TOL = 1e-4;
    
    static const Matrix3d sMX(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    static const Matrix3d sMY(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    static const Matrix3d sMZ(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    
    if (normInf(matrix - sMX) < TOL ||
        normInf(matrix - sMY) < TOL ||
        normInf(matrix - sMZ) < TOL ||
        normInf(matrix) < TOL)
    {
        return true;
    }
    
    /*
    Vector3d diag = matrix.diagonal();
    
    if (normInf(diag - Vector3d::unit(0)) < TOL ||
        normInf(diag - Vector3d::unit(1)) < TOL ||
        normInf(diag - Vector3d::unit(2)) < TOL ||
        normInf(diag) == 0)
    {
        return true;
    }
    */
    
    return false;
}

void VoxelizedPartition::
initPermittivity(
    const std::vector<Precision::RationalFunction> & materials,
    const std::vector<std::vector<SimpleMesh::Triangle> > & meshTriangles,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices)
{
    Rect3d realBounds(
        Vector3d(gridDescription()->dxyz()) * Vector3d(extents().constantsHalfCells().p1)/2.0,
        Vector3d(gridDescription()->dxyz()) * Vector3d(extents().constantsHalfCells().p2 + 1)/2.0);
    VoxelizePolyhedra vp(extents().constantsHalfCells(), realBounds);
    vector<DynamicRLE3<double> > halfCellFillFactors;
    vector<map<unsigned int, DynamicRLE3<Vector3d> > > halfCellFillJacobians;
    vp.fillFactors(meshTriangles, halfCellFillFactors);
    vp.fillFactorJacobians(meshTriangles, halfCellFillJacobians);
    
    ofstream DepsilonStream;
    if (UserPreferences::defines("sensitivity"))
        DepsilonStream.open(OutputDirectory::path("Depsilon").c_str(), ios::binary);
    
    OrientationReport orientationReport;
    OrientationSensitivityReport orientationSensitivityReport;
    if (UserPreferences::defines("saveOrientation"))
        orientationReport.open("orientations", extents().physicalYeeCells());
    if (UserPreferences::defines("saveOrientationSensitivity"))
        orientationSensitivityReport.open("orientationSensitivity.m",
            extents().physicalYeeCells());
    
    FillFactorReport fillFactorReport;
    FillFactorSensitivityReport fillFactorSensitivityReport;
    if (UserPreferences::defines("saveFillFactors"))
        fillFactorReport.open("fillFactors", extents().physicalYeeCells());
    if (UserPreferences::defines("saveFillFactorSensitivity"))
        fillFactorSensitivityReport.open("fillFactorSensitivity.m",
            extents().physicalYeeCells());
    
    for (int jj = 0; jj < 3; jj++)
    for (int ii = 0; ii < 3; ii++)
    {
        int octant = inversePermittivity().tensorOctant(ii, jj);
        OrientedVoxels voxels;
        
        LOG << "Voxelizing octant " << octant << ".\n";
        
        if (ii == jj)
        {
            initVoxels(voxels, meshTriangles,
                halfCellFillFactors, halfCellFillJacobians, octant,
                fillFactorReport, fillFactorSensitivityReport,
                orientationReport, orientationSensitivityReport);
            
            if (UserPreferences::defines("excludeCorners"))
            {
                cerr << "Excluding corners!\n";
                SupportRegion3 axisAlignedSurfaces;
                transform(voxels.orientations(), axisAlignedSurfaces, sIsAxisAlignedSurface);
                mSensitiveCells[ii][jj] = voxels.sensitiveCells(ii, jj) * axisAlignedSurfaces;
            }
            else
            {
                mSensitiveCells[ii][jj] = voxels.sensitiveCells(ii, jj);
            }
            
//            cerr << "ii jj:\n" << mSensitiveCells[ii][jj] << "\n";
            
            TensorGridConstants tensors(ii, jj, inversePermittivity(ii, jj),
                voxels.voxelBounds(),
                gridDescription()->backgroundPermittivity());
            // TODO: Right here I need to have fill factors and orientations
            // with the right dimensionality.
            tensors.calculateConstitutives(materials, voxels.fillFactors(),
                voxels.orientations(ii, jj));
            
            if (UserPreferences::defines("sensitivity"))
            {
                // This just needs the number of CVs and the free directions,
                // and if I did away with free directions then I would not
                // need this controlVertices array at all.
                tensors.writeSensitivities(materials, controlVertices,
                    voxels.fillFactors(), voxels.orientations(ii,jj),
                    voxels.fillFactorJacobians(),
                    voxels.orientationJacobians(ii,jj),
                    mSensitiveCells[ii][jj],
                    DepsilonStream);
            }
        }
        else if (UserPreferences::defines("fullTensor"))
        {
            initVoxels(voxels, meshTriangles,
                halfCellFillFactors, halfCellFillJacobians, octant,
                fillFactorReport, fillFactorSensitivityReport,
                orientationReport, orientationSensitivityReport);
            
            mSensitiveCells[ii][jj] = voxels.sensitiveCells(ii,jj);
            
            TensorGridConstants tensors(ii, jj, inversePermittivity(ii,jj),
                voxels.voxelBounds(),
                gridDescription()->backgroundPermittivity());
            tensors.calculateConstitutives(materials, voxels.fillFactors(),
                voxels.orientations(ii,jj));
            
            if (UserPreferences::defines("sensitivity"))
                throw(std::logic_error("Can't do sensitivity calculation on"
                    " full tensor."));
            
//            if (UserPreferences::defines("sensitivity"))
//            {
//                tensors.writeSensitivities(materials, controlVertices,
//                    voxels.fillFactors(), voxels.orientations(ii,jj),
//                    voxels.fillFactorJacobians(),
//                    voxels.orientationJacobians(ii,jj),
//                    mSensitiveCells[ii][jj],
//                    Depsilon);
//            }
        }
        
    }
    
    if (UserPreferences::defines("saveOrientation"))
        orientationReport.close();
    if (UserPreferences::defines("saveOrientationSensitivity"))
        orientationSensitivityReport.close();
    if (UserPreferences::defines("saveFillFactors"))
        fillFactorReport.close();
    if (UserPreferences::defines("saveFillFactorSensitivity"))
        fillFactorSensitivityReport.close();
    
    if (UserPreferences::defines("sensitivity"))
        DepsilonStream.close();
    
    if (UserPreferences::defines("savePermittivity"))
    {
        mInversePermittivity.writeBinary("permittivity",
            extents().physicalYeeCells());
    }
}

void VoxelizedPartition::
initPermeability(
    const std::vector<Precision::RationalFunction> & materials,
    const std::vector<std::vector<SimpleMesh::Triangle> > & meshTriangles,
    const std::vector<SimpleMesh::ControlVertex> & controlVertices)
{
    Rect3d realBounds(
        Vector3d(gridDescription()->dxyz()) * Vector3d(extents().constantsHalfCells().p1)/2.0,
        Vector3d(gridDescription()->dxyz()) * Vector3d(extents().constantsHalfCells().p2 + 1)/2.0);
    VoxelizePolyhedra vp(extents().constantsHalfCells(), realBounds);
    vector<DynamicRLE3<double> > halfCellFillFactors;
    //vector<map<unsigned int, DynamicRLE3<Vector3d> > > halfCellFillJacobians;
    vp.fillFactors(meshTriangles, halfCellFillFactors);
    //vp.fillFactorJacobians(polyhedra, controlVertices, halfCellFillJacobians);
    
    // Calculate diagonal elements of pemeability tensor
    for (int xyz = 0; xyz < 3; xyz++)
    {
        int octant = inversePermeability().tensorOctant(xyz, xyz);
        Rect3i voxelBounds = halfToYee(extents().constantsHalfCells(), octant);
        Rect3d realBounds(YeeUtilities::realBounds(voxelBounds, octant));
        realBounds.p1 *= Vector3d(gridDescription()->dxyz());
        realBounds.p2 *= Vector3d(gridDescription()->dxyz());
        Rect3d clipBounds(
            Vector3d(gridDescription()->dxyz()) * Vector3d(extents().physicalYeeCells().p1),
            Vector3d(gridDescription()->dxyz()) * Vector3d(extents().physicalYeeCells().p2) + 1.0);
        
        OrientedVoxels voxels(voxelBounds, realBounds, clipBounds);
        
        voxels.downsampleFillFactors(halfCellFillFactors, octant);
//        voxels.downsampleFillJacobians(halfCellFillJacobians, octantH(xyz));
        voxels.initOrientations(meshTriangles); // and Jacobians
//        mSensitiveCells[octantH(xyz)] = voxels.sensitiveCells(xyz, xyz);
        
        TensorGridConstants tensors(xyz, xyz, inversePermeability(xyz, xyz),
            voxels.voxelBounds(), gridDescription()->backgroundPermeability());
        tensors.calculateConstitutives(materials, voxels.fillFactors(),
            voxels.orientations(xyz, xyz));
//        tensors.calculateSensitivities(materials, controlVertices,
//            voxels.fillFactors(), voxels.orientations(xyz,xyz),
//            voxels.fillFactorJacobians(), voxels.orientationJacobians(xyz,xyz),
//            mSensitiveCells[octantH(xyz)]);
    }
}

void VoxelizedPartition::
initVoxels(OrientedVoxels & voxels,
    const std::vector<std::vector<SimpleMesh::Triangle> > & meshTriangles,
    const std::vector<RLE::DynamicRLE3<double> > & halfCellFillFactors,
    const std::vector<std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > >&
        halfCellFillJacobians,
    int octant,
    FillFactorReport & ffReport, FillFactorSensitivityReport & ffsReport,
    OrientationReport & oReport, OrientationSensitivityReport & osReport)
{
    Rect3i voxelBounds = halfToYee(extents().constantsHalfCells(), octant);
    Rect3d realBounds(YeeUtilities::realBounds(voxelBounds, octant));
    realBounds.p1 *= Vector3d(gridDescription()->dxyz());
    realBounds.p2 *= Vector3d(gridDescription()->dxyz());
    Rect3d clipBounds(
        Vector3d(gridDescription()->dxyz()) * Vector3d(extents().physicalYeeCells().p1),
        Vector3d(gridDescription()->dxyz()) * Vector3d(extents().physicalYeeCells().p2) + 1.0);
    
    voxels = OrientedVoxels(voxelBounds, realBounds, clipBounds);
    
    if (UserPreferences::defines("oSmoothing"))
        voxels.setOrientationSmoothing(UserPreferences::valueAs<int>("oSmoothing"));
    if (UserPreferences::defines("fSmoothing"))
        voxels.setFillFactorSmoothing(UserPreferences::valueAs<int>("fSmoothing"));
    
    voxels.downsampleFillFactors(halfCellFillFactors, octant);
    voxels.downsampleFillJacobians(halfCellFillJacobians, octant);
    voxels.initOrientations(meshTriangles); // and Jacobians
    
//    LOG << "Orientation (oct " << octant << "):\n" << voxels.orientations() << "\n";
    
    if (UserPreferences::defines("saveOrientation"))
        oReport.write(voxels);
    if (UserPreferences::defines("saveOrientationSensitivity"))
    {
        for (int jj = 0; jj < 3; jj++)
        for (int ii = 0; ii < 3; ii++)
            osReport.write(voxels, ii, jj, octant);
    }
    if (UserPreferences::defines("saveFillFactors"))
        ffReport.write(voxels);
    if (UserPreferences::defines("saveFillFactorSensitivity"))
        ffsReport.write(voxels, octant);
}

void VoxelizedPartition::
addAdjointBoundaryOutputs()
{
    // E, D
    vector<Duration> durations(1, Duration(0, gridDescription()->numTimesteps()-1));
    
    if (false == UserPreferences::defines("sensitivity"))
        return;
    
    for (int ii = 0; ii < 3; ii++)
    {
        int jj = ii;
    //for (int jj = 0; jj < 3; jj++)
    if (sensitiveCells(ii,jj).numRuns() > 0)
    {
//        int octant = octantE(ii,jj);
        vector<Region> regions(sensitiveCells(ii,jj).numRuns());
        SupportRegion3::ConstIterator itr;
        int rr;
        for (rr = 0, itr = sensitiveCells(ii,jj).begin();
            itr != sensitiveCells(ii,jj).end();
            rr++, itr.nextMarkedRun())
        {
            Rect3i yeeCells;
            sensitiveCells(ii,jj).cartesianCoordinates(
                itr.runStart(), yeeCells.p1[0], yeeCells.p1[1], yeeCells.p1[2]);
            sensitiveCells(ii,jj).cartesianCoordinates(
                itr.runEnd(), yeeCells.p2[0], yeeCells.p2[1], yeeCells.p2[2]);
            regions[rr] = Region(yeeCells);
        }
        
        ostringstream fields;
        fields << "e" << char('x'+ii) << char('x'+jj);
//        if (ii == jj)
//            fields << " d" << char('x'+ii);
        string filename;
        filename = string("boundary_adjoint_e") + char('x'+ii) + char('x'+jj);
        OutputDescPtr out(new OutputDescription(fields.str(), filename,
            regions, durations));
        gridDescription()->pushOutput(out);
    }
    }
}

void VoxelizedPartition::
addForwardBoundaryOutputs()
{
    // E, D
    vector<Duration> durations(1, Duration(0, gridDescription()->numTimesteps()-1));
    
    if (false == UserPreferences::defines("sensitivity"))
        return;
    
    for (int ii = 0; ii < 3; ii++)
    {
        int jj = ii;
    //for (int jj = 0; jj < 3; jj++)
    //if (voxelGrid.sensitiveCells(ii,jj).numRuns() > 0)
    {
//        int octant = octantE(ii,jj);
        vector<Region> regions(sensitiveCells(ii,jj).numRuns());
        SupportRegion3::ConstIterator itr;
        int rr;
        for (rr = 0, itr = sensitiveCells(ii,jj).begin();
            itr != sensitiveCells(ii,jj).end();
            rr++, itr.nextMarkedRun())
        {
            Rect3i yeeCells;
            sensitiveCells(ii,jj).cartesianCoordinates(
                itr.runStart(), yeeCells.p1[0], yeeCells.p1[1], yeeCells.p1[2]);
            sensitiveCells(ii,jj).cartesianCoordinates(
                itr.runEnd(), yeeCells.p2[0], yeeCells.p2[1], yeeCells.p2[2]);
            regions[rr] = Region(yeeCells);
        }
        
        ostringstream fields;
        fields << "e" << char('x'+ii) << " d" << char('x'+ii);
        //fields << "e" << char('x'+ii) << char('x'+jj);
        //fields << " d" << char('x'+ii) << char('x'+jj); // an interpolated field
        string filename;
        filename = string("boundary_de") + char('x'+ii) + char('x'+jj);
        OutputDescPtr out(new OutputDescription(fields.str(), filename,
            regions, durations));
        gridDescription()->pushOutput(out);
    }
    }
}

void VoxelizedPartition::
indexFields()
{
    mDispersiveMaterialIndices.index(inversePermittivity(),
        inversePermeability());
    mFieldIndices.index(mExtents, inversePermittivity(),
        inversePermeability());
    mPMLFieldIndices.index(mExtents, inversePermittivity(),
        inversePermeability(), pmlE(), pmlH());
    
    // Make sure we store D, B, E, and H wherever we need outputs as well.
    
    SupportRegion3 outSupp;
    
    if (UserPreferences::defines("adjoint"))
    {
        for (int xyz = 0; xyz < 3; xyz++)
        for (int nn = 0; nn < gridDescription()->outputs().size(); nn++)
        {
            const OutputDescription & output = *(gridDescription()->outputs()[nn]);
            
            if (output.includes(Field(kE, xyz)))
            {
                outSupp = support(output.regions(), extents().physicalYeeCells().num());
                
                if (mFieldIndices.e(xyz).numRuns() > 0)
                {
                    mFieldIndices.e(xyz, outSupp +
                        SupportRegion3(mFieldIndices.e(xyz)));
                }
                else
                    mFieldIndices.e(xyz, outSupp);
            }
            
            if (output.includes(Field(kH, xyz)))
            {
                outSupp = support(output.regions(), extents().physicalYeeCells().num());
                
                if (mFieldIndices.h(xyz).numRuns() > 0)
                {
                    mFieldIndices.h(xyz, outSupp +
                        SupportRegion3(mFieldIndices.h(xyz)));
                }
                else
                    mFieldIndices.h(xyz, outSupp);
            }
        }
    }
    else
    {
        for (int xyz = 0; xyz < 3; xyz++)
        for (int nn = 0; nn < gridDescription()->outputs().size(); nn++)
        {
            const OutputDescription & output = *(gridDescription()->outputs()[nn]);
            
            if (output.includes(Field(kD, xyz)))
            {
                outSupp = support(output.regions(), extents().physicalYeeCells().num());
                
                if (mFieldIndices.d(xyz).numRuns() > 0)
                {
                    mFieldIndices.d(xyz, outSupp
                        + SupportRegion3(mFieldIndices.d(xyz)));
                }
                else
                    mFieldIndices.d(xyz, outSupp);
            }
            
            if (output.includes(Field(kB, xyz)))
            {
                outSupp = support(output.regions(), extents().physicalYeeCells().num());
                
                if (mFieldIndices.b(xyz).numRuns() > 0)
                {
                    mFieldIndices.b(xyz, outSupp
                        + SupportRegion3(mFieldIndices.b(xyz)));
                }
                else
                    mFieldIndices.b(xyz, outSupp);
            }
        }
    }
        
    // What does this do?
    // It iterates through the grid denizens and determines which cells of the
    // grid need to store J or M (currents).  It should also check all the
    // outputs and determine the extra locations on the grid where I need to
    // store D and B.
    calculateCurrentIndices(*gridDescription());
    
//    for (int xyz = 0; xyz < 3; xyz++)
//    {
//        LOG << "Direction " << char('X'+xyz) << ":\n";
//        LOG << "E " << mFieldIndices.e(xyz) << "\n";
//        LOG << "H " << mFieldIndices.h(xyz) << "\n";
//        LOG << "D " << mFieldIndices.d(xyz) << "\n";
//        LOG << "B " << mFieldIndices.b(xyz) << "\n";
//    }
}

void VoxelizedPartition::
calculateCurrentIndices(const GridDescription & desc)
{
    Vector3i dims = mExtents.allocatedYeeCells().num();
    
    // J, M
    for (int xyz = 0; xyz < 3; xyz++)
    {
        SupportRegion3 supportJ(dims[0] > 1, dims[1] > 1, dims[2] > 1);
        SupportRegion3 supportM(dims[0] > 1, dims[1] > 1, dims[2] > 1);
        SupportRegion3 supportJE(dims[0] > 1, dims[1] > 1, dims[2] > 1);
        SupportRegion3 supportMH(dims[0] > 1, dims[1] > 1, dims[2] > 1);
        for (int nn = 0; nn < desc.currentSources().size(); nn++)
        {
            SupportRegion3 srJ = support(*desc.currentSources()[nn],
                mExtents.allocatedYeeCells().num(), octantE(xyz));
            supportJ += srJ;
            
            // TODO: consider only actual JE sources, not any ol' J
            SupportRegion3 srJE = support(*desc.currentSources()[nn],
                mExtents.allocatedYeeCells().num(), octantE(xyz));
            supportJE += srJE;
            
            SupportRegion3 srM = support(*desc.currentSources()[nn],
                mExtents.allocatedYeeCells().num(), octantH(xyz));
            supportM += srM;
            
            // TODO: consider only actual MH sources, not any ol' M
            SupportRegion3 srMH = support(*desc.currentSources()[nn],
                mExtents.allocatedYeeCells().num(), octantH(xyz));
            supportMH += srMH;
        }
        for (int nn = 0; nn < desc.huygensSurfaces().size(); nn++)
        {
            // TODO: make sure adjoint Huygens surfaces use JE and MH (right?)
            SupportRegion3 srJ = support(*desc.huygensSurfaces()[nn],
                mExtents.allocatedYeeCells().num(), octantE(xyz));
            supportJ += srJ;
            SupportRegion3 srM = support(*desc.huygensSurfaces()[nn],
                mExtents.allocatedYeeCells().num(), octantH(xyz));
            supportM += srM;
        }
        
        // Include the PML as currents as well, for now
        {
            SupportRegion3 pmlJ = support(
                halfToYee(extents().calcHalfCells(), octantE(xyz)),
                Vector3b(dims[0] > 1, dims[1] > 1, dims[2] > 1));
            pmlJ = pmlJ - support(
                halfToYee(extents().nonPMLHalfCells(), octantE(xyz)),
                Vector3b(dims[0] > 1, dims[1] > 1, dims[2] > 1));
            
            if (!UserPreferences::defines("adjoint"))
                supportJ += pmlJ;
            else
                supportJE += pmlJ;
        }

        {
            SupportRegion3 pmlM = support(
                halfToYee(extents().calcHalfCells(), octantH(xyz)),
                Vector3b(dims[0] > 1, dims[1] > 1, dims[2] > 1));
            pmlM = pmlM - support(
                halfToYee(extents().nonPMLHalfCells(), octantH(xyz)),
                Vector3b(dims[0] > 1, dims[1] > 1, dims[2] > 1));
            
            if (!UserPreferences::defines("adjoint"))
                supportM += pmlM;
            else
                supportMH += pmlM;
        }
            
        
        mFieldIndices.j(xyz, supportJ);
        mFieldIndices.je(xyz, supportJE);
        mFieldIndices.m(xyz, supportM);
        mFieldIndices.mh(xyz, supportMH);
    }
}

std::ostream &
operator<< (std::ostream & out, const VoxelizedPartition & grid)
{
    out << "[VoxelizedPartition]";
	return out;
}






//                          ..:,MMM...N ..MMDM.          ....................  .. 
//                      ,..M..,,     MM,N,D.MM          .O:::~,,Z:.O:,O~:.Z,.Z:,. 
//                 .,:MM..,      MMMMMNNM,,.. ,M.       .Z:,::,.D:.Z:.O=,,Z:.Z,.  
//                 ,M. .          ,.M M . .     M..      8..:,..O..8..O..,O..O..  
//                MM             .N.MMM          M ,     ..........,............  
//             .. N.             .MM,MM   .  ..  ,MM     ZOO,ZZ.OO,$OODOO,OO:ZZZ  
//       .  .  .N,.              MNMM M  ... .M  ...N,   ZZ,:ZZ.OZ.OZ:OOO,OO.ZO.  
//       .MMMMMMNN....          ,MMMMM...M.,.. N.  ,.M.. O..,:..8..Z..8,.,O..O..  
//           MM  .,MMMMN,        .MM. :,M.  . ,  :. N.M  .................. ....  
//       M. N .     ,.MM.       MMMN ,M..   .M. . N ~MM. ..$O,8Z88Z$,ZZ:ZZOOZZ..  
//       M.M..    .   .N.       M M,.N.     .MM    M ,MM ..Z$,O$,OZD,ZO.$O,OZZ..  
//       NM,.    ..NMMMMM.,     ,M .M.       .. ... M  M ..Z..O.,O,.,O..Z..O....  
//      .M    ..MM      :MN.        M               ,...M,    ....................
//      NN     M.        ,MM:.       M,                 MM.   .~ZO,OZ~OOOZ..,O.,,.
//    M.,M    MNMM,       ,M,M       .M .                MN   .:OO.O8.ZO~Z.,:O,::.
//   .M.M    ,....N        .          ,M .          .N   M N. .:..,Z..O,.:$O::8O,.
//   .M M      MMM.        ,N   .   ..,,,M.:.,    .,M     MM:.....................
//   .MM,     M          .,M ..NM   . ...M,,,MMMMMM       M.N.,...::...~,.,:....  
//   .MM      .           M..M,.D,  .MM.M .MM   . .       .MN .:,.OO.,:O,,:O.::.  
//   .MM      ,M.      . M .N.  :N ....,,.    .M..MN.     .MM .,::Z=::,=~8,~:O,.  
//   .MM        M      M   N:D.     M.           ...M,     .M.  ..................
//.  .MN         ,....    :MMM .:M                  :N     .N.  ...,.,.....,...,..
//   .MM                  ..  .:.,                    M    ,N.  .O.Z.Z.Z.O.Z.Z.Z..
//.  .MM                 ,N,M.            ..         .,M.   N.  .O.Z.Z.O.$.O.Z.O..
//   .MM~                .N         ..  . D.,M,, .     ..   M  ...................
//   .M,N                .M        ..MD.,     .M.,MMM  ..   M................,..  
//     ,M                 .        M...        NN.. .M~ .   M..:..OO,.:O,::O.::.. 
//     .N.               M.      .M :    ...,... M    .M.   M..,,,OO,:,O,,:O.::.. 
//     M N              .M       .M.   .MN.M M,..M. ...M.  .....,,,..,.......,..  
//     N:D              ,        .    .M.   N..M.....,       .. ..................
//     .M,N.           .M         . M.,N,    M .            .N. .Z.Z.Z.O.O.Z.Z.Z..
//      ,MM,           M.       .M.,,D,M..NM,              .NM  .Z.Z.Z.Z.Z.Z.Z.Z..
//       .M..         .N .......M    .,.M, .              .MM   ..:,..:...:,.,~...
//         M,:        N.. M. M.M. .   M,.                 .M, .,....... ........  
//        ,..M        N.M,.   ,M ,NM..                    M.. .,:.Z::.Z::.O:,,Z.  
//         M,M  ..    ...,.MMDMNM...                     .M   .,:.Z::.Z::,Z:,,O...
//        .M... M .                                     .M.   ..~:,.=~:.,,..::....
//          M.    .M,,                                  ,M    ..................  
//           M     ::.                                  M.                        
//                                                                 GlassGiant.com
        
