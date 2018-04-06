/*
 *  ForwardGridOperations.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 7/7/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "ForwardGridOperations.h"

#include "operations/ForwardEH.h"
#include "operations/ForwardEH_0.h"
#include "operations/ForwardEH_i.h"
#include "operations/ForwardEH_j.h"
#include "operations/ForwardEH_templated.h"
#include "operations/ForwardNondispersiveEH.h"

#include "operations/ForwardDB.h"

#include "Log.h"
#include "rle/MergeRunlines.h"
#include "Support.h"
#include "YeeUtilities.h"
#include "PhysicalConstants.h"
#include "RLEOperations.h"
#include "UserPreferences.h"

#include "EncodeOutput.h"

#include "ForwardRunlineCallbacks.h"
#include "UniversalRunlineCallbacks.h"

#include "TimeWrapper.h"

using namespace std;
using namespace RLE;
using namespace YeeUtilities;
using namespace TimeWrapper;

ForwardGridOperations::
ForwardGridOperations(GridDescPtr gridDesc,
    const map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids) :
    GridOperations(gridDesc, voxelizedGrids)
{
    VoxelizedPartitionPtr voxelGrid = voxelizedGrids.find(gridDesc)->second;
    
//    if (voxelGrid->gridType() == VoxelizedPartition::kUserGrid)
//    {
//        LOG << "Add boundary outputs\n";
//        addBoundaryOutputs(*gridDesc, *voxelGrid);
//    }
    
    LOG << "Preparing outputs\n";
    // All outputs are potentially filters
    // Encode what to write?  How hard can this be?
    encodeOutputs(*gridDesc, *voxelGrid);
    
    LOG << "Preparing sources\n";
    // These are strictly hard sources... run length encode them too.  (-:
    encodeSources(*gridDesc, *voxelGrid);
    
    LOG << "Preparing current sources\n";
    // Same deal.  They write to J.
    encodeCurrentSources(*gridDesc, *voxelGrid);
    
    LOG << "Preparing periodic boundaries\n";
    encodePeriodicBCs(*voxelGrid);
    
    if (false == UserPreferences::defines("noPMLE"))
    {
        LOG << "PMLJ\n";
        encodePMLJ(*voxelGrid);
    }
    
    LOG << "Preparing Huygens surfaces\n";
    // Use the support of the WRITE region as a mask in the READ grid.
    // Anyway this includes E and H, both.
    encodeHuygens(*gridDesc, *voxelGrid, voxelizedGrids);
    
    LOG << "Preparing D updates\n";
    // Separate D with J and D without J
    // Runlength encode for each of x y z
    encodeD(*voxelGrid);
    
    LOG << "Preparing fast nondispersive isotropic E updates\n";
    encodeNondispersiveE(*voxelGrid);
    
    LOG << "Preparing E updates\n";
    // Diagonal D -> E
    // Off-diagonal: average D, filter, add to diagonal E
    // We save off-diagonal E because it has to add to more than one diagonal E
    // (OR not save it?  huh.)
    //
    // Encode:
    //  - Diagonal D to E, same for all diagonals
    //  - Off-diagonal D to E, which requires some averaging of D
    //  - Summation of anisotropic E to diagonal E (SAME AS D TO E?)
    // This is several things.  I guess it's three different updates...
    // HEY it looks like I do not need to save Exx etc. because they actually
    // can write straight into the diagonal Es.
    encodeEi(*voxelGrid);
    encodeEii(*voxelGrid);
    
    string oct = UserPreferences::valueAs<string>("offDiagonalOctant");
    if (oct == "0")
        encodeEij_0(*voxelGrid);
    else if (oct == "i")
        encodeEij_i(*voxelGrid);
    else if (oct == "j")
        encodeEij_j(*voxelGrid);
    else
        cerr << "Off-diagonal octant 7 is not implemented yet.\n";
    
    if (false == UserPreferences::defines("noPMLH"))
    {
        LOG << "Preparing PML\n";
        encodePMLM(*voxelGrid);
    }
    
    LOG << "Preparing B updates\n";
    // Separate B with and without M
    // Runlength encode for each of x y z
    encodeB(*voxelGrid);
    
    LOG << "Preparing fast nondispersive isotropic H updates\n";
    encodeNondispersiveH(*voxelGrid);
    
    LOG << "Preparing H updates\n";
    // Diagonal B -> H
    // Off-diagonal: average B, filter, add to diagonal H
    //
    // Encode:
    //  - diagonal B to H, same for all diagonals
    //  - Off-diagonal B to H, which requires some averaging of B
    //  - Summation of anisotropic H to diagonal H (SAME AS B TO H?)
    encodeH(*voxelGrid);
    
//    for (int ii = 0; ii < mUpdatesE.size(); ii++)
//    {
//        std::cout << "UpdateE " << ii << ": " << mUpdatesE.at(ii)->name() << ".\n";
//    }
    
    initPerformance();
}

void ForwardGridOperations::
initPerformance()
{
    mTimeClearJM = 0.0;
    mTimeE = vector<double>(mUpdatesE.size(), 0.0);
    mTimeH = vector<double>(mUpdatesH.size(), 0.0);
    mTimeSourceJM = vector<double>(mSourceJM.size(), 0.0);
    mTimeSourceEH = vector<double>(mSourceEH.size(), 0.0);
    mTimeBCE = vector<double>(mPeriodicBCE.size(), 0.0);
    mTimeBCH = vector<double>(mPeriodicBCH.size(), 0.0);
    mTimeOutputs = vector<double>(mOutputs.size(), 0.0);
}

void ForwardGridOperations::
printPerformance(double numT) const
{
    cout << "Clear JM: " << mTimeClearJM/1e6 << " sec\n";
    
    double per_us = 1.0/numT/1e6;
    
    // microseconds:
    double totalE = 0, totalH = 0, totalSourceJM=0, totalSourceEH=0, totalPeriodicBC=0, totalOutput=0;
    
    cout << "Update E:\n";
    for (int nn = 0; nn < mUpdatesE.size(); nn++)
    {
        cout << "\t" << mUpdatesE[nn]->name() << " " << mTimeE[nn]*per_us << " sec/timestep, ";
        cout << mTimeE[nn]*1e3/mUpdatesE[nn]->numCells()/numT << " ns/cell/timestep\n";
        totalE += mTimeE[nn]*1e-6;
    }
    
    cout << "Update H:\n";
    for (int nn = 0; nn < mUpdatesH.size(); nn++)
    {
        cout << "\t" << mUpdatesH[nn]->name() << " " << mTimeH[nn]*per_us << " sec/timestep, ";
        cout << mTimeH[nn]*1e3/mUpdatesH[nn]->numCells()/numT << " ns/cell/timestep\n";
        totalH += mTimeH[nn]*1e-6;
    }
    
    cout << "Source JM:\n";
    for (int nn = 0; nn < mSourceJM.size(); nn++)
    {
        cout << "\t" << mSourceJM[nn]->name() << " " << mTimeSourceJM[nn]*per_us << " sec/timestep, ";
        cout << mTimeSourceJM[nn]*1e3/mSourceJM[nn]->numCells()/numT << " ns/cell/timestep\n";
        totalSourceJM += mTimeSourceJM[nn]*1e-6;
    }
    
    cout << "Source EH:\n";
    for (int nn = 0; nn < mSourceEH.size(); nn++)
    {
        cout << "\t" << mSourceEH[nn]->name() << " " << mTimeSourceEH[nn]*per_us << " sec/timestep, ";
        cout << mTimeSourceEH[nn]*1e3/mSourceEH[nn]->numCells()/numT << " ns/cell/timestep\n";
        totalSourceEH += mTimeSourceEH[nn]*1e-6;
    }
    
    cout << "Periodic BC E:\n";
    for (int nn = 0; nn < mPeriodicBCE.size(); nn++)
    {
        cout << "\t" << mPeriodicBCE[nn]->name() << " " << mTimeBCE[nn]*per_us << " sec/timestep, ";
        cout << mTimeBCE[nn]*1e3/mPeriodicBCE[nn]->numCells()/numT << " ns/cell/timestep\n";
        totalPeriodicBC += mTimeBCE[nn]*1e-6;
    }
    
    cout << "Periodic BC H:\n";
    for (int nn = 0; nn < mPeriodicBCH.size(); nn++)
    {
        cout << "\t" << mPeriodicBCH[nn]->name() << " " << mTimeBCH[nn]*per_us << " sec/timestep, ";
        cout << mTimeBCH[nn]*1e3/mPeriodicBCH[nn]->numCells()/numT << " ns/cell/timestep\n";
        totalPeriodicBC += mTimeBCH[nn]*1e-6;
    }
    
    cout << "Outputs:\n";
    for (int nn = 0; nn < mOutputs.size(); nn++)
    {
        cout << "\t" << mOutputs[nn]->fileName() << " " << mTimeOutputs[nn]*per_us << " sec/timestep\n";
        totalOutput += mTimeOutputs[nn]*1e-6;
    }
    
    cout << "TOTALS:\n";
    cout << "\tE           " << totalE << "\n";
    cout << "\tH           " << totalH << "\n";
    cout << "\tSource EH   " << totalSourceEH << "\n";
    cout << "\tSource JM   " << totalSourceJM << "\n";
    cout << "\tPeriodic BC " << totalPeriodicBC << "\n";
    cout << "\tOutput      " << totalOutput << "\n";
}

void ForwardGridOperations::
handleSetPointers(map<int, Pointer<GridFields> > & fields)
{
    assert(fields.count(gridDescription()->id()));
    GridFields & gridFields = *fields[gridDescription()->id()];
    
    for (int nn = 0; nn < mSourceJM.size(); nn++)
        mSourceJM[nn]->setPointers(gridFields, fields);
    
    for (int nn = 0; nn < mSourceEH.size(); nn++)
        mSourceEH[nn]->setPointers(gridFields, fields);
    
    // March 2011: sending GridFields in to each output, because doing it out
    // here is going to be way too complicated.
    for (int nn = 0; nn < mOutputs.size(); nn++)
        mOutputs[nn]->setPointers(gridFields);
    
    for (int nn = 0; nn < mUpdatesE.size(); nn++)
        mUpdatesE[nn]->setPointers(gridFields, fields);
    for (int nn = 0; nn < mUpdatesH.size(); nn++)
        mUpdatesH[nn]->setPointers(gridFields, fields);
    
    for (int nn = 0; nn < mPeriodicBCE.size(); nn++)
        mPeriodicBCE[nn]->setPointers(gridFields, fields);
    for (int nn = 0; nn < mPeriodicBCH.size(); nn++)
        mPeriodicBCH[nn]->setPointers(gridFields, fields);
    
}

void ForwardGridOperations::
allocate()
{
    for (int nn = 0; nn < mUpdatesE.size(); nn++)
        mUpdatesE[nn]->allocate();
    for (int nn = 0; nn < mUpdatesH.size(); nn++)
        mUpdatesH[nn]->allocate();
    for (int nn = 0; nn < mSourceJM.size(); nn++)
        mSourceJM[nn]->allocate();
}

void ForwardGridOperations::
firstHalfStep(long timestep, Precision::Float dt)
{
    clearJ();
    updateJSource(timestep, dt);
    
    TimePoint t0;
    
    // sourceE not handled yet!
    for (int nn = 0; nn < mUpdatesE.size(); nn++)
    {
        t0 = now();
        mUpdatesE[nn]->apply(timestep, dt);
        mTimeE[nn] += elapsedMicroseconds(t0, now());
    }
    
    for (int nn = 0; nn < mSourceEH.size(); nn++)
    {
        t0 = now();
        mSourceEH[nn]->calcE(timestep);
        mTimeSourceEH[nn] += elapsedMicroseconds(t0, now());
    }
        
    for (int nn = 0; nn < mPeriodicBCE.size(); nn++)
    {
        t0 = now();
        mPeriodicBCE[nn]->apply(timestep, dt);
        mTimeBCE[nn] += elapsedMicroseconds(t0, now());
    }
}

void ForwardGridOperations::
secondHalfStep(long timestep, Precision::Float dt)
{
    clearM();
    updateMSource(timestep, dt);
    
    TimePoint t0;
    
    for (int nn = 0; nn < mUpdatesH.size(); nn++)
    {
        t0 = now();
        mUpdatesH[nn]->apply(timestep, dt);
        mTimeH[nn] += elapsedMicroseconds(t0, now());
    }
        
    for (int nn = 0; nn < mSourceEH.size(); nn++)
    {
        t0 = now();
        mSourceEH[nn]->calcH(timestep);
        mTimeSourceEH[nn] += elapsedMicroseconds(t0, now());
    }
    
    for (int nn = 0; nn < mPeriodicBCH.size(); nn++)
    {
        t0 = now();
        mPeriodicBCH[nn]->apply(timestep, dt);
        mTimeBCH[nn] += elapsedMicroseconds(t0, now());
    }
}

void ForwardGridOperations::
output(long timestep)
{
    TimePoint t0;
    for (int oo = 0; oo < mOutputs.size(); oo++)
    {
        t0 = now();
        mOutputs[oo]->calc(timestep);
        mTimeOutputs[oo] += elapsedMicroseconds(t0, now());
    }
}

void ForwardGridOperations::
updateSourceE(long timestep, Precision::Float dt)
{
    TimePoint t0;
    for (int nn = 0; nn < mSourceEH.size(); nn++)
    {
        t0 = now();
        mSourceEH[nn]->calcE(timestep);
        mTimeSourceEH[nn] += elapsedMicroseconds(t0, now());
    }
}

void ForwardGridOperations::
clearM()
{
    TimePoint t0 = now();
    fields().clearM();
    mTimeClearJM += elapsedMicroseconds(t0, now());
}

void ForwardGridOperations::
updateMSource(long timestep, Precision::Float dt)
{
    TimePoint t0;
    for (int nn = 0; nn < mSourceJM.size(); nn++)
    {
        t0 = now();
        mSourceJM[nn]->calcM(timestep);
        mTimeSourceJM[nn] += elapsedMicroseconds(t0, now());
    }
}


void ForwardGridOperations::
updateSourceH(long timestep, Precision::Float dt)
{
    TimePoint t0;
    for (int nn = 0; nn < mSourceEH.size(); nn++)
    {
        t0 = now();
        mSourceEH[nn]->calcH(timestep);
        mTimeSourceEH[nn] += elapsedMicroseconds(t0, now());
    }
}

void ForwardGridOperations::
clearJ()
{
    TimePoint t0 = now();
    fields().clearJ();
    mTimeClearJM += elapsedMicroseconds(t0, now());
}

void ForwardGridOperations::
updateJSource(long timestep, Precision::Float dt)
{
    TimePoint t0;
    for (int nn = 0; nn < mSourceJM.size(); nn++)
    {
        t0 = now();
        mSourceJM[nn]->calcJ(timestep);
        mTimeSourceJM[nn] += elapsedMicroseconds(t0, now());
    }
}

long ForwardGridOperations::
bytes() const
{
    long totalBytes = 0;
    
    for (int nn = 0; nn < mUpdatesE.size(); nn++)
        totalBytes += mUpdatesE[nn]->bytes();
    for (int nn = 0; nn < mUpdatesH.size(); nn++)
        totalBytes += mUpdatesH[nn]->bytes();
    
    for (int nn = 0; nn < mSourceJM.size(); nn++)
        totalBytes += mSourceJM[nn]->bytes();
    for (int nn = 0; nn < mSourceEH.size(); nn++)
        totalBytes += mSourceEH[nn]->bytes();
    
    for (int nn = 0; nn < mPeriodicBCE.size(); nn++)
        totalBytes += mPeriodicBCE[nn]->bytes();
    for (int nn = 0; nn < mPeriodicBCH.size(); nn++)
        totalBytes += mPeriodicBCH[nn]->bytes();
    
    return totalBytes;
}

void ForwardGridOperations::
printRunlines() const
{
    cout << "\tUpdate D and E:\n";
    for (int nn = 0; nn < mUpdatesE.size(); nn++)
        mUpdatesE[nn]->printRunlines(cout);
    
    cout << "\tUpdate B and H:\n";
    for (int nn = 0; nn < mUpdatesH.size(); nn++)
        mUpdatesH[nn]->printRunlines(cout);
    
    cout << "\tSource J and M:\n";
    for (int nn = 0; nn < mSourceJM.size(); nn++)
        mSourceJM[nn]->printRunlines(cout);
    
    cout << "\tSource E and H:\n";
    for (int nn = 0; nn < mSourceEH.size(); nn++)
        mSourceEH[nn]->printRunlines(cout);
    
    cout << "\tOutputs:\n";
    for (int nn = 0; nn < mOutputs.size(); nn++)
        mOutputs[nn]->printRunlines(cout);
}




void ForwardGridOperations::
addBoundaryOutputs(GridDescription & gridDesc,
    const VoxelizedPartition & voxelGrid)
{
    // E, D
    vector<Duration> durations(1, Duration(0, gridDesc.numTimesteps()-1));
    
    if (false == UserPreferences::defines("sensitivity"))
        return;
    
    for (int ii = 0; ii < 3; ii++)
    {
        if (false == voxelGrid.modeUsesField(Field::e(ii)))
        {
            continue;
        }
        
        int jj = ii;
    //for (int jj = 0; jj < 3; jj++)
    //if (voxelGrid.sensitiveCells(ii,jj).numRuns() > 0)
    {
//        int octant = octantE(ii,jj);
        vector<Region> regions(voxelGrid.sensitiveCells(ii,jj).numRuns());
        SupportRegion3::ConstIterator itr;
        int rr;
        for (rr = 0, itr = voxelGrid.sensitiveCells(ii,jj).begin();
            itr != voxelGrid.sensitiveCells(ii,jj).end();
            rr++, itr.nextMarkedRun())
        {
            Rect3i yeeCells;
            voxelGrid.sensitiveCells(ii,jj).cartesianCoordinates(
                itr.runStart(), yeeCells.p1[0], yeeCells.p1[1], yeeCells.p1[2]);
            voxelGrid.sensitiveCells(ii,jj).cartesianCoordinates(
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
        gridDesc.pushOutput(out);
    }
    }
}


void ForwardGridOperations::
encodeOutputs(const GridDescription & gridDesc, const VoxelizedPartition & vp)
{
    EncodeOutput encoder(gridDesc, vp);
    
    for (int nn = 0; nn < gridDesc.outputs().size(); nn++)
    {
        LOG << "Output " << nn << "\n";
        
        Pointer<OutputEHDB> output = encoder.encode(*gridDesc.outputs()[nn]);
        mOutputs.push_back(output);
    }
}

void ForwardGridOperations::
encodeSources(const GridDescription & gridDesc,
    const VoxelizedPartition & voxelGrid)
{
    for (int nn = 0; nn < gridDesc.sources().size(); nn++)
    {
        LOG << "Source " << nn << "\n";
        const SourceDescription & sourceDesc =
            *(gridDesc.sources())[nn];
        
        Pointer<SourceEH> src(new SourceEH(sourceDesc));
        mSourceEH.push_back(src);
        
        // E sources
        for (int xyz = 0; xyz < 3; xyz++)
        if (sourceDesc.sourceFields().whichE()[xyz])
        if (voxelGrid.modeUsesField(Field::e(xyz)))
        {
            SupportRegion3 supp = support(sourceDesc,
                voxelGrid.extents().nodeYeeCells().num(),
                octantE(xyz));
            IndexArray3 indicesE;
            restriction(voxelGrid.fieldIndices().e(xyz), supp, indicesE);
            
            IndexArray3 const* arrays[] = { &indicesE };
            
            vector<RunlineSource> runlines;
            CallbackSource callback(runlines);
            RLE::merge(arrays, 1, callback);
            src->runlinesE(xyz, runlines);
        }
        
        // H sources
        for (int xyz = 0; xyz < 3; xyz++)
        if (sourceDesc.sourceFields().whichH()[xyz])
        if (voxelGrid.modeUsesField(Field::h(xyz)))
        {
            SupportRegion3 supp = support(sourceDesc,
                voxelGrid.extents().nodeYeeCells().num(),
                octantH(xyz));
            IndexArray3 indicesH;
            restriction(voxelGrid.fieldIndices().h(xyz), supp, indicesH);
            
            IndexArray3 const* arrays[] = { &indicesH };
            
            vector<RunlineSource> runlines;
            CallbackSource callback(runlines);
            RLE::merge(arrays, 1, callback);
            src->runlinesH(xyz, runlines);
        }
    }
}

void ForwardGridOperations::
encodeCurrentSources(const GridDescription & gridDesc,
    const VoxelizedPartition & voxelGrid)
{
    for (int nn = 0; nn < gridDesc.currentSources().size(); nn++)
    {
        LOG << "Current source " << nn << "\n";
        const CurrentSourceDescription & sourceDesc =
            *(gridDesc.currentSources())[nn];
        
        Pointer<SourceJM> src(new SourceJM(sourceDesc));
        mSourceJM.push_back(src);
        
        // J sources
        for (int xyz = 0; xyz < 3; xyz++)
        if (sourceDesc.sourceCurrents().whichJ()[xyz])
        if (voxelGrid.modeUsesField(Field::j(xyz)))
        {
            SupportRegion3 supp = support(sourceDesc,
                voxelGrid.extents().nodeYeeCells().num(),
                octantE(xyz));
            IndexArray3 indicesJ;
            restriction(voxelGrid.fieldIndices().j(xyz), supp, indicesJ);
            
            IndexArray3 const* arrays[] = { &indicesJ };
            
            vector<RunlineSource> runlines;
            CallbackSource callback(runlines);
            RLE::merge(arrays, 1, callback);
            src->runlinesJ(xyz, runlines);
        }
        
        // M sources
        for (int xyz = 0; xyz < 3; xyz++)
        if (sourceDesc.sourceCurrents().whichM()[xyz])
        if (voxelGrid.modeUsesField(Field::m(xyz)))
        {
            SupportRegion3 supp = support(sourceDesc,
                voxelGrid.extents().nodeYeeCells().num(),
                octantH(xyz));
            IndexArray3 indicesM;
            restriction(voxelGrid.fieldIndices().m(xyz), supp, indicesM);
            
            IndexArray3 const* arrays[] = { &indicesM };
            
            vector<RunlineSource> runlines;
            CallbackSource callback(runlines);
            RLE::merge(arrays, 1, callback);
            src->runlinesM(xyz, runlines);
        }
    }
}

static void sCalculateHuygensRunlines(
    vector<RunlineHuygens> & outRunlines,
    const HuygensSurfaceDescription & huyg,
    const IndexArray3 & destCurrent,
    const IndexArray3 & sourceField,
    int face, int xyz, int oct)
{
    Rect3i sideCells = huyg.yeeCells(face, oct);
    
    SupportRegion3 maskDest(support(sideCells)), maskSource;
    if (huyg.yeeCells().encloses(sideCells))
    {
//        LOG << "Side " << face << " e" << char('x'+xyz)
//            << " is total-field.\n";
        translate(maskDest, maskSource, cardinal(face).asArray());
    }
    else
    {
//        LOG << "Side " << face << " e" << char('x'+xyz)
//            << " is scattered-field.\n";
        translate(maskDest, maskSource, (-cardinal(face)).asArray());
    }
    //LOG << "Sign: " << signPermutation*signFace*signEH << "\n";
    
//    cout << "Source field:\n" << sourceField << "\n";
//    cout << "Source mask:\n" << maskSource << "\n";
    
    IndexArray3 indicesJM, indicesHE;
    restriction(destCurrent, maskDest, indicesJM);
    restriction(sourceField, maskSource, indicesHE);
    
//    cout << "Source indices:\n" << indicesHE << "\n";
//    cout << "Dest indices:\n" << indicesJM << "\n";
    
    // That's less trivial than it looks.  Use the mask for H in the
    // main grid to extract some of the H field indices in the aux grid
    // and save the aux grid indices.
    IndexArray3 const* arrays[] = { &indicesJM, &indicesHE };
    
    CallbackHuygensSurface callback(outRunlines);
    RLE::merge(arrays, 2, callback);
//    RLE::merge(indices, strides, lengths, 2, callback);
}
        
void ForwardGridOperations::
encodeHuygens(const GridDescription & gridDesc,
    const VoxelizedPartition & voxelGrid,
    const map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids)
{
    for (int nn = 0; nn < gridDesc.huygensSurfaces().size(); nn++)
    {
        const HuygensSurfaceDescription & huyg = 
            *(gridDesc.huygensSurfaces()[nn]);
        assert(voxelizedGrids.count(huyg.sourceGrid()));
        
        const VoxelizedPartition & sourceVoxels =
            *voxelizedGrids.find(huyg.sourceGrid())->second;
        
        LOG << "TF cells " << huyg.yeeCells() << "\n";
        
        for (int face = 0; face < 6; face++)
        {
            if (huyg.omitsSide(face))
            {
                continue;
            }
            
            for (int fieldXYZ = 0; fieldXYZ < 3; fieldXYZ++)
            {
                if (fieldXYZ == face/2)
                {
                    // e.g. Ex is not modified on the +x and -x faces
                    continue;
                }
                
                int srcXYZ = 3 - (fieldXYZ + face/2); // perp. to E and to the face normal.
                assert(srcXYZ != fieldXYZ);
                assert(srcXYZ != face/2);
                
                if (voxelGrid.modeUsesField(Field::e(fieldXYZ)) &&
                    voxelGrid.modeUsesField(Field::h(srcXYZ)))
                {
                    // Sign factors:
                    double signEH = 1; // +1 E, -1 H
                    double signPermute = srcXYZ == (fieldXYZ+1)%3 ? 1 : -1;
                    double signFace = face%2 ? 1 : -1;
                    
                    const IndexArray3 & destCurrent = voxelGrid.fieldIndices().j(fieldXYZ);
                    const IndexArray3 & sourceField = sourceVoxels.fieldIndices().h(srcXYZ);
                    
                    vector<RunlineHuygens> runlines;
                    sCalculateHuygensRunlines(runlines, huyg, destCurrent, sourceField,
                        face, fieldXYZ, octantE(fieldXYZ));
                    
    //                if (runlines.size() > 0)
    //                {
    //                    LOG << "Huygens E" << char('x'+fieldXYZ) << " has "
    //                        << runlines.size() << " output runs.\n";
    //                    for (int nn = 0; nn < runlines.size(); nn++)
    //                        cout << runlines.at(nn) << " ";
    //                }
                    
                    pushUpdateE(new HuygensSurface(
                        huyg.sourceGrid(),
                        signEH*signPermute*signFace/gridDescription()->dxyz()[face/2],
                        Field::j(fieldXYZ), Field::h(srcXYZ),
                        runlines));
                }
                
                if (voxelGrid.modeUsesField(Field::h(fieldXYZ)) &&
                    voxelGrid.modeUsesField(Field::e(srcXYZ)))
                {
                    // Sign factors:
                    double signEH = -1; // +1 E, -1 H
                    double signPermute = srcXYZ == (fieldXYZ+1)%3 ? 1 : -1;
                    double signFace = face%2 ? 1 : -1;
                    
                    const IndexArray3 & destCurrent = voxelGrid.fieldIndices().m(fieldXYZ);
                    const IndexArray3 & sourceField = sourceVoxels.fieldIndices().e(srcXYZ);
                    
                    vector<RunlineHuygens> runlines;
                    sCalculateHuygensRunlines(runlines, huyg, destCurrent, sourceField,
                        face, fieldXYZ, octantH(fieldXYZ));
                    
    //                if (runlines.size() > 0)
    //                {
    //                    LOG << "Huygens H" << char('x'+fieldXYZ) << " has "
    //                        << runlines.size() << " output runs.\n";
    //                    for (int nn = 0; nn < runlines.size(); nn++)
    //                        cout << runlines.at(nn) << " ";
    //                }
                    
                    pushUpdateH(new HuygensSurface(
                        huyg.sourceGrid(),
                        signEH*signPermute*signFace/gridDescription()->dxyz()[face/2],
                        Field::m(fieldXYZ), Field::e(srcXYZ),
                        runlines));
                }
            }
        }
    }
}

void ForwardGridOperations::
encodePeriodicBCs(const VoxelizedPartition & voxelGrid)
{
    // Find facing pairs of periodic boundaries.
    // Hook 'em up!
    
    const NodeExtents & ext = voxelGrid.extents();
    
    for (int readSide = 0; readSide < 6; readSide++)
    if (ext.nodeBoundaries()[readSide] == kPeriodicBoundary)
    {
        int writeSide;
        if (readSide%2 == 0)
            writeSide = readSide + 1;
        else
            writeSide = readSide - 1;
        
        Rect3i readHalfCells = ext.ghostReadHalfCells(readSide);
        Rect3i writeHalfCells = ext.ghostWriteHalfCells(writeSide);
        
        if (ext.numGhostHalfCells(readSide) > 1)
            throw(std::logic_error("Fat ghost regions are not implemented."));
        
        for (int fieldXYZ = 0; fieldXYZ < 3; fieldXYZ++)
        {
            if (false == voxelGrid.modeUsesField(Field::e(fieldXYZ)))
            {
                continue;
            }
            
            // Encode E:
            Rect3i readYee = halfToYee(readHalfCells, octantE(fieldXYZ));
            Rect3i writeYee = halfToYee(writeHalfCells, octantE(fieldXYZ));
            
            if (readYee.count() <= 0 || writeYee.count() <= 0)
                continue;
            
            //cout << char('x'+fieldXYZ) << ": Copy " << readYee << " to " << writeYee << "\n";
            
            SupportRegion3 suppRead = support(readYee,
                voxelGrid.extents().dimensions());
            SupportRegion3 suppWrite = support(writeYee,
                voxelGrid.extents().dimensions());
            
            IndexArray3 indRead, indWrite;
            
            restriction(voxelGrid.fieldIndices().e(fieldXYZ), suppRead, indRead);
            restriction(voxelGrid.fieldIndices().e(fieldXYZ), suppWrite, indWrite);
            
//            cout << "Copy:\n" << indRead;
//            cout << "\nPaste:\n" << indWrite << "\n";
            
            IndexArray3 const* arrays[] = { &indRead, &indWrite };
            vector<RunlinePeriodicBC> runlines;
            CallbackPeriodicBC callback(runlines);
            RLE::merge(arrays, 2, callback);
            
            mPeriodicBCE.push_back(Pointer<PeriodicBC>(new PeriodicBC(
                Field::e(fieldXYZ),
                runlines)));
        }
        
        for (int fieldXYZ = 0; fieldXYZ < 3; fieldXYZ++)
        {
            if (false == voxelGrid.modeUsesField(Field::h(fieldXYZ)))
            {
                continue;
            }
            // Encode H:
            Rect3i readYee = halfToYee(readHalfCells, octantH(fieldXYZ));
            Rect3i writeYee = halfToYee(writeHalfCells, octantH(fieldXYZ));
            
            if (readYee.count() <= 0 || writeYee.count() <= 0)
                continue;
            
            //cout << char('x'+fieldXYZ) << ": Copy " << readYee << " to " << writeYee << "\n";
            
            SupportRegion3 suppRead = support(readYee,
                voxelGrid.extents().dimensions());
            SupportRegion3 suppWrite = support(writeYee,
                voxelGrid.extents().dimensions());
            
            IndexArray3 indRead, indWrite;
            
            restriction(voxelGrid.fieldIndices().h(fieldXYZ), suppRead, indRead);
            restriction(voxelGrid.fieldIndices().h(fieldXYZ), suppWrite, indWrite);
            
            IndexArray3 const* arrays[] = { &indRead, &indWrite };
            vector<RunlinePeriodicBC> runlines;
            CallbackPeriodicBC callback(runlines);
            RLE::merge(arrays, 2, callback);
            
            mPeriodicBCH.push_back(Pointer<PeriodicBC>(new PeriodicBC(
                Field::h(fieldXYZ),
                runlines)));
        }
    }
}

/*
    Here I need some cooperation between the D and E updates if I want to 
    optimize these routines.  My idea is to NOT store D and NOT split the
    updates, provided that I don't need to do the fancy off-diagonal stuff
    and provided that I don't need to save the D field.
    
    In point of fact, I might just decide that if I don't have any anisotropic
    materials that I consequently ALWAYS try to (say) combine the E and D field
    updates.  That might be an ok optimization in all cases.
    
    There are two things I might want to do.
    1. Economize on memory by not always storing D.
    2. Speed things up by combining the D and E updates
    
    #2 is a little more straightforward.  I can handle this with a suitable
    runline structure.  The basic idea is this: the D = curl H update is run
    alone only for off-diagonal permittivities; elsewhere, D and E can update
    simultaneously.  This will automatically run as fast as possible for my
    simulations now (Feb 2012) because I don't do any off-diagonal epsilons.
*/


// TODO: Use curiously recurring template pattern to permit return-by-value
// for restriction(), translation(), etc.  (Will this work?)
void ForwardGridOperations::
encodeD(const VoxelizedPartition & voxelGrid)
{
    Matrix3i m = voxelGrid.projectionMatrix();
    
    // Without current sources
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (false == voxelGrid.modeUsesField(Field::d(xyz)))
        {
            continue;
        }
        Vector3i shiftY = -m*Vector3i::unit( (xyz+1)%3 );
        Vector3i shiftZ = -m*Vector3i::unit( (xyz+2)%3 );
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 maskD(voxelGrid.fieldIndices().d(xyz));
        maskD -= voxelGrid.fieldIndices().j(xyz);
        maskD *= calcCells;
        SupportRegion3 maskHjPlus, maskHjMinus, maskHkPlus, maskHkMinus;
        maskHjPlus = maskD;
        translate(maskD, maskHjMinus, shiftZ[0], shiftZ[1], shiftZ[2]);
        maskHkPlus = maskD;
        translate(maskD, maskHkMinus, shiftY[0], shiftY[1], shiftY[2]);
        
        IndexArray3 indicesD, indicesHjPlus, indicesHjMinus,
            indicesHkPlus, indicesHkMinus;
        restriction(voxelGrid.fieldIndices().d(xyz), maskD, indicesD);
        restriction(voxelGrid.fieldIndices().h((xyz+1)%3), maskHjMinus,
            indicesHjMinus);
        restriction(voxelGrid.fieldIndices().h((xyz+1)%3), maskHjPlus,
            indicesHjPlus);
        restriction(voxelGrid.fieldIndices().h((xyz+2)%3), maskHkMinus,
            indicesHkMinus);
        restriction(voxelGrid.fieldIndices().h((xyz+2)%3), maskHkPlus,
            indicesHkPlus);
        
//        cerr << char('x'+xyz) << ":\n"
//            << "calc: " << calcCells << "\n"
//            << "mask: " << maskD << "\n";
        
        IndexArray3 const* arrays[] = { &indicesD,
            &indicesHjMinus, &indicesHjPlus,
            &indicesHkMinus, &indicesHkPlus, };
        
        vector<RunlineDB> runlines;
        CallbackDB callback(runlines);
        RLE::merge(arrays, 5, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new ForwardDB(
                Field::d(xyz), Field::h((xyz+1)%3), Field::h((xyz+2)%3),
                gridDescription()->dt()/gridDescription()->dxyz()[(xyz+1)%3],
                -gridDescription()->dt()/gridDescription()->dxyz()[(xyz+2)%3],
                runlines));
        }
    }
    
    // With current sources
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (false == voxelGrid.modeUsesField(Field::d(xyz)))
        {
            continue;
        }
        Vector3i shiftY = -m*Vector3i::unit( (xyz+1)%3 );
        Vector3i shiftZ = -m*Vector3i::unit( (xyz+2)%3 );
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 maskD(voxelGrid.fieldIndices().d(xyz));
        maskD *= voxelGrid.fieldIndices().j(xyz);
        maskD *= calcCells;
        SupportRegion3 maskHjPlus, maskHjMinus, maskHkPlus, maskHkMinus;
        maskHjPlus = maskD;
        translate(maskD, maskHjMinus, shiftZ[0], shiftZ[1], shiftZ[2]);
        maskHkPlus = maskD;
        translate(maskD, maskHkMinus, shiftY[0], shiftY[1], shiftY[2]);
        
        IndexArray3 indicesD, indicesHjPlus, indicesHjMinus,
            indicesHkPlus, indicesHkMinus, indicesJ;
        restriction(voxelGrid.fieldIndices().d(xyz), maskD, indicesD);
        restriction(voxelGrid.fieldIndices().h((xyz+1)%3), maskHjMinus,
            indicesHjMinus);
        restriction(voxelGrid.fieldIndices().h((xyz+1)%3), maskHjPlus,
            indicesHjPlus);
        restriction(voxelGrid.fieldIndices().h((xyz+2)%3), maskHkMinus,
            indicesHkMinus);
        restriction(voxelGrid.fieldIndices().h((xyz+2)%3), maskHkPlus,
            indicesHkPlus);
        restriction(voxelGrid.fieldIndices().j(xyz), maskD, indicesJ);
        
//        cerr << char('x'+xyz) << ":\n"
//            << "calc: " << calcCells << "\n"
//            << "mask: " << maskD << "\n";
        
        IndexArray3 const* arrays[] = { &indicesD,
            &indicesHjMinus, &indicesHjPlus,
            &indicesHkMinus, &indicesHkPlus,
            &indicesJ };
        
        vector<RunlineDB_JM> runlines;
        CallbackDB_JM callback(runlines);
        RLE::merge(arrays, 6, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new ForwardDB_JM(
                Field::d(xyz), Field::h((xyz+1)%3), Field::h((xyz+2)%3),
                Field::j(xyz),
                gridDescription()->dt()/gridDescription()->dxyz()[(xyz+1)%3],
                -gridDescription()->dt()/gridDescription()->dxyz()[(xyz+2)%3],
                -gridDescription()->dt(),
                runlines));
        }
    }
}

void ForwardGridOperations::
encodeNondispersiveE(const VoxelizedPartition & voxelGrid)
{
    Matrix3i m = voxelGrid.projectionMatrix();
    
    // Without current sources
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (false == voxelGrid.modeUsesField(Field::e(xyz)))
        {
            continue;
        }
        
        Vector3i shiftY = -m*Vector3i::unit( (xyz+1)%3 );
        Vector3i shiftZ = -m*Vector3i::unit( (xyz+2)%3 );
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
        
        SupportRegion3 maskNondispersive(voxelGrid.inversePermittivity().support(xyz, xyz, 0, 0));
        
        SupportRegion3 maskE(voxelGrid.inversePermittivity().support(xyz, xyz, 0, 0) - voxelGrid.fieldIndices().d(xyz));
        
        maskE -= voxelGrid.fieldIndices().j(xyz);
        maskE *= maskNondispersive;
        maskE *= calcCells;
        
        SupportRegion3 maskHjPlus, maskHjMinus, maskHkPlus, maskHkMinus;
        maskHjPlus = maskE;
        translate(maskE, maskHjMinus, shiftZ[0], shiftZ[1], shiftZ[2]);
        maskHkPlus = maskE;
        translate(maskE, maskHkMinus, shiftY[0], shiftY[1], shiftY[2]);
        
        IndexArray3 permittivityIndices(voxelGrid.inversePermittivity().filter(
            xyz, xyz, 0, 0));
        
        IndexArray3 indicesE, indicesHjPlus, indicesHjMinus,
            indicesHkPlus, indicesHkMinus, indicesEps;
        restriction(voxelGrid.fieldIndices().e(xyz), maskE, indicesE);
        restriction(voxelGrid.fieldIndices().h((xyz+1)%3), maskHjMinus,
            indicesHjMinus);
        restriction(voxelGrid.fieldIndices().h((xyz+1)%3), maskHjPlus,
            indicesHjPlus);
        restriction(voxelGrid.fieldIndices().h((xyz+2)%3), maskHkMinus,
            indicesHkMinus);
        restriction(voxelGrid.fieldIndices().h((xyz+2)%3), maskHkPlus,
            indicesHkPlus);
        restriction(permittivityIndices, maskE, indicesEps);
        
//        cerr << char('x'+xyz) << ":\n"
//            << "calc: " << calcCells << "\n"
//            << "mask: " << maskE << "\n";
        
        IndexArray3 const* arrays[] = { &indicesE,
            &indicesHjMinus, &indicesHjPlus,
            &indicesHkMinus, &indicesHkPlus,
            &indicesEps};
        
        vector<RunlineNondispersiveEH> runlines;
        CallbackNondispersiveEH callback(runlines);
        RLE::merge(arrays, 6, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new ForwardNondispersiveEH(
                Field::e(xyz), Field::h((xyz+1)%3), Field::h((xyz+2)%3),
                gridDescription()->dt()/gridDescription()->dxyz()[(xyz+1)%3],
                -gridDescription()->dt()/gridDescription()->dxyz()[(xyz+2)%3],
                runlines,
                voxelGrid.inversePermittivity().filter(xyz,xyz,0,0).values(),
                Constants::eps0));
        }
    }
    
    // With current sources
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (false == voxelGrid.modeUsesField(Field::e(xyz)))
        {
            continue;
        }
        
        Vector3i shiftY = -m*Vector3i::unit( (xyz+1)%3 );
        Vector3i shiftZ = -m*Vector3i::unit( (xyz+2)%3 );
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
        
        SupportRegion3 maskNondispersive(voxelGrid.inversePermittivity()
            .support(xyz, xyz, 0, 0));
        
        SupportRegion3 maskE(voxelGrid.inversePermittivity()
            .support(xyz, xyz, 0, 0) - voxelGrid.fieldIndices().d(xyz));
        
        maskE *= voxelGrid.fieldIndices().j(xyz);
        maskE *= maskNondispersive;
        maskE *= calcCells;
        
        SupportRegion3 maskHjPlus, maskHjMinus, maskHkPlus, maskHkMinus;
        maskHjPlus = maskE;
        translate(maskE, maskHjMinus, shiftZ[0], shiftZ[1], shiftZ[2]);
        maskHkPlus = maskE;
        translate(maskE, maskHkMinus, shiftY[0], shiftY[1], shiftY[2]);
        
        IndexArray3 permittivityIndices(voxelGrid.inversePermittivity().filter(
            xyz, xyz, 0, 0));
        
        IndexArray3 indicesE, indicesHjPlus, indicesHjMinus,
            indicesHkPlus, indicesHkMinus, indicesEps, indicesJ;
        restriction(voxelGrid.fieldIndices().e(xyz), maskE, indicesE);
        restriction(voxelGrid.fieldIndices().h((xyz+1)%3), maskHjMinus,
            indicesHjMinus);
        restriction(voxelGrid.fieldIndices().h((xyz+1)%3), maskHjPlus,
            indicesHjPlus);
        restriction(voxelGrid.fieldIndices().h((xyz+2)%3), maskHkMinus,
            indicesHkMinus);
        restriction(voxelGrid.fieldIndices().h((xyz+2)%3), maskHkPlus,
            indicesHkPlus);
        restriction(permittivityIndices, maskE, indicesEps);
        restriction(voxelGrid.fieldIndices().j(xyz), maskE, indicesJ);
        
//        cerr << char('x'+xyz) << ":\n"
//            << "calc: " << calcCells << "\n"
//            << "mask: " << maskE << "\n";
        
        IndexArray3 const* arrays[] = { &indicesE,
            &indicesHjMinus, &indicesHjPlus,
            &indicesHkMinus, &indicesHkPlus,
            &indicesEps, &indicesJ};
        
        vector<RunlineNondispersiveEH_JM> runlines;
        CallbackNondispersiveEH_JM callback(runlines);
        RLE::merge(arrays, 7, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new ForwardNondispersiveEH_JM(
                Field::e(xyz), Field::h((xyz+1)%3), Field::h((xyz+2)%3),
                Field::j(xyz),
                gridDescription()->dt()/gridDescription()->dxyz()[(xyz+1)%3],
                -gridDescription()->dt()/gridDescription()->dxyz()[(xyz+2)%3],
                -gridDescription()->dt(),
                runlines,
                voxelGrid.inversePermittivity().filter(xyz,xyz,0,0).values(),
                Constants::eps0));
        }
    }
}

// Ei is the plain-vanilla E update

void ForwardGridOperations::
encodeEi(const VoxelizedPartition & voxelGrid)
{
    LOG << "Encoding diagonal E\n";
    for (int numerOrder = 0; numerOrder < 10; numerOrder++)
    for (int denomOrder = 0; denomOrder < 10; denomOrder++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (false == voxelGrid.modeUsesField(Field::e(xyz)))
        {
            continue;
        }
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
        
        SupportRegion3 suppEii(voxelGrid.fieldIndices().ee(xyz, xyz));
        SupportRegion3 suppDi(voxelGrid.fieldIndices().d(xyz));
        
        SupportRegion3 maskD(voxelGrid.inversePermittivity().support(xyz, xyz, numerOrder, denomOrder));
        maskD *= calcCells;
        maskD *= suppDi;
        maskD -= suppEii;
        
        IndexArray3 permittivityIndices(
            voxelGrid.inversePermittivity().filter(
            xyz, xyz, numerOrder, denomOrder));
        IndexArray3 indicesD, indicesE, indicesEps, indicesAux(maskD);
        restriction(voxelGrid.fieldIndices().d(xyz), maskD, indicesD);
        restriction(voxelGrid.fieldIndices().e(xyz), maskD, indicesE);
        restriction(permittivityIndices, maskD, indicesEps);
        
        IndexArray3 const* arrays[] =
            { &indicesE, &indicesD, &indicesAux, &indicesEps };
        vector<RunlineEH> runlines;
        CallbackEH callback(runlines);
        RLE::merge(arrays, 4, callback);
        
        if (runlines.size() > 0)
        {
            assert(numerOrder == denomOrder);
            int order = numerOrder;
            
            pushUpdateE(FwdEHFactory::newForwardEH(
                Field::e(xyz), Field::d(xyz),
                order,
                runlines,
                voxelGrid.inversePermittivity().filter(xyz, xyz, numerOrder, denomOrder).values(),
                Constants::eps0));
        }
    }
}

void ForwardGridOperations::
encodeEii(const VoxelizedPartition & voxelGrid)
{
    LOG << "Encoding diagonal EE\n";
    for (int numerOrder = 0; numerOrder < 10; numerOrder++)
    for (int denomOrder = 0; denomOrder < 10; denomOrder++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
        
        SupportRegion3 suppEii(voxelGrid.fieldIndices().ee(xyz, xyz));
        SupportRegion3 suppDi(voxelGrid.fieldIndices().d(xyz));
        
        SupportRegion3 maskD(voxelGrid.inversePermittivity()
            .support(xyz, xyz, numerOrder, denomOrder));
        maskD *= calcCells;
        maskD *= suppDi;
        maskD *= suppEii;
        
        IndexArray3 permittivityIndices(
            voxelGrid.inversePermittivity().filter(
            xyz, xyz, numerOrder, denomOrder));
        IndexArray3 indicesD, indicesE, indicesEps, indicesAux(maskD);
        restriction(voxelGrid.fieldIndices().d(xyz), maskD, indicesD);
        restriction(voxelGrid.fieldIndices().ee(xyz, xyz), maskD, indicesE);
        restriction(permittivityIndices, maskD, indicesEps);
        
        IndexArray3 const* arrays[] =
            { &indicesE, &indicesD, &indicesAux, &indicesEps };
        vector<RunlineEH> runlines;
        CallbackEH callback(runlines);
        RLE::merge(arrays, 4, callback);
        
        if (runlines.size() > 0)
        {
            assert(numerOrder == denomOrder);
            int order = numerOrder;
            pushUpdateE(FwdEHFactory::newForwardEH(
                Field::ee(xyz,xyz), Field::d(xyz),
                order,
                runlines,
                voxelGrid.inversePermittivity()
                    .filter(xyz, xyz, numerOrder, denomOrder).values(),
                Constants::eps0));
        }
    }
}

void ForwardGridOperations::
encodeEij_0(const VoxelizedPartition & voxelGrid)
{
    Matrix3i m = voxelGrid.projectionMatrix();
    
    LOG << "Encoding off-diagonal E (method 0)\n";
    for (int numerOrder = 0; numerOrder < 10; numerOrder++)
    for (int denomOrder = 0; denomOrder < 10; denomOrder++)
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    if (ii != jj)
    {
        const int OCTANT_0 = 0;
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), OCTANT_0),
            voxelGrid.extents().dimensions());
        
        SupportRegion3 suppEij(voxelGrid.fieldIndices().ee(ii, jj));
        
        SupportRegion3 permittivityMask(voxelGrid.inversePermittivity()
            .support(ii,jj, numerOrder, denomOrder));
        
        SupportRegion3 maskD0, maskD1;
        maskD0 = permittivityMask * calcCells;
        
        Vector3i shift = -m*Vector3i::unit(jj);
        translate(maskD0, maskD1, shift[0], shift[1], shift[2]);
        
        IndexArray3 permittivityIndices(
            voxelGrid.inversePermittivity().filter(
            ii, jj, numerOrder, denomOrder));
        IndexArray3 indicesD0, indicesD1, indicesE, indicesEps,
            indicesAux(maskD0);
        restriction(voxelGrid.fieldIndices().d(jj), maskD0, indicesD0);
        restriction(voxelGrid.fieldIndices().d(jj), maskD1, indicesD1);
        restriction(voxelGrid.fieldIndices().ee(ii,jj), maskD0, indicesE);
        restriction(permittivityIndices, maskD0, indicesEps);
        
        IndexArray3 const* arrays[] =
            { &indicesE, &indicesD0, &indicesD1, &indicesAux, &indicesEps };
        vector<RunlineAnisotropicEH_0> runlines;
        CallbackAnisotropicEH_0 callback(runlines);
        RLE::merge(arrays, 5, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new UpdateAnisotropicEH_0(
                Field::ee(ii,jj), Field::d(jj),
                denomOrder,
                numerOrder,
                runlines,
                voxelGrid.inversePermittivity()
                    .filter(ii, jj, numerOrder, denomOrder).values(),
                Constants::eps0));
        }
    }
    
    LOG << "Encoding summation of E (method 0)\n";
    for (int xyz = 0; xyz < 3; xyz++)
    {
        int j = (xyz+1)%3, k = (xyz+2)%3;
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 suppEii(voxelGrid.fieldIndices().ee(xyz,xyz));
        
        SupportRegion3 sumRegion = suppEii*calcCells;
        
        SupportRegion3 sumRegionJ, sumRegionK;
        Vector3i dj(Vector3i::unit(j)), dk(Vector3i::unit(k));
        translate(sumRegion, sumRegionJ, dj[0], dj[1], dj[2]);
        translate(sumRegion, sumRegionK, dk[0], dk[1], dk[2]);
        
        IndexArray3 indicesEi, indicesEii, indicesEij0, indicesEij1, 
            indicesEik0, indicesEik1;
        restriction(voxelGrid.fieldIndices().e(xyz), sumRegion, indicesEi);
        restriction(voxelGrid.fieldIndices().ee(xyz,xyz), sumRegion, indicesEii);
        restriction(voxelGrid.fieldIndices().ee(xyz,j), sumRegion, indicesEij0);
        restriction(voxelGrid.fieldIndices().ee(xyz,j), sumRegionJ, indicesEij1);
        restriction(voxelGrid.fieldIndices().ee(xyz,k), sumRegion, indicesEik0);
        restriction(voxelGrid.fieldIndices().ee(xyz,k), sumRegionK, indicesEik1);
        
        IndexArray3 const* arrays[] = { &indicesEi, &indicesEii,
            &indicesEij0, &indicesEij1, &indicesEik0, &indicesEik1 };
        vector<RunlineSumEH_0> runlines;
        CallbackSumEH_0 callback(runlines);
        RLE::merge(arrays, 6, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new SumEH_0(Field::e(xyz), 
                Field::ee(xyz,xyz), Field::ee(xyz,j), Field::ee(xyz,k),
                runlines));
        }
    }
}

void ForwardGridOperations::
encodeEij_i(const VoxelizedPartition & voxelGrid)
{    
    Matrix3i m = voxelGrid.projectionMatrix();
    
    LOG << "Encoding off-diagonal E (method i)\n";
    for (int numerOrder = 0; numerOrder < 10; numerOrder++)
    for (int denomOrder = 0; denomOrder < 10; denomOrder++)
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    if (ii != jj)
    {
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(ii)),
            voxelGrid.extents().dimensions());
        
        SupportRegion3 suppEij(voxelGrid.fieldIndices().ee(ii, jj));
        
        SupportRegion3 permittivityMask(voxelGrid.inversePermittivity()
            .support(ii,jj, numerOrder, denomOrder));
        
        SupportRegion3 maskD0, maskD1, maskD2, maskD3;
        maskD0 = permittivityMask * calcCells;
        
        Vector3i shiftJ = -m*Vector3i::unit(jj);
        Vector3i shiftI = m*Vector3i::unit(ii);
        Vector3i shiftIJ = shiftJ + shiftI;
        translate(maskD0, maskD1, shiftI[0], shiftI[1], shiftI[2]);
        translate(maskD0, maskD2, shiftJ[0], shiftJ[1], shiftJ[2]);
        translate(maskD0, maskD3, shiftIJ[0], shiftIJ[1], shiftIJ[2]);
        
        IndexArray3 permittivityIndices(
            voxelGrid.inversePermittivity().filter(
            ii, jj, numerOrder, denomOrder));
        IndexArray3 indicesD0, indicesD1, indicesD2, indicesD3,
            indicesE, indicesEps, indicesAux(maskD0);
        restriction(voxelGrid.fieldIndices().d(jj), maskD0, indicesD0);
        restriction(voxelGrid.fieldIndices().d(jj), maskD1, indicesD1);
        restriction(voxelGrid.fieldIndices().d(jj), maskD2, indicesD2);
        restriction(voxelGrid.fieldIndices().d(jj), maskD3, indicesD3);
        restriction(voxelGrid.fieldIndices().ee(ii,jj), maskD0, indicesE);
        restriction(permittivityIndices, maskD0, indicesEps);
        
        IndexArray3 const* arrays[] =
            { &indicesE, &indicesD0, &indicesD1, &indicesD2, &indicesD3,
            &indicesAux, &indicesEps };
        vector<RunlineAnisotropicEH_i> runlines;
        CallbackAnisotropicEH_i callback(runlines);
        RLE::merge(arrays, 7, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new UpdateAnisotropicEH_i(
                Field::ee(ii,jj), Field::d(jj),
                denomOrder,
                numerOrder,
                runlines,
                voxelGrid.inversePermittivity()
                    .filter(ii, jj, numerOrder, denomOrder).values(),
                Constants::eps0));
        }
    }
    
    LOG << "Encoding summation of E (method i)\n";
    for (int xyz = 0; xyz < 3; xyz++)
    {
        int j = (xyz+1)%3, k = (xyz+2)%3;
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 suppEii(voxelGrid.fieldIndices().ee(xyz,xyz));
        
        SupportRegion3 sumRegion = suppEii*calcCells;
        
        IndexArray3 indicesEi, indicesEii, indicesEij, indicesEik;
        restriction(voxelGrid.fieldIndices().e(xyz), sumRegion, indicesEi);
        restriction(voxelGrid.fieldIndices().ee(xyz,xyz), sumRegion, indicesEii);
        restriction(voxelGrid.fieldIndices().ee(xyz,j), sumRegion, indicesEij);
        restriction(voxelGrid.fieldIndices().ee(xyz,k), sumRegion, indicesEik);
        
        IndexArray3 const* arrays[] = { &indicesEi, &indicesEii,
            &indicesEij, &indicesEik };
        vector<RunlineSumEH_i> runlines;
        CallbackSumEH_i callback(runlines);
        RLE::merge(arrays, 4, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new SumEH_i(Field::e(xyz), 
                Field::ee(xyz,xyz), Field::ee(xyz,j), Field::ee(xyz,k),
                runlines));
        }
    }
}
void ForwardGridOperations::
encodeEij_j(const VoxelizedPartition & voxelGrid)
{
    LOG << "Encoding off-diagonal E (method j)\n";
    for (int numerOrder = 0; numerOrder < 10; numerOrder++)
    for (int denomOrder = 0; denomOrder < 10; denomOrder++)
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    if (ii != jj)
    {
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(jj)),
            voxelGrid.extents().dimensions());
        
        SupportRegion3 suppEij(voxelGrid.fieldIndices().ee(ii, jj));
        
        SupportRegion3 permittivityMask(voxelGrid.inversePermittivity()
            .support(ii,jj, numerOrder, denomOrder));
        
        SupportRegion3 maskD;
        maskD = permittivityMask * calcCells;
        
        IndexArray3 permittivityIndices(
            voxelGrid.inversePermittivity().filter(
            ii, jj, numerOrder, denomOrder));
        IndexArray3 indicesD, indicesE, indicesEps, indicesAux(maskD);
        restriction(voxelGrid.fieldIndices().d(jj), maskD, indicesD);
        restriction(voxelGrid.fieldIndices().ee(ii,jj), maskD, indicesE);
        restriction(permittivityIndices, maskD, indicesEps);
        
        IndexArray3 const* arrays[] =
            { &indicesE, &indicesD, &indicesAux, &indicesEps };
        vector<RunlineAnisotropicEH_j> runlines;
        CallbackAnisotropicEH_j callback(runlines);
        RLE::merge(arrays, 4, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new UpdateAnisotropicEH_j(
                Field::ee(ii,jj), Field::d(jj),
                denomOrder,
                numerOrder,
                runlines,
                voxelGrid.inversePermittivity()
                    .filter(ii, jj, numerOrder, denomOrder).values(),
                Constants::eps0));
        }
    }
    
    LOG << "Encoding summation of E (method j)\n";
    for (int xyz = 0; xyz < 3; xyz++)
    {
        int j = (xyz+1)%3, k = (xyz+2)%3;
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 suppEii(voxelGrid.fieldIndices().ee(xyz,xyz));
        
        SupportRegion3 sumRegion = suppEii*calcCells;
        
        SupportRegion3 sumRegionI, sumRegionJ, sumRegionK,
            sumRegionIJ, sumRegionIK;
        
        Vector3i di(Vector3i::unit(xyz)),
            dj(-Vector3i::unit(j)), dk(-Vector3i::unit(k)),
            dij(di+dj), dik(di+dk);
        
        translate(sumRegion, sumRegionI, di[0], di[1], di[2]);
        translate(sumRegion, sumRegionJ, dj[0], dj[1], dj[2]);
        translate(sumRegion, sumRegionK, dk[0], dk[1], dk[2]);
        translate(sumRegion, sumRegionIJ, dij[0], dij[1], dij[2]);
        translate(sumRegion, sumRegionIK, dik[0], dik[1], dik[2]);
        
        IndexArray3 indicesEi, indicesEii,
            indicesEij0, indicesEij1, indicesEij2, indicesEij3,
            indicesEik0, indicesEik1, indicesEik2, indicesEik3;
        
        restriction(voxelGrid.fieldIndices().e(xyz), sumRegion, indicesEi);
        restriction(voxelGrid.fieldIndices().ee(xyz,xyz), sumRegion, indicesEii);
        restriction(voxelGrid.fieldIndices().ee(xyz,j), sumRegion, indicesEij0);
        restriction(voxelGrid.fieldIndices().ee(xyz,j), sumRegionI, indicesEij1);
        restriction(voxelGrid.fieldIndices().ee(xyz,j), sumRegionJ, indicesEij2);
        restriction(voxelGrid.fieldIndices().ee(xyz,j), sumRegionIJ, indicesEij3);
        restriction(voxelGrid.fieldIndices().ee(xyz,k), sumRegion, indicesEik0);
        restriction(voxelGrid.fieldIndices().ee(xyz,k), sumRegionI, indicesEik1);
        restriction(voxelGrid.fieldIndices().ee(xyz,k), sumRegionK, indicesEik2);
        restriction(voxelGrid.fieldIndices().ee(xyz,k), sumRegionIK, indicesEik3);
        
        IndexArray3 const* arrays[] = { &indicesEi, &indicesEii,
            &indicesEij0, &indicesEij1, &indicesEij2, &indicesEij3,
            &indicesEik0, &indicesEik1, &indicesEik2, &indicesEik3 };
        vector<RunlineSumEH_j> runlines;
        CallbackSumEH_j callback(runlines);
        RLE::merge(arrays, 10, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new SumEH_j(Field::e(xyz), 
                Field::ee(xyz,xyz), Field::ee(xyz,j), Field::ee(xyz,k),
                runlines));
        }
    }
}

void ForwardGridOperations::
encodeB(const VoxelizedPartition & voxelGrid)
{
    Matrix3i m = voxelGrid.projectionMatrix();
    
    // Without current source
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (false == voxelGrid.modeUsesField(Field::b(xyz)))
        {
            continue;
        }
        
        Vector3i shiftY = m*Vector3i::unit( (xyz+1)%3 );
        Vector3i shiftZ = m*Vector3i::unit( (xyz+2)%3 );
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantH(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 maskB(voxelGrid.fieldIndices().b(xyz));
        maskB -= voxelGrid.fieldIndices().m(xyz);
        maskB *= calcCells;
        SupportRegion3 maskEjPlus, maskEjMinus, maskEkPlus, maskEkMinus;
        maskEjMinus = maskB;
        translate(maskB, maskEjPlus, shiftZ[0], shiftZ[1], shiftZ[2]);
        maskEkMinus = maskB;
        translate(maskB, maskEkPlus, shiftY[0], shiftY[1], shiftY[2]);
        
        IndexArray3 indicesB, indicesEjPlus, indicesEjMinus,
            indicesEkPlus, indicesEkMinus;
        restriction(voxelGrid.fieldIndices().b(xyz), maskB, indicesB);
        restriction(voxelGrid.fieldIndices().e((xyz+1)%3), maskEjMinus,
            indicesEjMinus);
        restriction(voxelGrid.fieldIndices().e((xyz+1)%3), maskEjPlus,
            indicesEjPlus);
        restriction(voxelGrid.fieldIndices().e((xyz+2)%3), maskEkMinus,
            indicesEkMinus);
        restriction(voxelGrid.fieldIndices().e((xyz+2)%3), maskEkPlus,
            indicesEkPlus);
        
        IndexArray3 const* arrays[] = { &indicesB,
            &indicesEjMinus, &indicesEjPlus,
            &indicesEkMinus, &indicesEkPlus };
        
        vector<RunlineDB> runlines;
        CallbackDB callback(runlines);
        RLE::merge(arrays, 5, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateH(new ForwardDB(
                Field::b(xyz), Field::e((xyz+1)%3), Field::e((xyz+2)%3),
                -gridDescription()->dt()/gridDescription()->dxyz()[(xyz+1)%3],
                gridDescription()->dt()/gridDescription()->dxyz()[(xyz+2)%3],
                runlines));
        }
    }
    
    // With current source M
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (false == voxelGrid.modeUsesField(Field::b(xyz)))
        {
            continue;
        }
        
        Vector3i shiftY = m*Vector3i::unit( (xyz+1)%3 );
        Vector3i shiftZ = m*Vector3i::unit( (xyz+2)%3 );
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantH(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 maskB(voxelGrid.fieldIndices().b(xyz));
        maskB *= voxelGrid.fieldIndices().m(xyz);
        maskB *= calcCells;
        SupportRegion3 maskEjPlus, maskEjMinus, maskEkPlus, maskEkMinus;
        maskEjMinus = maskB;
        translate(maskB, maskEjPlus, shiftZ[0], shiftZ[1], shiftZ[2]);
        maskEkMinus = maskB;
        translate(maskB, maskEkPlus, shiftY[0], shiftY[1], shiftY[2]);
        
        IndexArray3 indicesB, indicesEjPlus, indicesEjMinus,
            indicesEkPlus, indicesEkMinus, indicesM;
        restriction(voxelGrid.fieldIndices().b(xyz), maskB, indicesB);
        restriction(voxelGrid.fieldIndices().e((xyz+1)%3), maskEjMinus,
            indicesEjMinus);
        restriction(voxelGrid.fieldIndices().e((xyz+1)%3), maskEjPlus,
            indicesEjPlus);
        restriction(voxelGrid.fieldIndices().e((xyz+2)%3), maskEkMinus,
            indicesEkMinus);
        restriction(voxelGrid.fieldIndices().e((xyz+2)%3), maskEkPlus,
            indicesEkPlus);
        restriction(voxelGrid.fieldIndices().m(xyz), maskB, indicesM);
        
        IndexArray3 const* arrays[] = { &indicesB,
            &indicesEjMinus, &indicesEjPlus,
            &indicesEkMinus, &indicesEkPlus,
            &indicesM };
        
        vector<RunlineDB_JM> runlines;
        CallbackDB_JM callback(runlines);
        RLE::merge(arrays, 6, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateH(new ForwardDB_JM(
                Field::b(xyz), Field::e((xyz+1)%3), Field::e((xyz+2)%3),
                Field::m(xyz),
                -gridDescription()->dt()/gridDescription()->dxyz()[(xyz+1)%3],
                gridDescription()->dt()/gridDescription()->dxyz()[(xyz+2)%3],
                -gridDescription()->dt(),
                runlines));
        }
    }
}


void ForwardGridOperations::
encodeNondispersiveH(const VoxelizedPartition & voxelGrid)
{
    Matrix3i m = voxelGrid.projectionMatrix();
    
    // Without current sources
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (false == voxelGrid.modeUsesField(Field::h(xyz)))
        {
            continue;
        }
        
        Vector3i shiftY = m*Vector3i::unit( (xyz+1)%3 );
        Vector3i shiftZ = m*Vector3i::unit( (xyz+2)%3 );
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantH(xyz)),
            voxelGrid.extents().dimensions());
        
        SupportRegion3 maskNondispersive(voxelGrid.inversePermeability()
            .support(xyz, xyz, 0, 0));
        
        SupportRegion3 maskH(voxelGrid.inversePermeability()
            .support(xyz, xyz, 0, 0) - voxelGrid.fieldIndices().b(xyz));
        
        maskH -= voxelGrid.fieldIndices().m(xyz);
        maskH *= maskNondispersive;
        maskH *= calcCells;
        
        SupportRegion3 maskEjPlus, maskEjMinus, maskEkPlus, maskEkMinus;
        maskEjMinus = maskH;
        translate(maskH, maskEjPlus, shiftZ[0], shiftZ[1], shiftZ[2]);
        maskEkMinus = maskH;
        translate(maskH, maskEkPlus, shiftY[0], shiftY[1], shiftY[2]);
        
        IndexArray3 permeabilityIndices(voxelGrid.inversePermeability().filter(
            xyz, xyz, 0, 0));
        
        IndexArray3 indicesH, indicesEjPlus, indicesEjMinus,
            indicesEkPlus, indicesEkMinus, indicesMu;
        restriction(voxelGrid.fieldIndices().h(xyz), maskH, indicesH);
        restriction(voxelGrid.fieldIndices().e((xyz+1)%3), maskEjMinus,
            indicesEjMinus);
        restriction(voxelGrid.fieldIndices().e((xyz+1)%3), maskEjPlus,
            indicesEjPlus);
        restriction(voxelGrid.fieldIndices().e((xyz+2)%3), maskEkMinus,
            indicesEkMinus);
        restriction(voxelGrid.fieldIndices().e((xyz+2)%3), maskEkPlus,
            indicesEkPlus);
        restriction(permeabilityIndices, maskH, indicesMu);
        
//        cout << "H: " << indicesH << "\n";
//        cout << "Ej1: " << indicesEjMinus << "\n";
//        cout << "Ej2: " << indicesEjPlus << "\n";
//        cout << "Ek1: " << indicesEkMinus << "\n";
//        cout << "Ek2: " << indicesEkPlus << "\n";
        
        IndexArray3 const* arrays[] = { &indicesH,
            &indicesEjMinus, &indicesEjPlus,
            &indicesEkMinus, &indicesEkPlus,
            &indicesMu};
        
        vector<RunlineNondispersiveEH> runlines;
        CallbackNondispersiveEH callback(runlines);
        RLE::merge(arrays, 6, callback);
        
//        for (int rr = 0; rr < runlines.size(); rr++)
//            cout << runlines[rr] << "\n";
        
        if (runlines.size() > 0)
        {
            pushUpdateH(new ForwardNondispersiveEH(
                Field::h(xyz), Field::e((xyz+1)%3), Field::e((xyz+2)%3),
                -gridDescription()->dt()/gridDescription()->dxyz()[(xyz+1)%3],
                gridDescription()->dt()/gridDescription()->dxyz()[(xyz+2)%3],
                runlines,
                voxelGrid.inversePermeability().filter(xyz,xyz,0,0).values(),
                Constants::mu0));
        }
    }
    
    // With current sources
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (false == voxelGrid.modeUsesField(Field::h(xyz)))
        {
            continue;
        }
        
        Vector3i shiftY = m*Vector3i::unit( (xyz+1)%3 );
        Vector3i shiftZ = m*Vector3i::unit( (xyz+2)%3 );
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantH(xyz)),
            voxelGrid.extents().dimensions());
        
        SupportRegion3 maskNondispersive(voxelGrid.inversePermeability()
            .support(xyz, xyz, 0, 0));
        
        SupportRegion3 maskH(voxelGrid.inversePermeability()
            .support(xyz, xyz, 0, 0) - voxelGrid.fieldIndices().b(xyz));
        
        maskH *= voxelGrid.fieldIndices().m(xyz);
        maskH *= maskNondispersive;
        maskH *= calcCells;
        
        SupportRegion3 maskEjPlus, maskEjMinus, maskEkPlus, maskEkMinus;
        maskEjMinus = maskH;
        translate(maskH, maskEjPlus, shiftZ[0], shiftZ[1], shiftZ[2]);
        maskEkMinus = maskH;
        translate(maskH, maskEkPlus, shiftY[0], shiftY[1], shiftY[2]);
        
        IndexArray3 permeabilityIndices(voxelGrid.inversePermeability().filter(
            xyz, xyz, 0, 0));
        
        IndexArray3 indicesH, indicesEjPlus, indicesEjMinus,
            indicesEkPlus, indicesEkMinus, indicesMu, indicesM;
        restriction(voxelGrid.fieldIndices().h(xyz), maskH, indicesH);
        restriction(voxelGrid.fieldIndices().e((xyz+1)%3), maskEjMinus,
            indicesEjMinus);
        restriction(voxelGrid.fieldIndices().e((xyz+1)%3), maskEjPlus,
            indicesEjPlus);
        restriction(voxelGrid.fieldIndices().e((xyz+2)%3), maskEkMinus,
            indicesEkMinus);
        restriction(voxelGrid.fieldIndices().e((xyz+2)%3), maskEkPlus,
            indicesEkPlus);
        restriction(permeabilityIndices, maskH, indicesMu);
        restriction(voxelGrid.fieldIndices().m(xyz), maskH, indicesM);
        
//        cout << "H: " << indicesH << "\n";
//        cout << "Ej1: " << indicesEjMinus << "\n";
//        cout << "Ej2: " << indicesEjPlus << "\n";
//        cout << "Ek1: " << indicesEkMinus << "\n";
//        cout << "Ek2: " << indicesEkPlus << "\n";
//        cout << "M: " << indicesM << "\n";
        
        IndexArray3 const* arrays[] = { &indicesH,
            &indicesEjMinus, &indicesEjPlus,
            &indicesEkMinus, &indicesEkPlus,
            &indicesMu, &indicesM};
        
        vector<RunlineNondispersiveEH_JM> runlines;
        CallbackNondispersiveEH_JM callback(runlines);
        RLE::merge(arrays, 7, callback);
        
//        for (int rr = 0; rr < runlines.size(); rr++)
//            cout << runlines[rr] << "\n";
        
        if (runlines.size() > 0)
        {
            pushUpdateH(new ForwardNondispersiveEH_JM(
                Field::h(xyz), Field::e((xyz+1)%3), Field::e((xyz+2)%3),
                Field::m(xyz),
                -gridDescription()->dt()/gridDescription()->dxyz()[(xyz+1)%3],
                gridDescription()->dt()/gridDescription()->dxyz()[(xyz+2)%3],
                -gridDescription()->dt(),
                runlines,
                voxelGrid.inversePermeability().filter(xyz,xyz,0,0).values(),
                Constants::mu0));
        }
    }
}



void ForwardGridOperations::
encodeH(const VoxelizedPartition & voxelGrid)
{
    LOG << "Encoding diagonal H, crudely\n";
    for (int numerOrder = 0; numerOrder < 10; numerOrder++)
    for (int denomOrder = 0; denomOrder < 10; denomOrder++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (false == voxelGrid.modeUsesField(Field::h(xyz)))
        {
            continue;
        }
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantH(xyz)),
            voxelGrid.extents().dimensions());
        
        SupportRegion3 suppHii(voxelGrid.fieldIndices().hh(xyz, xyz));
        
        SupportRegion3 suppB(voxelGrid.fieldIndices().b(xyz));
        SupportRegion3 maskB(voxelGrid.inversePermeability()
            .support(xyz, xyz, numerOrder, denomOrder));
        maskB *= calcCells;
        maskB *= suppB;
        maskB -= suppHii;
        
        IndexArray3 permeabilityIndices(voxelGrid
            .inversePermeability().filter(xyz, xyz, numerOrder, denomOrder));
        IndexArray3 indicesB, indicesH, indicesMu, indicesAux(maskB);
        restriction(voxelGrid.fieldIndices().b(xyz), maskB, indicesB);
        restriction(voxelGrid.fieldIndices().h(xyz), maskB, indicesH);
        restriction(permeabilityIndices, maskB, indicesMu);
        
        IndexArray3 const* arrays[] =
            { &indicesH, &indicesB, &indicesAux, &indicesMu };
        vector<RunlineEH> runlines;
        CallbackEH callback(runlines);
        RLE::merge(arrays, 4, callback);
        
        if (runlines.size() > 0)
        {
            assert(numerOrder == denomOrder);
            int order = numerOrder;
            pushUpdateH(FwdEHFactory::newForwardEH(
                Field::h(xyz),
                Field::b(xyz),
                order,
                runlines,
                voxelGrid.inversePermeability()
                    .filter(xyz, xyz, numerOrder, denomOrder).values(),
                Constants::mu0));
            
//            pushUpdateH(new ForwardEH(
//                Field::h(xyz), Field::b(xyz),
//                numerOrder,
//                denomOrder,
//                runlines,
//                voxelGrid.inversePermeability()
//                    .filter(xyz, xyz, numerOrder, denomOrder).values(),
//                Constants::mu0));
        }
    }
    /*
    LOG << "Encoding summation of H\n";
    for (int xyz = 0; xyz < 3; xyz++)
    {
        int j = (xyz+1)%3, k = (xyz+2)%3;
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantH(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 suppHii(voxelGrid.fieldIndices().hh(xyz,xyz));
        
        SupportRegion3 sumRegion = suppHii*calcCells;
        
        IndexArray3 indicesHi, indicesHii, indicesHij, indicesHik;
        restriction(voxelGrid.fieldIndices().h(xyz), sumRegion, indicesHi);
        restriction(voxelGrid.fieldIndices().hh(xyz,xyz), sumRegion, indicesHii);
        restriction(voxelGrid.fieldIndices().hh(xyz,j), sumRegion, indicesHij);
        restriction(voxelGrid.fieldIndices().hh(xyz,k), sumRegion, indicesHik);
        
        IndexArray3 const* arrays[] = { &indicesHi, &indicesHii,
            &indicesHij, &indicesHik };
        vector<RunlineSumEH_0> runlines;
        CallbackSumEH_0 callback(runlines);
        RLE::merge(arrays, 4, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateH(new SumEH_0(Field::h(xyz),
                Field::hh(xyz,xyz), Field::hh(xyz,j), Field::hh(xyz,k),
                runlines));
        }
    }
    */
}

void ForwardGridOperations::
encodePMLJ(const VoxelizedPartition & voxelGrid)
{
    Matrix3i m = voxelGrid.projectionMatrix();
    
    LOG << "Encoding PML J updates\n";
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (false == voxelGrid.modeUsesField(Field::j(xyz)))
        {
            continue;
        }
        
        for (int mm = 1; mm <= 2; mm++)
        {
            Precision::Float leadingSign = (mm == 1 ? -1.0f : 1.0f);
            int absorbXYZ = (xyz+mm)%3;
            int hXYZ = (xyz+3-mm)%3;
            assert(absorbXYZ != xyz);
            assert(hXYZ != xyz);
            assert(absorbXYZ != hXYZ);
            
            if (false == voxelGrid.modeUsesField(Field::h(hXYZ)))
            {
                continue;
            }
            
            Vector3i shift = -m*Vector3i::unit(absorbXYZ);
            
            SupportRegion3 maskPML(voxelGrid.pmlIndices().e(xyz, absorbXYZ));
            SupportRegion3 maskHPlus = maskPML;
            SupportRegion3 maskHMinus;
            translate(maskPML, maskHMinus, shift[0], shift[1], shift[2]);
            
            IndexArray3 indicesPML(maskPML); // could get from pmlIndices()...
            IndexArray3 indicesJ, indicesHPlus, indicesHMinus, indicesConsts;
            restriction(voxelGrid.fieldIndices().j(xyz), maskPML, indicesJ);
            restriction(voxelGrid.fieldIndices().h(hXYZ), maskHMinus,
                indicesHMinus);
            restriction(voxelGrid.fieldIndices().h(hXYZ), maskHPlus,
                indicesHPlus);
            indicesConsts.segment(0);
            
//            IndexArray3 rawIndices(voxelGrid.pmlE(xyz, absorbXYZ));
//            restriction(rawIndices, maskPML, indicesConsts);
            restriction(
                IndexArray3(voxelGrid.pmlE(xyz, absorbXYZ)),
                maskPML, indicesConsts);
                
//            cerr << "\nPML J " << char('x'+xyz) << char('x'+mm) << ":\n";
//            cerr << "Mask PML " << maskPML << "\n";
//            cerr << "indicesPML " << indicesPML << "\n";
//            cerr << "indicesJ " << indicesJ << "\n";
//            cerr << "indicesH+: " << indicesHPlus << "\n";
//            cerr << "indicesH-: " << indicesHMinus << "\n";
//            cerr << "indices consts: " << indicesConsts << "\n";
            
            IndexArray3 const* arrays[] = { &indicesJ,
                &indicesHPlus, &indicesHMinus,
                &indicesPML,
                &indicesConsts };
            
            vector<RunlineJMPML> runlines;
            CallbackJMPML callback(runlines);
            RLE::merge(arrays, 5, callback);
            
            if (runlines.size() > 0)
            {
                pushUpdateE(new UpdateJMPML(
                    Field::j(xyz), Field::h(hXYZ),
                    leadingSign,
                    runlines,
                    voxelGrid.pmlE(xyz, absorbXYZ).values(),
                    voxelGrid.gridDescription()->dt(),
                    voxelGrid.gridDescription()->dxyz()[absorbXYZ]));
            }
        }
    }
}

void ForwardGridOperations::
encodePMLM(const VoxelizedPartition & voxelGrid)
{
    Matrix3i m = voxelGrid.projectionMatrix();
    
    LOG << "Encoding PML M updates\n";
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (false == voxelGrid.modeUsesField(Field::m(xyz)))
        {
            continue;
        }
        
        for (int mm = 1; mm <= 2; mm++)
        {
            Precision::Float leadingSign = (mm == 1 ? 1.0f : -1.0f);
            int absorbXYZ = (xyz+mm)%3;
            int eXYZ = (xyz+3-mm)%3;
            assert(absorbXYZ != xyz);
            assert(eXYZ != xyz);
            assert(absorbXYZ != eXYZ);
            
            if (false == voxelGrid.modeUsesField(Field::e(eXYZ)))
            {
                continue;
            }
            
            Vector3i shift = m*Vector3i::unit(absorbXYZ);
            
            SupportRegion3 maskPML(voxelGrid.pmlIndices().h(xyz, absorbXYZ));
            SupportRegion3 maskEPlus;
            translate(maskPML, maskEPlus, shift[0], shift[1], shift[2]);
            SupportRegion3 maskEMinus = maskPML;
            
            IndexArray3 indicesPML(maskPML); // could get from pmlIndices()...
            IndexArray3 indicesM, indicesEPlus, indicesEMinus, indicesConsts;
            restriction(voxelGrid.fieldIndices().m(xyz), maskPML, indicesM);
            restriction(voxelGrid.fieldIndices().e(eXYZ), maskEMinus,
                indicesEMinus);
            restriction(voxelGrid.fieldIndices().e(eXYZ), maskEPlus,
                indicesEPlus);
            restriction(
                IndexArray3(voxelGrid.pmlH(xyz, absorbXYZ)),
                maskPML, indicesConsts);
            
//            cerr << "\nPML M " << char('x'+xyz) << char('x'+mm) << ":\n";
//            cerr << "Mask PML " << maskPML << "\n";
//            cerr << "indicesPML " << indicesPML << "\n";
//            cerr << "indicesM " << indicesM << "\n";
//            cerr << "indicesE: " << indicesEPlus << "\n";
//            cerr << indicesEMinus << "\n";
//            cerr << "consts: " << indicesConsts << "\n";
//            cerr << "const vals: " << voxelGrid.pmlH(xyz, absorbXYZ) << "\n";
            
            IndexArray3 const* arrays[] = { &indicesM,
                &indicesEPlus, &indicesEMinus,
                &indicesPML, &indicesConsts };
            
            vector<RunlineJMPML> runlines;
            CallbackJMPML callback(runlines);
            RLE::merge(arrays, 5, callback);
            
            if (runlines.size() > 0)
            {
                pushUpdateH(new UpdateJMPML(
                    Field::m(xyz), Field::e(eXYZ),
                    leadingSign,
                    runlines,
                    voxelGrid.pmlH(xyz, absorbXYZ).values(),
                    voxelGrid.gridDescription()->dt(),
                    voxelGrid.gridDescription()->dxyz()[absorbXYZ]));
            }
        }
    }
}

