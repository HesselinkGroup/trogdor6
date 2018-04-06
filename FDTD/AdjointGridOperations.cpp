/*
 *  AdjointGridOperations.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 7/7/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "AdjointGridOperations.h"

#include "Log.h"
#include "rle/MergeRunlines.h"
#include "Support.h"
#include "YeeUtilities.h"
#include "PhysicalConstants.h"
#include "TimeWrapper.h"
#include "RLEOperations.h"
#include "UserPreferences.h"

#include "EncodeOutput.h"

#include "AdjointRunlineCallbacks.h"

#include "UniversalRunlineCallbacks.h"

using namespace std;
using namespace RLE;
using namespace YeeUtilities;
using namespace TimeWrapper;

AdjointGridOperations::
AdjointGridOperations(GridDescPtr gridDesc,
    const map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids) :
    GridOperations(gridDesc, voxelizedGrids)
{
    VoxelizedPartitionPtr voxelGrid = voxelizedGrids.find(gridDesc)->second;
    
    if (voxelGrid->gridType() == VoxelizedPartition::kUserGrid)
    {
        LOG << "Add adjoint outputs\n";
        addBoundaryOutputs(*gridDesc, *voxelGrid);
    }
    
    LOG << "Preparing outputs\n";
    // All outputs are potentially filters
    // Encode what to write?  How hard can this be?
    encodeOutputs(*gridDesc, *voxelGrid);
    
    LOG << "NOT preparing sources (what are they?)\n";
    
    LOG << "Preparing current sources\n";
    // Same deal.  They write to J.
    encodeCurrentSources(*gridDesc, *voxelGrid);
    
    LOG << "Preparing periodic boundaries\n";
    encodePeriodicBCs(*voxelGrid);
    
    LOG << "Preparing PML\n";
    encodePMLJ(*voxelGrid);
    
    LOG << "Preparing adjoint Huygens surfaces\n";
    // Use the support of the WRITE region as a mask in the READ grid.
    // Index that, PUNK!  wahaha
    // Anyway this includes E and H, both.
    encodeHuygens(*gridDesc, *voxelGrid, voxelizedGrids);
    
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
    encodeE(*voxelGrid);
    
    LOG << "Preparing D updates\n";
    // Separate D with J and D without J
    // Runlength encode for each of x y z
    encodeD(*voxelGrid);
    
    
    LOG << "Preparing PML M updates\n";
    encodePMLM(*voxelGrid);
    
    LOG << "Preparing H updates\n";
    // Diagonal B -> H
    // Off-diagonal: average B, filter, add to diagonal H
    //
    // Encode:
    //  - diagonal B to H, same for all diagonals
    //  - Off-diagonal B to H, which requires some averaging of B
    //  - Summation of anisotropic H to diagonal H (SAME AS B TO H?)
    encodeH(*voxelGrid);
    
    LOG << "Preparing B updates\n";
    // Separate B with and without M
    // Runlength encode for each of x y z
    encodeB(*voxelGrid);
    
    initPerformance();
}

void AdjointGridOperations::
initPerformance()
{
    mTimeClearJM = 0.0;
    mTimeE = vector<double>(mUpdatesE.size(), 0.0);
    mTimeH = vector<double>(mUpdatesH.size(), 0.0);
    mTimeSourceJM = vector<double>(mSourceJM.size(), 0.0);
    mTimeBCD = vector<double>(mPeriodicBCD.size(), 0.0);
    mTimeBCB = vector<double>(mPeriodicBCB.size(), 0.0);
    mTimeOutputs = vector<double>(mOutputs.size(), 0.0);
}

void AdjointGridOperations::
printPerformance(double numT) const
{
}

void AdjointGridOperations::
handleSetPointers(map<int, Pointer<GridFields> > & fields)
{
    assert(fields.count(gridDescription()->id()));
    GridFields & gridFields = *fields[gridDescription()->id()];
    
    for (int nn = 0; nn < mSourceJM.size(); nn++)
        mSourceJM[nn]->setPointers(gridFields, fields);

    for (int nn = 0; nn < mOutputs.size(); nn++)
        mOutputs[nn]->setPointers(gridFields);
    
    for (int nn = 0; nn < mUpdatesE.size(); nn++)
        mUpdatesE[nn]->setPointers(gridFields, fields);
    for (int nn = 0; nn < mUpdatesH.size(); nn++)
        mUpdatesH[nn]->setPointers(gridFields, fields);
    
    for (int nn = 0; nn < mPeriodicBCD.size(); nn++)
        mPeriodicBCD[nn]->setPointers(gridFields, fields);
    for (int nn = 0; nn < mPeriodicBCB.size(); nn++)
        mPeriodicBCB[nn]->setPointers(gridFields, fields);
}

void AdjointGridOperations::
allocate()
{
    for (int nn = 0; nn < mUpdatesE.size(); nn++)
        mUpdatesE[nn]->allocate();
    for (int nn = 0; nn < mUpdatesH.size(); nn++)
        mUpdatesH[nn]->allocate();
    for (int nn = 0; nn < mSourceJM.size(); nn++)
        mSourceJM[nn]->allocate();
}

void AdjointGridOperations::
firstHalfStep(long timestep, Precision::Float dt)
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
    
    for (int nn = 0; nn < mPeriodicBCB.size(); nn++)
    {
        t0 = now();
        mPeriodicBCB[nn]->apply(timestep, dt);
        mTimeH[nn] += elapsedMicroseconds(t0, now());
    }
}

void AdjointGridOperations::
secondHalfStep(long timestep, Precision::Float dt)
{
    clearJ();
    updateJSource(timestep, dt);
    
    TimePoint t0;
    
    for (int nn = 0; nn < mUpdatesE.size(); nn++)
    {
        t0 = now();
        mUpdatesE[nn]->apply(timestep, dt);
        mTimeE[nn] += elapsedMicroseconds(t0, now());
    }
    
    for (int nn = 0; nn < mPeriodicBCD.size(); nn++)
    {
        t0 = now();
        mPeriodicBCD[nn]->apply(timestep, dt);
        mTimeBCD[nn] += elapsedMicroseconds(t0, now());
    }
}

void AdjointGridOperations::
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

void AdjointGridOperations::
clearM()
{
    TimePoint t0 = now();
    fields().clearM();
    mTimeClearJM += elapsedMicroseconds(t0, now());
}

void AdjointGridOperations::
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

void AdjointGridOperations::
clearJ()
{
    TimePoint t0 = now();
    fields().clearJ();
    mTimeClearJM += elapsedMicroseconds(t0, now());
}

void AdjointGridOperations::
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

long AdjointGridOperations::
bytes() const
{
    long totalBytes = 0;
    
    for (int nn = 0; nn < mUpdatesE.size(); nn++)
        totalBytes += mUpdatesE[nn]->bytes();
    for (int nn = 0; nn < mUpdatesH.size(); nn++)
        totalBytes += mUpdatesH[nn]->bytes();
    
    return totalBytes;
}

void AdjointGridOperations::
printRunlines() const
{
    cout << "E fields:\n";
    for (int nn = 0; nn < mUpdatesE.size(); nn++)
    {
        cout << mUpdatesE[nn]->name() << ":\n";
        mUpdatesE[nn]->printRunlines(cout);
    }
    cout << "H fields:\n";
    for (int nn = 0; nn < mUpdatesH.size(); nn++)
    {
        cout << mUpdatesH[nn]->name() << ":\n";
        mUpdatesH[nn]->printRunlines(cout);
    }
//    for (int xyz = 0; xyz < 3; xyz++)
//    {
//        for (int nn = 0; nn < mUpdateD[xyz].size(); nn++)
//        {
//            cout << "\tUpdate D" << char('x'+xyz) << ":\n";
//            mUpdateD[xyz][0].printRunlines(cout);
//        }
//    }
//    
////    for (int xyz = 0; xyz < 3; xyz++)
////    {
////        cout << "\tUpdate DJ" << char('x'+xyz) << ":\n";
////        mUpdateDJ[xyz][0].printRunlines(cout);
////        assert(mUpdateDJ[xyz].size() == 1);
////    }
//    
//    for (int xyz = 0; xyz < 3; xyz++)
//    {
//        for (int nn = 0; nn < mUpdateE[xyz].size(); nn++)
//        {
//            cout << nn << "\tUpdate E" << char('x'+xyz) << ":\n";
//            mUpdateE[xyz][nn].printRunlines(cout);
//        }
//    }
//    
//    for (int xyz = 0; xyz < 3; xyz++)
//    {
//        for (int nn = 0; nn < mUpdateB[xyz].size(); nn++)
//        {
//            cout << "\tUpdate B" << char('x'+xyz) << ":\n";
//            mUpdateB[xyz][0].printRunlines(cout);
//        }
//    }
//    
////    for (int xyz = 0; xyz < 3; xyz++)
////    {
////        cout << "\tUpdate BM" << char('x'+xyz) << ":\n";
////        mUpdateBM[xyz][0].printRunlines(cout);
////        assert(mUpdateBM[xyz].size() == 1);
////    }
//    
//    for (int xyz = 0; xyz < 3; xyz++)
//    {
//        for (int nn = 0; nn < mUpdateH[xyz].size(); nn++)
//        {
//            cout << nn << "\tUpdate H" << char('x'+xyz) << ":\n";
//            mUpdateH[xyz][nn].printRunlines(cout);
//        }
//    }
//    
//    for (int xyz = 0; xyz < 3; xyz++)
//    {
//        for (int absorb = 0; absorb <= 1; absorb++)
//        {
//            cout << "PML J" << char('x'+xyz) << " along "
//                << char('x'+ (absorb+xyz+1)%3 )
//                << ":\n";
//            mUpdateJPML[xyz][absorb].printRunlines(cout);
//        }
//    }
//    
//    for (int xyz = 0; xyz < 3; xyz++)
//    {
//        for (int absorb = 0; absorb < 1; absorb++)
//        {
//            cout << "PML M" << char('x'+xyz) << " along "
//                << char('x'+ (absorb+xyz+1)%3)
//                << ":\n";
//            mUpdateMPML[xyz][absorb].printRunlines(cout);
//        }
//    }
}


void AdjointGridOperations::
addBoundaryOutputs(GridDescription & gridDesc,
    const VoxelizedPartition & voxelGrid)
{
    // E, D
    vector<Duration> durations(1, Duration(0, gridDesc.numTimesteps()-1));
    
    if (false == UserPreferences::defines("sensitivity"))
        return;
    
    for (int ii = 0; ii < 3; ii++)
    {
        int jj = ii;
    //for (int jj = 0; jj < 3; jj++)
    if (voxelGrid.sensitiveCells(ii,jj).numRuns() > 0)
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
        fields << "e" << char('x'+ii) << char('x'+jj);
//        if (ii == jj)
//            fields << " d" << char('x'+ii);
        string filename;
        filename = string("boundary_adjoint_e") + char('x'+ii) + char('x'+jj);
        OutputDescPtr out(new OutputDescription(fields.str(), filename,
            regions, durations));
        gridDesc.pushOutput(out);
    }
    }
}

void AdjointGridOperations::
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

void AdjointGridOperations::
encodeCurrentSources(const GridDescription & gridDesc,
    const VoxelizedPartition & voxelGrid)
{
    for (int nn = 0; nn < gridDesc.currentSources().size(); nn++)
    {
        LOG << "Current source " << nn << "\n";
        const CurrentSourceDescription & outputDesc =
            *(gridDesc.currentSources())[nn];
        
        Pointer<SourceJM> src(new SourceJM(outputDesc));
        mSourceJM.push_back(src);
        
        // J sources
        for (int xyz = 0; xyz < 3; xyz++)
        {
            if (outputDesc.sourceCurrents().whichJ()[xyz])
            {
                SupportRegion3 supp = support(outputDesc,
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
            if (outputDesc.sourceCurrents().whichJE()[xyz])
            {
                SupportRegion3 supp = support(outputDesc,
                    voxelGrid.extents().nodeYeeCells().num(),
                    octantE(xyz));
                IndexArray3 indicesJE;
                restriction(voxelGrid.fieldIndices().je(xyz), supp, indicesJE);
                
                IndexArray3 const* arrays[] = { &indicesJE };

                vector<RunlineSource> runlines;
                CallbackSource callback(runlines);
                RLE::merge(arrays, 1, callback);
                src->runlinesJE(xyz, runlines);
            }
        }
        
        // M sources
        for (int xyz = 0; xyz < 3; xyz++)
        {
            if (outputDesc.sourceCurrents().whichM()[xyz])
            {
                SupportRegion3 supp = support(outputDesc,
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
            if (outputDesc.sourceCurrents().whichMH()[xyz])
            {
                SupportRegion3 supp = support(outputDesc,
                    voxelGrid.extents().nodeYeeCells().num(),
                    octantH(xyz));
                IndexArray3 indicesMH;
                restriction(voxelGrid.fieldIndices().mh(xyz), supp, indicesMH);
                
                IndexArray3 const* arrays[] = { &indicesMH };
                
                vector<RunlineSource> runlines;
                CallbackSource callback(runlines);
                RLE::merge(arrays, 1, callback);
                src->runlinesMH(xyz, runlines);
            }
        }
    }
}
        
void AdjointGridOperations::
encodeHuygens(const GridDescription & gridDesc,
    const VoxelizedPartition & voxelGrid,
    const map<GridDescPtr, VoxelizedPartitionPtr> & voxelizedGrids)
{
    return;
    cerr << "This does not calculate the right thing...\n";
    for (int nn = 0; nn < gridDesc.huygensSurfaces().size(); nn++)
    {
        LOG << "Huygens surface #" << nn << ", source runs:\n";
        const HuygensSurfaceDescription & huyg = 
            *(gridDesc.huygensSurfaces()[nn]);
        
        // Encode E
        for (int xyz = 0; xyz < 3; xyz++)
        {
            SupportRegion3 supp = support(huyg,
                voxelGrid.extents().nodeYeeCells().num(),
                octantE(xyz));
            IndexArray3 indicesE;
            restriction(voxelGrid.fieldIndices().e(xyz), supp, indicesE);
            
            IndexArray3 const* arrays[] = { &indicesE };
            
            int numRuns;
            RLECallback callback(&numRuns);
            
            RLE::merge(arrays, 1, callback);
            if (numRuns > 0)
            {
                LOG << "Huygens E" << char('x'+xyz) << " has " << numRuns
                    << " source runs.\n";
            }
        }
        
        // Encode H
        for (int xyz = 0; xyz < 3; xyz++)
        {
            SupportRegion3 supp = support(huyg,
                voxelGrid.extents().nodeYeeCells().num(),
                octantH(xyz));
            IndexArray3 indicesH;
            restriction(voxelGrid.fieldIndices().h(xyz), supp, indicesH);
            
            IndexArray3 const* arrays[] = { &indicesH };
            
            int numRuns;
            RLECallback callback(&numRuns);
            
            RLE::merge(arrays, 1, callback);
            if (numRuns > 0)
            {
                LOG << "Huygens H" << char('x'+xyz) << " has " << numRuns
                    << " source runs.\n";
            }
        }
    }
}


void AdjointGridOperations::
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
            // Encode D:
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
            
            restriction(voxelGrid.fieldIndices().d(fieldXYZ), suppRead, indRead);
            restriction(voxelGrid.fieldIndices().d(fieldXYZ), suppWrite, indWrite);
            
//            cout << "Copy:\n" << indRead;
//            cout << "\nPaste:\n" << indWrite << "\n";
            
            IndexArray3 const* arrays[] = { &indRead, &indWrite };
            vector<RunlinePeriodicBC> runlines;
            CallbackPeriodicBC callback(runlines);
            RLE::merge(arrays, 2, callback);
            
            mPeriodicBCD.push_back(Pointer<PeriodicBC>(new PeriodicBC(
                Field::d(fieldXYZ),
                runlines)));
        }
        
        for (int fieldXYZ = 0; fieldXYZ < 3; fieldXYZ++)
        {
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
            
            restriction(voxelGrid.fieldIndices().b(fieldXYZ), suppRead, indRead);
            restriction(voxelGrid.fieldIndices().b(fieldXYZ), suppWrite, indWrite);
            
            IndexArray3 const* arrays[] = { &indRead, &indWrite };
            vector<RunlinePeriodicBC> runlines;
            CallbackPeriodicBC callback(runlines);
            RLE::merge(arrays, 2, callback);
            
            mPeriodicBCB.push_back(Pointer<PeriodicBC>(new PeriodicBC(
                Field::b(fieldXYZ),
                runlines)));
        }
    }
}


// TODO: Use curiously recurring template pattern to permit return-by-value
// for restriction(), translation(), etc.  (Will this work?)
void AdjointGridOperations::
encodeD(const VoxelizedPartition & voxelGrid)
{    
    // Without current sources
    int denomOrder = -1; // matches anything
    
    for (int numerOrder = 0; numerOrder < 10; numerOrder++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 maskD(voxelGrid.inversePermittivity()
            .support(xyz, xyz, numerOrder, denomOrder));
        maskD -= voxelGrid.fieldIndices().j(xyz);
        maskD *= calcCells;
        //LOG << "D" << char('x'+xyz) << " numer " << numerOrder
        //    << " has " << maskD.numRuns() << " runs.\n";
        IndexArray3 permittivityIndices(
            voxelGrid.inversePermittivity().filter(
            xyz, xyz, numerOrder, denomOrder));

        IndexArray3 indicesD, indicesE, indicesEps, indicesAux(maskD);
        restriction(voxelGrid.fieldIndices().d(xyz), maskD, indicesD);
        restriction(voxelGrid.fieldIndices().e(xyz), maskD, indicesE);
        restriction(permittivityIndices, maskD, indicesEps);
        
        IndexArray3 const* arrays[] = { &indicesD,
            &indicesE, &indicesAux, &indicesEps};
        
        vector<RunlineAdjointDB> runlines;
        CallbackAdjointDB callback(runlines);
        RLE::merge(arrays, 4, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new AdjointDB(
                Field::d(xyz), Field::e(xyz),
                numerOrder,
                runlines,
                voxelGrid.inversePermittivity()
                    .filter(xyz, xyz, numerOrder, denomOrder).values(),
                gridDescription()->dt() ));
        }
    }
    
    // With current sources
    
    for (int numerOrder = 0; numerOrder < 10; numerOrder++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
//        SupportRegion3 maskD(voxelGrid.permittivitySupport(
//            xyz, xyz, numerOrder, denomOrder));
        SupportRegion3 maskD(voxelGrid.inversePermittivity()
            .support(xyz, xyz, numerOrder, denomOrder));
        maskD *= voxelGrid.fieldIndices().j(xyz);
        maskD *= calcCells;
        //LOG << "D" << char('x'+xyz) << " with source, numer " << numerOrder
        //    << " has " << maskD.numRuns() << " runs.\n";
        IndexArray3 permittivityIndices(
            voxelGrid.inversePermittivity().filter(
            xyz, xyz, numerOrder, denomOrder));

        IndexArray3 indicesD, indicesE, indicesEps, indicesAux(maskD);
        restriction(voxelGrid.fieldIndices().d(xyz), maskD, indicesD);
        restriction(voxelGrid.fieldIndices().e(xyz), maskD, indicesE);
        restriction(permittivityIndices, maskD, indicesEps);
        
        IndexArray3 indicesJ;
        restriction(voxelGrid.fieldIndices().j(xyz), maskD, indicesJ);
        
        IndexArray3 const* arrays[] = { &indicesD,
            &indicesE, &indicesAux, &indicesEps, &indicesJ };
        
        vector<RunlineAdjointDB_JM> runlines;
        CallbackAdjointDB_JM callback(runlines);
        RLE::merge(arrays, 5, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new AdjointDB_JM(
                Field::d(xyz), Field::e(xyz), Field::j(xyz),
                numerOrder,
                -gridDescription()->dt(),
                runlines,
                voxelGrid.inversePermittivity()
                    .filter(xyz, xyz, numerOrder, denomOrder).values(),
                gridDescription()->dt() ));
//                Constants::eps0));
        }
    }
}

void AdjointGridOperations::
encodeE(const VoxelizedPartition & voxelGrid)
{
    Matrix3i m = voxelGrid.projectionMatrix();
    
    LOG << "Encoding diagonal E, crudely\n";
    int numerOrder = -1;
    
    for (int denomOrder = 0; denomOrder < 10; denomOrder++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Vector3i shiftY = -m*Vector3i::unit( (xyz+1)%3 );
        Vector3i shiftZ = -m*Vector3i::unit( (xyz+2)%3 );
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 maskD(voxelGrid.inversePermittivity()
            .support(xyz, xyz, numerOrder, denomOrder));
        maskD -= voxelGrid.fieldIndices().je(xyz);
        maskD *= calcCells;
        //LOG << "E" << char('x'+xyz) << " denom " << denomOrder
        //    << " has " << maskD.numRuns() << " runs.\n";
        
        SupportRegion3 maskBjPlus, maskBjMinus, maskBkPlus, maskBkMinus;
        maskBjPlus = maskD;
        translate(maskD, maskBjMinus, shiftZ[0], shiftZ[1], shiftZ[2]);
        maskBkPlus = maskD;
        translate(maskD, maskBkMinus, shiftY[0], shiftY[1], shiftY[2]);
        
        IndexArray3 indicesBjPlus, indicesBjMinus,
            indicesBkPlus, indicesBkMinus, indicesAux(maskD);
        restriction(voxelGrid.fieldIndices().b((xyz+1)%3), maskBjMinus,
            indicesBjMinus);
        restriction(voxelGrid.fieldIndices().b((xyz+1)%3), maskBjPlus,
            indicesBjPlus);
        restriction(voxelGrid.fieldIndices().b((xyz+2)%3), maskBkMinus,
            indicesBkMinus);
        restriction(voxelGrid.fieldIndices().b((xyz+2)%3), maskBkPlus,
            indicesBkPlus);
        
        IndexArray3 permittivityIndices(
            voxelGrid.inversePermittivity().filter(
            xyz, xyz, numerOrder, denomOrder));
        IndexArray3 indicesE, indicesEps;
        restriction(voxelGrid.fieldIndices().e(xyz), maskD, indicesE);
        restriction(permittivityIndices, maskD, indicesEps);
        
        IndexArray3 const* arrays[] = { &indicesE,
            &indicesBjMinus, &indicesBjPlus,
            &indicesBkMinus, &indicesBkPlus,
            &indicesAux, &indicesEps };
        
        vector<RunlineAdjointEH> runlines;
        CallbackAdjointEH callback(runlines);
        RLE::merge(arrays, 7, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new AdjointEH(
                Field::e(xyz), Field::b((xyz+1)%3), Field::b((xyz+2)%3),
                denomOrder,
                -1.0/gridDescription()->dxyz()[(xyz+1)%3],
                1.0/gridDescription()->dxyz()[(xyz+2)%3],
                runlines,
                voxelGrid.inversePermittivity()
                    .filter(xyz, xyz, numerOrder, denomOrder).values(),
                Constants::eps0));
        }
    }
    
    LOG << "Diagonal E with sources...\n";
    
    for (int denomOrder = 0; denomOrder < 10; denomOrder++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Vector3i shiftY = -m*Vector3i::unit( (xyz+1)%3 );
        Vector3i shiftZ = -m*Vector3i::unit( (xyz+2)%3 );
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantE(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 maskD(voxelGrid.inversePermittivity()
            .support(xyz, xyz, numerOrder, denomOrder));
        maskD *= voxelGrid.fieldIndices().je(xyz);
        maskD *= calcCells;
        //LOG << "E" << char('x'+xyz) << " with source, denom " << denomOrder
        //    << " has " << maskD.numRuns() << " runs.\n";
        
        SupportRegion3 maskBjPlus, maskBjMinus, maskBkPlus, maskBkMinus;
        maskBjPlus = maskD;
        translate(maskD, maskBjMinus, shiftZ[0], shiftZ[1], shiftZ[2]);
        maskBkPlus = maskD;
        translate(maskD, maskBkMinus, shiftY[0], shiftY[1], shiftY[2]);
        
        IndexArray3 indicesBjPlus, indicesBjMinus,
            indicesBkPlus, indicesBkMinus, indicesAux(maskD);
        restriction(voxelGrid.fieldIndices().b((xyz+1)%3), maskBjMinus,
            indicesBjMinus);
        restriction(voxelGrid.fieldIndices().b((xyz+1)%3), maskBjPlus,
            indicesBjPlus);
        restriction(voxelGrid.fieldIndices().b((xyz+2)%3), maskBkMinus,
            indicesBkMinus);
        restriction(voxelGrid.fieldIndices().b((xyz+2)%3), maskBkPlus,
            indicesBkPlus);
        
        IndexArray3 permittivityIndices(
            voxelGrid.inversePermittivity()
                .filter(xyz, xyz, numerOrder, denomOrder));
        IndexArray3 indicesE, indicesEps;
        restriction(voxelGrid.fieldIndices().e(xyz), maskD, indicesE);
        restriction(permittivityIndices, maskD, indicesEps);
        
        IndexArray3 indicesJE;
        restriction(voxelGrid.fieldIndices().je(xyz), maskD, indicesJE);
        
        IndexArray3 const* arrays[] = { &indicesE,
            &indicesBjMinus, &indicesBjPlus,
            &indicesBkMinus, &indicesBkPlus,
            &indicesAux, &indicesEps, &indicesJE };
        
        vector<RunlineAdjointEH_JM> runlines;
        CallbackAdjointEH_JM callback(runlines);
        RLE::merge(arrays, 8, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateE(new AdjointEH_JM(
                Field::e(xyz), Field::b((xyz+1)%3), Field::b((xyz+2)%3),
                Field::je(xyz),
                denomOrder,
                -1.0/gridDescription()->dxyz()[(xyz+1)%3],
                1.0/gridDescription()->dxyz()[(xyz+2)%3],
                -1.0,
                runlines,
                voxelGrid.inversePermittivity()
                    .filter(xyz, xyz, numerOrder, denomOrder).values(),
                Constants::eps0));
        }
    }
}

void AdjointGridOperations::
encodeB(const VoxelizedPartition & voxelGrid)
{    
    // Without current sources
    int denomOrder = -1; // matches anything
    
    for (int numerOrder = 0; numerOrder < 10; numerOrder++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantH(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 maskB(voxelGrid.inversePermeability()
            .support(xyz, xyz, numerOrder, denomOrder));
        maskB -= voxelGrid.fieldIndices().m(xyz);
        maskB *= calcCells;
        IndexArray3 permeabilityIndices(
            voxelGrid.inversePermeability()
            .filter(xyz, xyz, numerOrder, denomOrder));
        
        IndexArray3 indicesB, indicesH, indicesMu, indicesAux(maskB);
        restriction(voxelGrid.fieldIndices().b(xyz), maskB, indicesB);
        restriction(voxelGrid.fieldIndices().h(xyz), maskB, indicesH);
        restriction(permeabilityIndices, maskB, indicesMu);
        
        IndexArray3 const* arrays[] = { &indicesB,
            &indicesH, &indicesAux, &indicesMu };
        
        vector<RunlineAdjointDB> runlines;
        CallbackAdjointDB callback(runlines);
        RLE::merge(arrays, 4, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateH(new AdjointDB(
                Field::b(xyz), Field::h(xyz),
                numerOrder,
                runlines,
                voxelGrid.inversePermeability()
                    .filter(xyz, xyz, numerOrder, denomOrder).values(),
                gridDescription()->dt()));
        }
    }
    
    // With current sources
    for (int numerOrder = 0; numerOrder < 10; numerOrder++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantH(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 maskB(voxelGrid.inversePermeability()
            .support(xyz, xyz, numerOrder, denomOrder));
        maskB *= voxelGrid.fieldIndices().m(xyz);
        maskB *= calcCells;
        IndexArray3 permeabilityIndices(
            voxelGrid.inversePermeability()
                .filter(xyz, xyz, numerOrder, denomOrder));

        IndexArray3 indicesB, indicesH, indicesMu, indicesAux(maskB);
        restriction(voxelGrid.fieldIndices().b(xyz), maskB, indicesB);
        restriction(voxelGrid.fieldIndices().h(xyz), maskB, indicesH);
        restriction(permeabilityIndices, maskB, indicesMu);
        
        IndexArray3 indicesM;
        restriction(voxelGrid.fieldIndices().m(xyz), maskB, indicesM);
        
        IndexArray3 const* arrays[] = { &indicesB,
            &indicesH, &indicesAux, &indicesMu, &indicesM };
        
        vector<RunlineAdjointDB_JM> runlines;
        CallbackAdjointDB_JM callback(runlines);
        RLE::merge(arrays, 5, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateH(new AdjointDB_JM(
                Field::b(xyz), Field::h(xyz), Field::m(xyz),
                numerOrder,
                -gridDescription()->dt(),
                runlines,
                voxelGrid.inversePermeability()
                    .filter(xyz, xyz, numerOrder, denomOrder).values(),
                gridDescription()->dt()));
        }
    }
}


void AdjointGridOperations::
encodeH(const VoxelizedPartition & voxelGrid)
{
    Matrix3i m = voxelGrid.projectionMatrix();
    
//    LOG << "Encoding diagonal H, crudely\n";
    int numerOrder = -1;
    
    for (int denomOrder = 0; denomOrder < 10; denomOrder++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Vector3i shiftY = m*Vector3i::unit( (xyz+1)%3 );
        Vector3i shiftZ = m*Vector3i::unit( (xyz+2)%3 );
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantH(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 maskB(voxelGrid.inversePermeability()
            .support(xyz, xyz, numerOrder, denomOrder));
        maskB -= voxelGrid.fieldIndices().mh(xyz);
        maskB *= calcCells;
        
        SupportRegion3 maskDjPlus, maskDjMinus, maskDkPlus, maskDkMinus;
        maskDjMinus = maskB;
        translate(maskB, maskDjPlus, shiftZ[0], shiftZ[1], shiftZ[2]);
        maskDkMinus = maskB;
        translate(maskB, maskDkPlus, shiftY[0], shiftY[1], shiftY[2]);
        
        IndexArray3 indicesDjPlus, indicesDjMinus,
            indicesDkPlus, indicesDkMinus, indicesAux(maskB);
        restriction(voxelGrid.fieldIndices().d((xyz+1)%3), maskDjMinus,
            indicesDjMinus);
        restriction(voxelGrid.fieldIndices().d((xyz+1)%3), maskDjPlus,
            indicesDjPlus);
        restriction(voxelGrid.fieldIndices().d((xyz+2)%3), maskDkMinus,
            indicesDkMinus);
        restriction(voxelGrid.fieldIndices().d((xyz+2)%3), maskDkPlus,
            indicesDkPlus);
        
        IndexArray3 permeabilityIndices(
            voxelGrid.inversePermeability()
                .filter(xyz, xyz, numerOrder, denomOrder));
        IndexArray3 indicesH, indicesMu;
        restriction(voxelGrid.fieldIndices().h(xyz), maskB, indicesH);
        restriction(permeabilityIndices, maskB, indicesMu);
        
        IndexArray3 const* arrays[] = { &indicesH,
            &indicesDjMinus, &indicesDjPlus,
            &indicesDkMinus, &indicesDkPlus,
            &indicesAux, &indicesMu };
        
        vector<RunlineAdjointEH> runlines;
        CallbackAdjointEH callback(runlines);
        RLE::merge(arrays, 7, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateH(new AdjointEH(
                Field::h(xyz), Field::d((xyz+1)%3), Field::d((xyz+2)%3),
                denomOrder,
                1.0/gridDescription()->dxyz()[(xyz+1)%3],
                -1.0/gridDescription()->dxyz()[(xyz+2)%3],
                runlines,
                voxelGrid.inversePermeability()
                    .filter(xyz, xyz, numerOrder, denomOrder).values(),
                Constants::mu0));
        }
    }
    
//    LOG << "Diagonal H with sources...\n";
    
    for (int denomOrder = 0; denomOrder < 10; denomOrder++)
    for (int xyz = 0; xyz < 3; xyz++)
    {
        Vector3i shiftY = m*Vector3i::unit( (xyz+1)%3 );
        Vector3i shiftZ = m*Vector3i::unit( (xyz+2)%3 );
        
        SupportRegion3 calcCells = support(
            halfToYee(voxelGrid.extents().calcHalfCells(), octantH(xyz)),
            voxelGrid.extents().dimensions());
        SupportRegion3 maskB(voxelGrid.inversePermeability()
            .support(xyz, xyz, numerOrder, denomOrder));
        maskB *= voxelGrid.fieldIndices().mh(xyz);
        maskB *= calcCells;
        
        SupportRegion3 maskDjPlus, maskDjMinus, maskDkPlus, maskDkMinus;
        maskDjMinus = maskB;
        translate(maskB, maskDjPlus, shiftZ[0], shiftZ[1], shiftZ[2]);
        maskDkMinus = maskB;
        translate(maskB, maskDkPlus, shiftY[0], shiftY[1], shiftY[2]);
        
        IndexArray3 indicesDjPlus, indicesDjMinus,
            indicesDkPlus, indicesDkMinus, indicesAux(maskB);
        restriction(voxelGrid.fieldIndices().d((xyz+1)%3), maskDjMinus,
            indicesDjMinus);
        restriction(voxelGrid.fieldIndices().d((xyz+1)%3), maskDjPlus,
            indicesDjPlus);
        restriction(voxelGrid.fieldIndices().d((xyz+2)%3), maskDkMinus,
            indicesDkMinus);
        restriction(voxelGrid.fieldIndices().d((xyz+2)%3), maskDkPlus,
            indicesDkPlus);
        
        IndexArray3 permeabilityIndices(
            voxelGrid.inversePermeability()
                .filter(xyz, xyz, numerOrder, denomOrder));
        IndexArray3 indicesH, indicesMu;
        restriction(voxelGrid.fieldIndices().h(xyz), maskB, indicesH);
        restriction(permeabilityIndices, maskB, indicesMu);
        
        IndexArray3 indicesMH;
        restriction(voxelGrid.fieldIndices().mh(xyz), maskB, indicesMH);
        
        IndexArray3 const* arrays[] = { &indicesH,
            &indicesDjMinus, &indicesDjPlus,
            &indicesDkMinus, &indicesDkPlus,
            &indicesAux, &indicesMu, &indicesMH };
        
        vector<RunlineAdjointEH_JM> runlines;
        CallbackAdjointEH_JM callback(runlines);
        RLE::merge(arrays, 8, callback);
        
        if (runlines.size() > 0)
        {
            pushUpdateH(new AdjointEH_JM(
                Field::h(xyz), Field::d((xyz+1)%3), Field::d((xyz+2)%3),
                Field::mh(xyz),
                denomOrder,
                1.0/gridDescription()->dxyz()[(xyz+1)%3],
                -1.0/gridDescription()->dxyz()[(xyz+2)%3],
                -1.0,
                runlines,
                voxelGrid.inversePermeability()
                    .filter(xyz, xyz, numerOrder, denomOrder).values(),
                Constants::mu0));
        }
    }
}

void AdjointGridOperations::
encodePMLJ(const VoxelizedPartition & voxelGrid)
{
    Matrix3i m = voxelGrid.projectionMatrix();
    
    LOG << "Encoding adjoint PML JE updates\n";
    for (int xyz = 0; xyz < 3; xyz++)
    {   
        for (int mm = 1; mm <= 2; mm++)
        {
            Precision::Float leadingSign = (mm == 1 ? 1.0f : -1.0f);
            int absorbXYZ = (xyz+mm)%3;
            int bXYZ = (xyz+3-mm)%3;
            assert(absorbXYZ != xyz);
            assert(bXYZ != xyz);
            assert(absorbXYZ != bXYZ);
            
            Vector3i shift = -m*Vector3i::unit(absorbXYZ);
            
            SupportRegion3 maskPML(voxelGrid.pmlIndices().e(xyz, absorbXYZ));
            SupportRegion3 maskBPlus = maskPML;
            SupportRegion3 maskBMinus;
            translate(maskPML, maskBMinus, shift[0], shift[1], shift[2]);
            
            IndexArray3 indicesPML(maskPML); // could get from pmlIndices()...
            IndexArray3 indicesJE, indicesBPlus, indicesBMinus, indicesConsts;
            restriction(voxelGrid.fieldIndices().je(xyz), maskPML, indicesJE);
            restriction(voxelGrid.fieldIndices().b(bXYZ), maskBMinus,
                indicesBMinus);
            restriction(voxelGrid.fieldIndices().b(bXYZ), maskBPlus,
                indicesBPlus);
            restriction(
                IndexArray3(voxelGrid.pmlE(xyz, absorbXYZ)),
                maskPML, indicesConsts);
            
            IndexArray3 const* arrays[] = { &indicesJE,
                &indicesBPlus, &indicesBMinus,
                &indicesPML, &indicesConsts };
            
            vector<RunlineJMPML> runlines;
            CallbackJMPML callback(runlines);
            RLE::merge(arrays, 5, callback);
            
            if (runlines.size() > 0)
            {
                pushUpdateE(new UpdateJMPML(
                    Field::je(xyz), Field::b(bXYZ),
                    leadingSign,
                    runlines,
                    voxelGrid.pmlE(xyz, absorbXYZ).values(),
                    gridDescription()->dt(),
                    gridDescription()->dxyz()[absorbXYZ]));
            }
        }
    }
}

void AdjointGridOperations::
encodePMLM(const VoxelizedPartition & voxelGrid)
{
    Matrix3i m = voxelGrid.projectionMatrix();
    
    LOG << "Encoding adjoint PML MH updates\n";
    for (int xyz = 0; xyz < 3; xyz++)
    {   
        for (int mm = 1; mm <= 2; mm++)
        {
            Precision::Float leadingSign = (mm == 1 ? -1.0f : 1.0f);
            int absorbXYZ = (xyz+mm)%3;
            int dXYZ = (xyz+3-mm)%3;
            assert(absorbXYZ != xyz);
            assert(dXYZ != xyz);
            assert(absorbXYZ != dXYZ);
            
            Vector3i shift = m*Vector3i::unit(absorbXYZ);
            
            SupportRegion3 maskPML(voxelGrid.pmlIndices().h(xyz, absorbXYZ));
            SupportRegion3 maskDPlus;
            translate(maskPML, maskDPlus, shift[0], shift[1], shift[2]);
            SupportRegion3 maskDMinus = maskPML;
            
            IndexArray3 indicesPML(maskPML); // could get from pmlIndices()...
            IndexArray3 indicesMH, indicesDPlus, indicesDMinus, indicesConsts;
            restriction(voxelGrid.fieldIndices().mh(xyz), maskPML, indicesMH);
            restriction(voxelGrid.fieldIndices().d(dXYZ), maskDMinus,
                indicesDMinus);
            restriction(voxelGrid.fieldIndices().d(dXYZ), maskDPlus,
                indicesDPlus);
            restriction(
                IndexArray3(voxelGrid.pmlH(xyz, absorbXYZ)),
                maskPML, indicesConsts);

            IndexArray3 const* arrays[] = { &indicesMH,
                &indicesDPlus, &indicesDMinus,
                &indicesPML, &indicesConsts };
            
            vector<RunlineJMPML> runlines;
            CallbackJMPML callback(runlines);
            RLE::merge(arrays, 5, callback);
            
            if (runlines.size() > 0)
            {
                pushUpdateH(new UpdateJMPML(
                    Field::mh(xyz), Field::d(dXYZ),
                    leadingSign,
                    runlines,
                    voxelGrid.pmlH(xyz, absorbXYZ).values(),
                    gridDescription()->dt(),
                    gridDescription()->dxyz()[absorbXYZ]));
            }
        }
    }
}
