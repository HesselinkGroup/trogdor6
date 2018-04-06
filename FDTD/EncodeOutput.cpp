/*
 *  EncodeOutput.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/23/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "EncodeOutput.h"

#include "SimulationDescription.h"
#include "VoxelizedPartition.h"
#include "operations/OutputEHDB.h"

#include "rle/SupportRegion3.h"
#include "rle/MergeRunlines.h"

#include "Support.h"

#include <vector>

using namespace std;
using namespace RLE;
using namespace YeeUtilities;

struct CallbackOutput
{
    CallbackOutput(vector<RunlineOutput> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        RunlineOutput runline(indices[0], length);
        mRunlines.push_back(runline);
    }
    
    vector<RunlineOutput> & mRunlines;
};
typedef CallbackOutput CallbackSource;

struct CallbackSparseOutput
{
    CallbackSparseOutput(vector<RunlineOutput> & runlines) : mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        
        RunlineOutput runline(indices[0], length);
        mRunlines.push_back(runline);
    }
    
    vector<RunlineOutput> & mRunlines;
};

struct CallbackOutputAverageTwo
{
    CallbackOutputAverageTwo(vector<RunlineAverageTwo> & runlines) :
        mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        RunlineAverageTwo runline(indices[0], indices[1], length);
        mRunlines.push_back(runline);
    }
    
    vector<RunlineAverageTwo> & mRunlines;
};

struct CallbackTensorialOutput
{
    CallbackTensorialOutput(
        vector<RunlineTensorialOutput> & runlines) :
        mRunlines(runlines) {}
    
    void operator()(const std::int64_t indices[],
        const std::uint64_t strides[],
        std::uint64_t length,
        std::uint64_t numStreams)
    {
        RunlineTensorialOutput runline(indices[0], indices[1], length);
        mRunlines.push_back(runline);
    }
    
    vector<RunlineTensorialOutput> & mRunlines;
};

EncodeOutput::
EncodeOutput(const GridDescription & gridDesc,
    const VoxelizedPartition & vp) :
    mGridDesc(gridDesc),
    mVoxelGrid(vp)
{
}


Pointer<OutputEHDB> EncodeOutput::
encode(const OutputDescription & outputDesc) const
{    
    Pointer<OutputEHDB> output(new OutputEHDB(outputDesc,
        gridDescription().dxyz(), gridDescription().dt(),
        gridDescription().origin()));
    
    SupportRegion3 supp, suppShifted;
    
    for (int ff = 0; ff < outputDesc.fields().size(); ff++)
    {
        const OutputField & field(outputDesc.fields().at(ff));
        
        if (field.field() == kEE || field.field() == kHH)
        {
            if (field.ii() == field.jj())
                encodeTensorialDiagonalOutput(*output, outputDesc, field);
            else
                encodeTensorialOffDiagonalOutput(*output, outputDesc, field);
        }
        else if (field.octant() == 0 || field.octant() == 7)
            encodeInterpolatedOutput(*output, outputDesc, field);
        else if (field.field() == kJ || field.field() == kM)
            encodeSparseOutput(*output, outputDesc, field);
        else
            encodeStandardOutput(*output, outputDesc, field);
    }
    
    return output;
}


void EncodeOutput::
encodeStandardOutput(OutputEHDB & output,
    const OutputDescription & outputDesc,
    const OutputField & field) const
{
    assert(field.ii() == field.jj());
    int xyz = field.ii();
    
    const IndexArray3 & fieldIndices = 
        voxelGrid().fieldIndices().field(field.field(), xyz);
    
    vector<RunlineOutput> runlines;
    for (int rr = 0; rr < outputDesc.regions().size(); rr++)
    {
        SupportRegion3 supp = support(outputDesc.regions().at(rr).yeeCells(),
            Vector3b(voxelGrid().extents().nodeYeeCells().num()));
        
        IndexArray3 outputIndices;
        
        restriction(fieldIndices, supp, outputIndices);
        
        IndexArray3 const* inds[] = { &outputIndices };
        
        CallbackOutput callback(runlines);
        RLE::mergeOverlapping(inds, 1, callback);
    }
    
    output.appendField(field, runlines);
}

void EncodeOutput::
encodeSparseOutput(OutputEHDB & output,
    const OutputDescription & outputDesc,
    const OutputField & field) const
{
    assert(field.ii() == field.jj());
    int xyz = field.ii();
    
    const IndexArray3 & fieldIndices = 
        voxelGrid().fieldIndices().field(field.field(), xyz);
    
    vector<RunlineOutput> runlines;
    for (int rr = 0; rr < outputDesc.regions().size(); rr++)
    {
        SupportRegion3 supp = support(outputDesc.regions().at(rr).yeeCells(),
            Vector3b(voxelGrid().extents().nodeYeeCells().num()));
        SupportRegion3 emptyCells = supp - fieldIndices;
        
        IndexArray3 outputIndices, emptyIndices(emptyCells);
        restriction(fieldIndices, supp, outputIndices);
        
        IndexArray3 const* inds[] = { &outputIndices, &emptyIndices };
        
        // this DOES work, using the standard output callback.
        CallbackOutput callback(runlines);
        RLE::mergeOverlapping(inds, 2, callback);
    }
    
    output.appendField(field, runlines);
}

void EncodeOutput::
encodeInterpolatedOutput(OutputEHDB & output,
    const OutputDescription & outputDesc,
    const OutputField & field) const
{
    assert(field.ii() == field.jj());
    int xyz = field.ii();
    
    const IndexArray3 & fieldIndices = 
        voxelGrid().fieldIndices().field(field.field(), xyz);
    vector<RunlineAverageTwo> runlines;
    
    for (int rr = 0; rr < outputDesc.regions().size(); rr++)
    {
        SupportRegion3 supp = support(outputDesc.regions().at(rr).yeeCells(),
            Vector3b(voxelGrid().extents().nodeYeeCells().num()));
        SupportRegion3 suppShifted;
        if (field.field() == kE || field.field() == kD || field.field() == kJ)
        {
            assert(field.octant() == 0);
            translate(supp, suppShifted, (-Vector3<std::int64_t>::unit(xyz)).asArray());
        }
        else if (field.field() == kH || field.field() == kB || field.field() == kM)
        {
            assert(field.octant() == 7);
            translate(supp, suppShifted, Vector3<std::int64_t>::unit(xyz).asArray());
        }
        else
            throw(std::logic_error("Bad field"));
        
        IndexArray3 outputIndices, outputIndicesShifted;
        restriction(fieldIndices, supp, outputIndices);
        restriction(fieldIndices, suppShifted, outputIndicesShifted);
        
        IndexArray3 const* inds[] = { &outputIndices,
            &outputIndicesShifted };
        CallbackOutputAverageTwo callback(runlines);
        RLE::mergeOverlapping(inds, 2, callback);
    }
    
    output.appendField(field, runlines);
}


void EncodeOutput::
encodeTensorialDiagonalOutput(OutputEHDB & output,
    const OutputDescription & outputDesc,
    const OutputField & field) const
{
    assert(field.ii() == field.jj());
    assert(field.field() == kEE || field.field() == kHH);
    int xyz = field.ii();
    
    // These are the fields we want to output, e.g. Exx.
    const IndexArray3 & tensorialIndices =
        voxelGrid().fieldIndices().field(field.field(), xyz, xyz);
    vector<RunlineTensorialOutput> runlines;
    
    for (int rr = 0; rr < outputDesc.regions().size(); rr++)
    {
        SupportRegion3 supp = support(outputDesc.regions().at(rr).yeeCells(),
            Vector3b(voxelGrid().extents().nodeYeeCells().num()));
        
        // These are the "backup fields."  For diagonal permittivities, the Exx
        // field is just ordinary Ex, since there are no Exy and Exz components.
        const IndexArray3 & vectorialIndices = (field.field() == kEE) ? 
            voxelGrid().fieldIndices().e(xyz) :
            voxelGrid().fieldIndices().h(xyz);
        
        IndexArray3 mainIndices, backupIndices;
        restriction(tensorialIndices, supp, mainIndices);
        restriction(vectorialIndices, supp, backupIndices);
        
        IndexArray3 const* inds[] = { &mainIndices,
            &backupIndices };
        CallbackTensorialOutput callback(runlines);
        RLE::mergeOverlapping(inds, 2, callback);
        
    }
    
    output.appendTensorialDiagonalField(field, runlines);
}


void EncodeOutput::
encodeTensorialOffDiagonalOutput(OutputEHDB & output,
    const OutputDescription & outputDesc,
    const OutputField & field) const
{
    assert(field.ii() != field.jj());
    assert(field.field() == kEE || field.field() == kHH);
    
    const IndexArray3 & tensorialIndices =
        voxelGrid().fieldIndices().field(field.field(), field.ii(), field.jj());
    vector<RunlineTensorialOutput> runlines;
    
    for (int rr = 0; rr < outputDesc.regions().size(); rr++)
    {    
        SupportRegion3 supp = support(outputDesc.regions().at(rr).yeeCells(),
            Vector3b(voxelGrid().extents().nodeYeeCells().num()));
        
        IndexArray3 outputIndices, emptyIndices(supp);
        restriction(tensorialIndices, supp, outputIndices);
        
        IndexArray3 const* inds[] = { &outputIndices,
            &emptyIndices };
        CallbackTensorialOutput callback(runlines);
        RLE::mergeOverlapping(inds, 2, callback);
    }
    
    output.appendTensorialOffDiagonalField(field, runlines);
}




