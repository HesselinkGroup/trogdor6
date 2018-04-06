/*
 *  TensorGridConstants.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 4/14/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _TENSORGRIDCONSTANTS_
#define _TENSORGRIDCONSTANTS_

#include "PrecisionRationalFunction.h"
#include "OrientedVoxels.h"

#include <iostream>

class TensorGridConstants
{
public:
    enum AveragingType
    {
        kHarmonicAveraging,
        kArithmeticAveraging,
        kMixedAveraging
    };
    
    // -------------- NEW STUFF
    
    TensorGridConstants(int i, int j,
        RLE::DynamicRLE3<Precision::RationalFunction> & invConstits,
        Rect3i voxelBounds,
        Precision::RationalFunction backgroundMaterial = 0.0);
    
    void calculateConstitutives(
        const std::vector<Precision::RationalFunction> & materials,
        const std::vector<RLE::DynamicRLE3<double> > & fillFactors,
        const RLE::DynamicRLE3<double> & orientations);
    
    // this needs the controlVertices array for two reasons:
    // 1. to get the # of CVs
    // 2. to get the free directions
    void writeSensitivities(
        const std::vector<Precision::RationalFunction> & materials,
        const std::vector<SimpleMesh::ControlVertex> & controlVertices,
        const std::vector<RLE::DynamicRLE3<double> > & fillFactors,
        const RLE::DynamicRLE3<double> & orientations,
        const std::vector<std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > > &
            fillFactorJacobians,
        const std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > & 
            orientationJacobians,
        const RLE::SupportRegion3 & savedFields,
        std::ostream & file = std::cout);
    
    // T5 mode
    RLE::DynamicRLE3<Precision::RationalFunction> winnerTakesAll(
        const std::vector<Precision::RationalFunction> & materials,
        const std::vector<RLE::DynamicRLE3<double> > & fillFactors);
    
    // tested simply
    RLE::DynamicRLE3<Precision::RationalFunction> diagonalElement(
        const std::vector<Precision::RationalFunction> & materials,
        const std::vector<RLE::DynamicRLE3<double> > & fillFactors,
        const RLE::DynamicRLE3<double> & orientation_ij);
    
    // tested simply
    RLE::DynamicRLE3<Precision::RationalFunction> offDiagonalElement(
        const std::vector<Precision::RationalFunction> & materials,
        const std::vector<RLE::DynamicRLE3<double> > & fillFactors,
        const RLE::DynamicRLE3<double> & orientation_ij);
    
    // Calculate the sensitivity of the material constants by a centered
    // difference.
    // tested for simple dielectric only, separately varying fill factor and
    // orientation.
    Pointer<RLE::DynamicRLE3<Precision::RationalFunction> >
    diagonalSensitivity(
        const std::vector<Precision::RationalFunction> & materials,
        const std::vector<RLE::DynamicRLE3<double> > & fillFactors,
        const RLE::DynamicRLE3<double> & orientations,
        const std::vector<RLE::DynamicRLE3<double> > & fillFactorSensitivities,
        const RLE::DynamicRLE3<double> & orientationSensitivities);
        
    Pointer<RLE::DynamicRLE3<Precision::RationalFunction> >
    offDiagonalSensitivity(
        const std::vector<Precision::RationalFunction> & materials,
        const std::vector<RLE::DynamicRLE3<double> > & fillFactors,
        const RLE::DynamicRLE3<double> & orientations,
        const std::vector<RLE::DynamicRLE3<double> > & fillFactorSensitivities,
        const RLE::DynamicRLE3<double> & orientationSensitivities);
    
    static void selectFillFactorSensitivity(
        const std::vector<std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > > &
            fillFactorJacobians,
        int vertexNumber, int xyz,
        std::vector<RLE::DynamicRLE3<double> > & out);
    
    static void selectOrientationSensitivity(
        const std::map<unsigned int, RLE::DynamicRLE3<Vector3d> > &
            orientationJacobians,
        int vertexNumber, int xyz,
        RLE::DynamicRLE3<double> & out);
    
    void calculateHarmonicMean(
        const std::vector<RLE::DynamicRLE3<double> > & yeeFillFactors,
        std::vector<Precision::RationalFunction> inverseConstitutives,
        RLE::DynamicRLE3<Precision::RationalFunction> & outHarmonicMean);
        
    void calculateArithmeticMean(
        const std::vector<RLE::DynamicRLE3<double> > & yeeFillFactors,
        std::vector<Precision::RationalFunction> inverseConstitutives,
        RLE::DynamicRLE3<Precision::RationalFunction> & outArithmeticMean);
    
    struct Differential
    {
        Differential(Precision::Float delta) : mFactor(1.0/delta) {}
        
        Precision::RationalFunction operator()(const Precision::RationalFunction & r1,
            const Precision::RationalFunction & r2) const
        {
            Precision::RationalFunction r = Precision::RationalFunction::irregular(
                mFactor*(r1.numerator()-r2.numerator()),
                mFactor*(r1.denominator()-r2.denominator()) );
            return r;
        }
        
        Precision::Float mFactor;
    };
    
    const Precision::RationalFunction & backgroundMaterial() const
        { return mBackground; }
    
    const Rect3i & voxels() const { return mVoxels; }
    
private:
    void writeSensitivity(
        Pointer<RLE::DynamicRLE3<Precision::RationalFunction> > sensitivity,
        const RLE::SupportRegion3 & savedFields,
        std::ostream & file) const;
    
    RLE::DynamicRLE3<Precision::RationalFunction> & mInverseConstitutives;
    
    Precision::RationalFunction mBackground;
    Rect3i mVoxels;
    int m_i, m_j;
};



#endif
