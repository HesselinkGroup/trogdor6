/*
 *  PMLParameters.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/5/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _PMLPARAMETERS_
#define _PMLPARAMETERS_

#include <string>
#include "Precision.h"

class PMLParameters
{
public:
    PMLParameters() :
        mSigma(0.0), mKappa(1.0), mAlpha(0.0)
    {
    }
    
    PMLParameters(Precision::Float sigma, Precision::Float kappa, Precision::Float alpha) :
        mSigma(sigma), mKappa(kappa), mAlpha(alpha)
    {
    }
    
    Precision::Float sigma() const { return mSigma; }
    Precision::Float kappa() const { return mKappa; }
    Precision::Float alpha() const { return mAlpha; }
    
    bool operator==(const PMLParameters & rhs) const
    {
        return (mSigma == rhs.mSigma && mKappa == rhs.mKappa &&
            mAlpha == rhs.mAlpha);
    }
    
private:
    Precision::Float mSigma;
    Precision::Float mKappa;
    Precision::Float mAlpha;
};

template <class STREAM>
STREAM & operator<<(STREAM & str, const PMLParameters & p)
{
    str << "(sigma = " << p.sigma() << " kappa = " << p.kappa()
        << " alpha = " << p.alpha() << ")";
    return str;
}

class CalculatePMLParameters
{
public:
    CalculatePMLParameters(double yeeCellSize, double gridDiagonalLength, 
        double eps0, double mu0);
    
    PMLParameters operator()(double depth) const;
    double kappa(double depth) const;
    double alpha(double depth) const;
    double sigma(double depth) const;
    
    void kappaEquation(std::string formula) { mKappaEquation = formula; }
    void alphaEquation(std::string formula) { mAlphaEquation = formula; }
    void sigmaEquation(std::string formula) { mSigmaEquation = formula; }
    
    std::string kappaEquation() const { return mKappaEquation; }
    std::string alphaEquation() const { return mAlphaEquation; }
    std::string sigmaEquation() const { return mSigmaEquation; }
private:
    double eval(const std::string & formula, double depth) const;
    
    double m_dx;
    double m_eps0;
    double m_mu0;
    double m_c;
    double m_L; // diagonal dimension of the whole grid
    std::string mSigmaEquation;
    std::string mAlphaEquation;
    std::string mKappaEquation;
};




#endif
