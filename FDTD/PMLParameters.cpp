/*
 *  PMLParameters.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 4/1/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "PMLParameters.h"

#include <cstdlib>
#include <cmath>

using namespace std;

CalculatePMLParameters::
CalculatePMLParameters(double dx, double l, double eps0, double mu0) :
    m_dx(dx),
    m_eps0(eps0),
    m_mu0(mu0),
    m_c(1.0 /sqrt(eps0*mu0)),
    m_L(l)
{
    mSigmaEquation = "(d^3)*0.8*4/(((mu0/eps0)^0.5)*dx)";
    mAlphaEquation = "(1-d)*c*eps0/L";
    mKappaEquation = "1 + (5-1)*(d^3)";
}

PMLParameters CalculatePMLParameters::
operator()(double depth) const
{
    if (depth > 0.0 && depth <= 1.0)
    {
        return PMLParameters(sigma(depth), kappa(depth), alpha(depth));
    }
    return PMLParameters();
}

double CalculatePMLParameters::
kappa(double depth) const
{
    double k = 1.0 + (5.0-1.0)*depth*depth*depth;
    return k;
}

double CalculatePMLParameters::
alpha(double depth) const
{
    double a = (1.0-depth)*m_c*m_eps0/m_L;
    return a;
}

double CalculatePMLParameters::
sigma(double depth) const
{
    double s = (depth*depth*depth)*0.8*4.0/(sqrt(m_mu0/m_eps0)*m_dx);
    return s;
}

