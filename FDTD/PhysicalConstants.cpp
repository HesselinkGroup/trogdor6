/*
 *  PhysicalConstants.cpp
 *  TROGDOR
 *
 *  Created by Paul Hansen on 7/14/06.
 *  Copyright 2007 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "PhysicalConstants.h"

#include <cmath>

//#define SI_UNITS

#ifdef SI_UNITS
const double Constants::c = 2.99792458e8;
const double Constants::eps0 = 8.854187817e-12;
const double Constants::mu0 = 4*M_PI*1e-7;
const double Constants::eta0 = 376.730313461;
#else
const double Constants::c = 1.0;
const double Constants::eps0 = 1.0;
const double Constants::mu0 = 1.0;
const double Constants::eta0 = 1.0;
#endif

