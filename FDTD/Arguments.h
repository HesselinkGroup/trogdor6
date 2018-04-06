/*
 *  Arguments.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 8/11/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _ARGUMENTS_
#define _ARGUMENTS_

#include <boost/program_options.hpp>

boost::program_options::variables_map
handleArguments(int argc, char* const argv[]);

#endif
