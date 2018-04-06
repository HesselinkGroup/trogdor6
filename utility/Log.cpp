/*
 *  Log.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 12/10/07.
 *  Copyright 2007 Stanford University. All rights reserved.
 *  This file is covered by the MIT license.  See LICENSE.txt.
 *
 */

#include "Log.h"

namespace Log
{

std::ofstream TrogLog::sLogfile("log.out.txt");
StreamTee TrogLog::sTee(std::cout, TrogLog::sLogfile);

std::string TrogLog::
stripArgs(const std::string & str)
{
    int c2 = str.find_first_of("(");  // find the start of the arguments block
    std::string noArgs = str.substr(0,c2);
    int c1 = noArgs.find_last_of(" ");  // find the beginning of the invocation
    if (c1 == std::string::npos)
        c1 = -1;
    return str.substr(c1+1, c2-c1) + ")";
}


}; // namespace Log
