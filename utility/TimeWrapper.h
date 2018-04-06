/*
 *  TimeWrapper.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 7/12/05.
 *  Copyright 2005 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _TIMEWRAPPER_
#define _TIMEWRAPPER_

#include <chrono>
#include <string>

namespace TimeWrapper
{

typedef std::chrono::time_point<std::chrono::steady_clock> TimePoint;

inline TimePoint now()
{
    return std::chrono::steady_clock::now();
}

double elapsedMicroseconds(TimePoint t0, TimePoint t1);
double elapsedSeconds(TimePoint t0, TimePoint t1);

std::string timestampNow();

}; // namespace TimeWrapper

#endif
