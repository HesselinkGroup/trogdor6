/*
 *  TimeWrapper.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 7/12/05.
 *  Copyright 2005 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 *
 */

#include "TimeWrapper.h"

#include <chrono>
#include <sstream>
#include <iomanip>

namespace TimeWrapper
{

double elapsedMicroseconds(TimePoint t0, TimePoint t1)
{
    long elapsedNanoseconds = (t1-t0)/std::chrono::nanoseconds(1);
    return elapsedNanoseconds * 1e-3;
}

double elapsedSeconds(TimePoint t0, TimePoint t1)
{
    return elapsedMicroseconds(t0, t1) * 1e-6;
}

std::string timestampNow()
{
    auto timePointNow = std::chrono::system_clock::now();
    std::time_t now = std::chrono::system_clock::to_time_t(timePointNow);
    
    std::ostringstream str;
    str << std::put_time(std::localtime(&now), "%F %T");
    return str.str();
}

}; // namespace TimeWrapper

