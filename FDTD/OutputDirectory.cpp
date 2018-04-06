/*
 *  OutputDirectory.cpp
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 5/19/12.
 *  Copyright 2012 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "OutputDirectory.h"

#include "UserPreferences.h"

namespace OutputDirectory
{

std::string path(std::string fileName)
{
    std::string str;
    
    if (UserPreferences::defines("outputDirectory"))
    {
        str = UserPreferences::valueAs<std::string>("outputDirectory")
            + "/" + fileName;
    }
    else
        str = fileName;
    
    return str;
}

};
