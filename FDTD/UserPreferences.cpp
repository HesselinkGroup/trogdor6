/*
 *  UserPreferences.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/28/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "UserPreferences.h"

UserPreferences UserPreferences::sInstance;

UserPreferences::
UserPreferences()
{
}

std::string UserPreferences::
myValue(const std::string & key) const
{
    if (mDictionary.count(key) != 0)
        return mDictionary.find(key)->second;
    return "";
}

void UserPreferences::
mySet(const std::string & key, const std::string & val)
{
    mDictionary[key] = val;
}

bool UserPreferences::
myDefines(const std::string & key) const
{
    return (mDictionary.count(key) != 0);
}
