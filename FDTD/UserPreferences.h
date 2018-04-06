/*
 *  UserPreferences.h
 *  Trogdor6
 *
 *  Created by Paul Hansen on 3/28/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _USERPREFERENCES_
#define _USERPREFERENCES_

#include <string>
#include <map>
#include <sstream>
#include <stdexcept>

class UserPreferences
{
public:
    static std::string value(const std::string & key)
    {
        return instance().myValue(key);
    }
    static void set(const std::string & key, const std::string & value = "1")
    {
        instance().mySet(key, value);
    }
    static bool defines(const std::string & key)
    {
        return instance().myDefines(key);
    }
    template<class T>
    static T valueAs(const std::string & key)
    {
        return instance().myValueAs<T>(key);
    }
    template<class T>
    static void set(const std::string & key, const T & value)
    {
        instance().mySet(key, value);
    }
    
    static const std::map<std::string, std::string> & dictionary()
    {
        return instance().myDictionary();
    }
    
private:
    static UserPreferences & instance() { return sInstance; }
    std::string myValue(const std::string & key) const;
    void mySet(const std::string & key, const std::string & value);
    bool myDefines(const std::string & key) const;
    
    const std::map<std::string, std::string> & myDictionary() const
    {
        return mDictionary;
    }
    
    template<class T>
    T myValueAs(const std::string & key)
    {
        if (!defines(key))
            throw std::logic_error("Undefined key");
        T returnValue;
        std::istringstream str(value(key));
        str >> returnValue;
        if (str.fail())
            throw std::logic_error("Type conversion failed");
        return returnValue;
    }
    
    template<class T>
    void mySet(const std::string & key, const T & value)
    {
        std::ostringstream str;
        str << value;
        if (str.fail())
            throw std::logic_error("String conversion failed");
        mySet(key, str.str());
    }
    
private:
    UserPreferences();
    static UserPreferences sInstance;
    std::map<std::string, std::string> mDictionary;
};


#endif
