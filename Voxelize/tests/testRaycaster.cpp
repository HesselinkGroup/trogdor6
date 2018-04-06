/*
 *  testRaycaster.cpp
 *
 *  Created by Paul C Hansen on 7/5/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testRaycaster

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Raycaster.h"
#include "SimpleMesh/SimpleMesh.h"

using namespace std;

//static const int X = 0, Y = 1, Z = 2;

BOOST_AUTO_TEST_CASE( Axes )
{
    using SimpleMesh::Triangle;
    
    Vector3d o(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    Triangle xy(SimpleMesh::Test::oneTri(o,x,y)), yz(SimpleMesh::Test::oneTri(o,y,z)), zx(SimpleMesh::Test::oneTri(o,z,x));
    
    vector<const Triangle*> triangles;
    triangles.push_back(&yz);
    triangles.push_back(&zx);
    triangles.push_back(&xy);
    
    vector<RayIntersection> intersections(10);
    
    Vector2d rayHalfWidth(1e-10,1e-10);
    
    for (int rayAxis = 0; rayAxis < 3; rayAxis++)
    {
        Raycaster cast(rayAxis, rayHalfWidth, triangles);
        // All triangles go through the origin, so we get three hits here.
        BOOST_CHECK_EQUAL(cast.boxHits(Vector2d(0,0), intersections), 3);
        
        BOOST_CHECK_EQUAL(cast.boxHits(Vector2d(0.1, 0.1), intersections), 1);
        BOOST_CHECK_EQUAL(intersections[0].triangle(), triangles[rayAxis]);
        intersections[0] = RayIntersection();
        
        // this ray is outside the triangle but inside the bounding box.
        BOOST_CHECK_EQUAL(cast.boxHits(Vector2d(0.9, 0.9), intersections), 1);
        BOOST_CHECK_EQUAL(intersections[0].triangle(), triangles[rayAxis]);
        intersections[0] = RayIntersection();
        
        BOOST_CHECK_EQUAL(cast.boxHits(Vector2d(-0.1, 0.1), intersections), 0);
    }
}
































