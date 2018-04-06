/*
 *  testVoxelize.cpp
 *  Trogdor6
 *
 *  Created by Paul C Hansen on 7/23/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */


#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testVoxelize

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Voxelize.h"


BOOST_AUTO_TEST_CASE( GetAndSet )
{
    Rect3d bounds(-1,-2,-3,1,2,3);
    Rect3i voxels(-1,-1,-1,1,1,1);
    
    Voxelize vox(bounds, voxels);
    BOOST_CHECK(vox.bounds() == bounds);
    BOOST_CHECK(vox.voxelBounds() == voxels);
    BOOST_CHECK(vox.voxelSize() == Vector3d(2.0/3.0, 4.0/3.0, 6.0/3.0));
    BOOST_CHECK(vox.voxelHalfWidth(0) == Vector2d(2.0/3.0, 1.0));
    BOOST_CHECK(vox.voxelHalfWidth(1) == Vector2d(1.0, 1.0/3.0));
    BOOST_CHECK(vox.voxelHalfWidth(2) == Vector2d(1.0/3.0, 2.0/3.0));
    BOOST_CHECK(vox.nonSymmetricDimensions() == Vector3b(true,true,true));
    
    BOOST_CHECK(vox.voxel(-1,-1,-1) == Rect3d(-1, -2, -3,
        -1+2.0/3, -2 + 4.0/3, -3 + 6.0/3));
}
