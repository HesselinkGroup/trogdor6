if ( COMMAND cmake_policy )
    cmake_policy( SET CMP0003 NEW )  
endif()

set(Boost_ADDITIONAL_VERSIONS "1.41" "1.41.0")
find_package( Boost COMPONENTS unit_test_framework REQUIRED  )
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

include_directories(
    /usr/local/include
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../../utility
)

add_executable(testClipTriangle
    testClipTriangle.cpp
)

add_executable(testClipTriangleJacobians
    testClipTriangleJacobians.cpp
)

add_executable(testRaycaster
    testRaycaster.cpp
)

add_executable(testRowVoxelizer
    testRowVoxelizer.cpp
)

add_executable(testVoxelizeBasics
    testVoxelize.cpp
)

add_executable(testFillFactors
    testFillFactors.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../RowVoxelizer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../RowVoxelizer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../Voxelize.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Voxelize.h
)

add_executable(testOrientation
    testOrientation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../RowVoxelizer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../RowVoxelizer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../Voxelize.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Voxelize.h
)

add_executable(testFillFactorJacobians
    testFillFactorJacobians.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../RowVoxelizer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../RowVoxelizer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../Voxelize.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Voxelize.h
)

add_executable(testOrientationJacobians
    testOrientationJacobians.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../RowVoxelizer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../RowVoxelizer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../Voxelize.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Voxelize.h
)

macro(setupUnitTest executableName)
    target_link_libraries(${executableName}
		${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
		SimpleMesh
		Voxelize)
	add_test(${executableName} ${executableName})
	# Uncomment to automatically run the tests during build.
	# If the test itself crashes then the executable is thrown out
	# and it becomes hard to debug the test.
    #add_custom_command( TARGET ${executableName} POST_BUILD
    #    COMMAND $<TARGET_FILE:${executableName}> --log_level=message)
endmacro()

setupUnitTest(testClipTriangle)
setupUnitTest(testClipTriangleJacobians)
setupUnitTest(testRaycaster)
setupUnitTest(testRowVoxelizer)
setupUnitTest(testVoxelizeBasics)
setupUnitTest(testFillFactors)
setupUnitTest(testOrientation)
setupUnitTest(testFillFactorJacobians)
setupUnitTest(testOrientationJacobians)

