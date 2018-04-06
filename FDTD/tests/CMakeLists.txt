
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
    cmake_policy( SET CMP0003 NEW )  
endif()

include_directories(
    ${FDTD_SOURCE_DIR}
    ${FDTD_SOURCE_DIR}/..                 # project root
    ${FDTD_SOURCE_DIR}/../utility
)

set(Boost_ADDITIONAL_VERSIONS "1.41" "1.41.0")
#set(Boost_USE_STATIC_LIBS ON)
find_package( Boost COMPONENTS unit_test_framework REQUIRED  )
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

include_directories(/usr/local/include)

if (Boost_FOUND)
    message("Boost found, include dirs ${Boost_INCLUDE_DIRS}.")
    message("Libraries ${Boost_LIBRARIES}.")
    message("Test library ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}.")
endif (Boost_FOUND)

add_executable(testPMLParameters
    testPMLParameters.cpp
    ${FDTD_SOURCE_DIR}/PMLParameters.cpp
    ${FDTD_SOURCE_DIR}/PMLParameters.h
    ${FDTD_SOURCE_DIR}/PhysicalConstants.cpp
    ${FDTD_SOURCE_DIR}/PhysicalConstants.h
)

add_executable(testPMLField
    testPMLField.cpp
    ${FDTD_SOURCE_DIR}/PMLField.cpp
    ${FDTD_SOURCE_DIR}/PMLField.h
    ${FDTD_SOURCE_DIR}/PMLParameters.cpp
    ${FDTD_SOURCE_DIR}/PMLParameters.h
    ${FDTD_SOURCE_DIR}/PhysicalConstants.cpp
    ${FDTD_SOURCE_DIR}/PhysicalConstants.h
    ${FDTD_SOURCE_DIR}/Extents.cpp
    ${FDTD_SOURCE_DIR}/Extents.h
    ${FDTD_SOURCE_DIR}/YeeUtilities.h
    ${FDTD_SOURCE_DIR}/YeeUtilities.cpp
)
target_link_libraries(testPMLField
    utility
)

add_executable(testExtents
    testExtents.cpp
    ${FDTD_SOURCE_DIR}/Extents.cpp
    ${FDTD_SOURCE_DIR}/Extents.h
    ${FDTD_SOURCE_DIR}/YeeUtilities.h
    ${FDTD_SOURCE_DIR}/YeeUtilities.cpp
)
target_link_libraries(testExtents utility)


add_executable(testConstitutiveTensorField
    testConstitutiveTensorField.cpp
    ${FDTD_SOURCE_DIR}/ConstitutiveTensorField.cpp
    ${FDTD_SOURCE_DIR}/ConstitutiveTensorField.h
    ${FDTD_SOURCE_DIR}/Extents.cpp
    ${FDTD_SOURCE_DIR}/Extents.h
    ${FDTD_SOURCE_DIR}/OutputDirectory.cpp
    ${FDTD_SOURCE_DIR}/OutputDirectory.h
    ${FDTD_SOURCE_DIR}/UserPreferences.cpp
    ${FDTD_SOURCE_DIR}/UserPreferences.h
    ${FDTD_SOURCE_DIR}/YeeUtilities.cpp
    ${FDTD_SOURCE_DIR}/YeeUtilities.h
)
target_link_libraries(testConstitutiveTensorField
    utility
)

add_executable(testOrientedVoxels
    testOrientedVoxels.cpp
    ${FDTD_SOURCE_DIR}/OrientedVoxels.cpp
    ${FDTD_SOURCE_DIR}/OrientedVoxels.h
    ${FDTD_SOURCE_DIR}/VoxelizePolyhedra.cpp
    ${FDTD_SOURCE_DIR}/VoxelizePolyhedra.h
    ${FDTD_SOURCE_DIR}/Extents.cpp
    ${FDTD_SOURCE_DIR}/Extents.h
    ${FDTD_SOURCE_DIR}/UserPreferences.cpp
    ${FDTD_SOURCE_DIR}/UserPreferences.h
    ${FDTD_SOURCE_DIR}/YeeUtilities.cpp
    ${FDTD_SOURCE_DIR}/YeeUtilities.h
)
target_link_libraries(testOrientedVoxels
    Voxelize
    utility
)

add_executable(testPermittivitySensitivity
    testPermittivitySensitivity.cpp
    ${FDTD_SOURCE_DIR}/Extents.cpp
    ${FDTD_SOURCE_DIR}/Extents.h
    ${FDTD_SOURCE_DIR}/OutputDirectory.cpp
    ${FDTD_SOURCE_DIR}/OutputDirectory.h
    ${FDTD_SOURCE_DIR}/PMLParameters.cpp
    ${FDTD_SOURCE_DIR}/PMLParameters.h
    ${FDTD_SOURCE_DIR}/TensorGridConstants.cpp
    ${FDTD_SOURCE_DIR}/TensorGridConstants.h
    ${FDTD_SOURCE_DIR}/UserPreferences.cpp
    ${FDTD_SOURCE_DIR}/UserPreferences.h
    ${FDTD_SOURCE_DIR}/YeeUtilities.cpp
    ${FDTD_SOURCE_DIR}/YeeUtilities.h
)
target_link_libraries(testPermittivitySensitivity
    Voxelize
    utility
)

add_executable(testTensorGridConstants
    testTensorGridConstants.cpp
    ${FDTD_SOURCE_DIR}/Precision.h
    ${FDTD_SOURCE_DIR}/PrecisionGeometry.h
    ${FDTD_SOURCE_DIR}/PrecisionRationalFunction.h
    ${FDTD_SOURCE_DIR}/TensorGridConstants.cpp
    ${FDTD_SOURCE_DIR}/TensorGridConstants.h
    ${FDTD_SOURCE_DIR}/UserPreferences.cpp
    ${FDTD_SOURCE_DIR}/UserPreferences.h
    ${FDTD_SOURCE_DIR}/YeeUtilities.cpp
    ${FDTD_SOURCE_DIR}/YeeUtilities.h
)

target_link_libraries(testTensorGridConstants
    utility
)

macro(setupUnitTest executableName)
    target_link_libraries(${executableName} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
	add_test(${executableName} ${executableName})
endmacro()

setupUnitTest(testPMLParameters)
setupUnitTest(testPMLField)
setupUnitTest(testConstitutiveTensorField)
setupUnitTest(testExtents)
setupUnitTest(testOrientedVoxels)
setupUnitTest(testTensorGridConstants)
setupUnitTest(testPermittivitySensitivity)


