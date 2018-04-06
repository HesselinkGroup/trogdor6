project( FDTD )
cmake_minimum_required(VERSION 2.4.5)

# ------------------ Boilerplate 
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

# ------------------ Set up the executable

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../utility
    ${CMAKE_CURRENT_SOURCE_DIR}/../tinyxml
)

add_executable( trogdor6
    AdjointGridOperations.cpp
    AdjointGridOperations.h
    AdjointRunlineCallbacks.h
    Arguments.cpp
    Arguments.h
    ConstitutiveTensorField.cpp
    ConstitutiveTensorField.h
    ElectromagneticFieldIndices.cpp
    ElectromagneticFieldIndices.h
    EncodeOutput.cpp
    EncodeOutput.h
    Extents.cpp
    Extents.h
    FieldEnum.h
    FillFactorReport.cpp
    FillFactorReport.h
    ForwardGridOperations.cpp
    ForwardGridOperations.h
    ForwardRunlineCallbacks.h
    GridFields.cpp
    GridFields.h
    GridOperations.cpp
    GridOperations.h
    main.cpp
    Operation.cpp # could be removed?
    Operation.h
    operations/HuygensSurface.cpp
    operations/HuygensSurface.h
    operations/OutputEHDB.cpp
    operations/OutputEHDB.h
    operations/Source.cpp
    operations/Source.h
    operations/AdjointDB.cpp
    operations/AdjointDB.h
    operations/AdjointEH.cpp
    operations/AdjointEH.h
    operations/ForwardDB.h
    operations/ForwardDB_2d.h
    operations/ForwardEH.cpp
    operations/ForwardEH.h
    operations/ForwardEH_0.cpp
    operations/ForwardEH_0.h
    operations/ForwardEH_i.cpp
    operations/ForwardEH_i.h
    operations/ForwardEH_j.cpp
    operations/ForwardEH_j.h
    operations/ForwardEH_templated.h
    operations/ForwardEH_templated-inl.h
    operations/ForwardNondispersiveEH.h
    operations/ForwardNondispersiveEH_2d.h
    operations/JMPML.h
    operations/PeriodicBC.h
    OrientationReport.cpp
    OrientationReport.h
    OrientedVoxels.cpp
    OrientedVoxels.h
    OutputDirectory.cpp
    OutputDirectory.h
    PhysicalConstants.cpp
    PhysicalConstants.h
    PMLField.cpp
    PMLField.h
    PMLParameters.cpp
    PMLParameters.h
    Precision.h
    PrecisionGeometry.h
    PrecisionRationalFunction.h
    SimulationDescription.cpp
    SimulationDescription.h
    Support.cpp
    Support.h
    TensorGridConstants.cpp
    TensorGridConstants.h
    UniversalRunlineCallbacks.h
    UserPreferences.cpp
    UserPreferences.h
    Version.h
    VoxelizedPartition.cpp
    VoxelizedPartition.h
    VoxelizePolyhedra.cpp
    VoxelizePolyhedra.h
    XMLParameterFile.cpp
    XMLParameterFile.h
    YeeUtilities.cpp
    YeeUtilities.h
)
add_custom_command(TARGET trogdor6 PRE_BUILD COMMAND touch Version.h)

find_package( Boost COMPONENTS program_options system filesystem)

if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message("\nBoost includes are at ${Boost_INCLUDE_DIRS}")
    target_link_libraries(trogdor6 ${Boost_LIBRARIES})
else (Boost_FOUND)
    message("\nCan't find Boost!")
endif (Boost_FOUND)

target_link_libraries(trogdor6
    ${Boost_LIBRARIES}
    tinyxml
    utility
    Voxelize
)

add_subdirectory(tests)
