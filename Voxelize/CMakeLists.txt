project(Voxelize)
# ------------------ Boilerplate

cmake_minimum_required(VERSION 2.4.5)
if (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif (COMMAND cmake_policy)

include_directories(
    /usr/local/include
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../utility
)

add_library( Voxelize
    ClipTriangle.cpp
    ClipTriangle.h
    KDTree.cpp
    KDTree.h
    Raycaster.cpp
    Raycaster.h
    RowVoxelizer.cpp
    RowVoxelizer.h
    Voxelize.cpp
    Voxelize.h
)

target_link_libraries( Voxelize
    utility
    SimpleMesh
)

add_subdirectory(tests)
