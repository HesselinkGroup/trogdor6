
project( Trogdor6TopLevel )

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

add_subdirectory( utility )
add_subdirectory( Voxelize )
add_subdirectory( tinyxml )
add_subdirectory( FDTD )
add_subdirectory( SimpleMesh )
