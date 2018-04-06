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

add_executable(testSimpleMesh
    testSimpleMesh.cpp
    ../SimpleMesh.h
    ../SimpleMesh.cpp
    ../Search.h
    ../Search.cpp
)

add_executable(testLine
    testLine.cpp
    ../SimpleMesh.h
    ../SimpleMesh.cpp
    ../Search.h
    ../Search.cpp
)

macro(setupUnitTest executableName)
    target_link_libraries(${executableName} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} utility)
	add_test(${executableName} ${executableName})
endmacro()

setupUnitTest(testSimpleMesh)
setupUnitTest(testLine)


