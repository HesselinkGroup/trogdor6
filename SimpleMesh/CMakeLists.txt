project( SimpleMesh )

add_library( SimpleMesh
    Search.cpp
    Search.h
    SimpleMesh.h
    SimpleMesh.cpp
	XML.h
	XML.cpp
)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/../utility
	${CMAKE_CURRENT_SOURCE_DIR}/../tinyxml
)

target_link_libraries( SimpleMesh
	utility
	tinyxml
)

add_subdirectory(tests)



add_executable( xmltest
	main.cpp
)

target_link_libraries( xmltest
	SimpleMesh
)