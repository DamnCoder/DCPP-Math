cmake_minimum_required(VERSION 3.1)

#[PROJECT_NAME]
PROJECT(DCPPMath)

#[BUILD_TYPE]
SET(CMAKE_BUILD_TYPE Debug)

# Set the output folders where the program will be created
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(SRC_DIR ${PROJECT_SOURCE_DIR}/project)

#[PRJ_INCLUDE]
INCLUDE_DIRECTORIES(${SRC_DIR}/include)
INCLUDE_DIRECTORIES(${SRC_DIR}/include/math)

#[PRJ_HEADER_FILES]
SET(HEADERS
	${SRC_DIR}/include/math/math.h
	${SRC_DIR}/include/math/matrix.inl
	${SRC_DIR}/include/math/matrix.h
	${SRC_DIR}/include/math/quaternion.inl
	${SRC_DIR}/include/math/quaternion.h
	${SRC_DIR}/include/math/vector2.inl
	${SRC_DIR}/include/math/vector2.h
	${SRC_DIR}/include/math/vector3.inl
	${SRC_DIR}/include/math/vector3.h
	${SRC_DIR}/include/math/vector4.inl
	${SRC_DIR}/include/math/vector4.h
	${SRC_DIR}/include/math/vector.h
	${SRC_DIR}/include/math/grid.h
	${SRC_DIR}/include/math/matrixstacks.h
)

SET(SOURCES
	${SRC_DIR}/src/math/matrixstacks.cpp
)

# Generate the static library from the sources
ADD_LIBRARY(${PROJECT_NAME} STATIC ${SOURCES} ${HEADERS})

SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)
 
# Set the location for library installation
# Use "sudo make install" to apply
INSTALL(TARGETS ${PROJECT_NAME} DESTINATION ${LIBRARY_OUTPUT_PATH})
INSTALL(DIRECTORY ${SRC_DIR}/include DESTINATION ${LIBRARY_OUTPUT_PATH})
