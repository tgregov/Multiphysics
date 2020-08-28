# find gmsh-sdk
# gmsh.h
FIND_PATH(GMSH_INCLUDE_DIRS NAMES "gmsh.h")
MESSAGE(STATUS "Found GMSH: " ${GMSH_INCLUDE_DIRS})
if(NOT GMSH_INCLUDE_DIRS)
    MESSAGE(FATAL_ERROR "gmsh.h not found!")
ENDIF()
INCLUDE_DIRECTORIES(SYSTEM ${GMSH_INCLUDE_DIRS})

# libgmsh.so
FIND_LIBRARY(GMSH_LIBRARIES gmsh)
MESSAGE(STATUS "Found GMSH: " ${GMSH_LIBRARIES})
IF(NOT GMSH_LIBRARIES)
    MESSAGE(FATAL_ERROR "gmsh library not found!")
ENDIF()

SET(LIBMESHSRCS
Mesh.hpp Mesh.inl Mesh.cpp
Entity.hpp PhysicalGroup.hpp ElementProperty.hpp Element.hpp Face.hpp
Entity.inl PhysicalGroup.inl Element.inl Face.inl)

ADD_LIBRARY(dgMesh SHARED ${LIBMESHSRCS})
GENERATE_EXPORT_HEADER(dgMesh
                       EXPORT_MACRO_NAME DG_MESH_API
					   EXPORT_FILE_NAME dgMesh_export.h)
TARGET_INCLUDE_DIRECTORIES(dgMesh SYSTEM
                           PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE ${GMSH_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(dgMesh PRIVATE ${GMSH_LIBRARIES})
SET_TARGET_PROPERTIES(dgMesh PROPERTIES COMPILE_FLAGS "-Wall -Wextra -pedantic-errors -Wold-style-cast -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wshadow")
