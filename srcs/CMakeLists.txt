SET(SRCS main.cpp readMesh.cpp readMesh.hpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})