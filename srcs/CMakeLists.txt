SET(SRCS main.cpp readMesh.cpp readMesh.hpp buildM.cpp buildM.hpp buildS.cpp buildS.hpp Mesh2D.cpp Mesh2D.hpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})