SET(SRCS main.cpp readMesh.cpp readMesh.hpp buildM.cpp buildM.hpp buildS.cpp buildS.hpp buildDM.cpp buildDM.hpp Mesh2D.cpp Mesh2D.hpp buildFlux.hpp buildFlux.cpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})
