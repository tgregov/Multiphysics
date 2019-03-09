SET(SRCS main.cpp readMesh.cpp readMesh.hpp Mesh2D.cpp Mesh2D.hpp buildM.cpp buildM.hpp buildS.cpp buildS.hpp buildDM.cpp buildDM.hpp buildFlux.hpp buildFlux.cpp displayMesh.cpp displayMesh.hpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})
