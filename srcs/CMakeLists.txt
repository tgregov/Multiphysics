SET(SRCS main.cpp Mesh2D.cpp Mesh2D.hpp buildM.cpp buildM.hpp buildS.cpp buildS.hpp buildFlux.hpp buildFlux.cpp displayMesh.cpp displayMesh.hpp timeInteg.cpp timeInteg.hpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})
