SET(SRCS
main.cpp
generateMesh.hpp generateMesh.cpp
./mesh/Mesh.cpp ./mesh/Mesh.hpp  ./mesh/displayMesh.cpp ./mesh/displayMesh.hpp
./matrices/buildM.cpp ./matrices/buildM.hpp ./matrices/buildS.cpp ./matrices/buildS.hpp ./matrices/buildMatrix.cpp ./matrices/buildMatrix.hpp ./matrices/matrix.hpp
./flux/buildFlux.cpp ./flux/buildFlux.hpp ./flux/flux.cpp ./flux/flux.hpp ./flux/phiPsi.cpp ./flux/phiPsi.hpp
./source/source.cpp ./source/source.hpp
./solver/timeInteg.cpp ./solver/timeInteg.hpp ./solver/field.hpp ./solver/rungeKutta.cpp ./solver/rungeKutta.hpp
./error/computeError.hpp ./error/computeNorm.cpp ./error/computeNorm.hpp ./error/computeError.cpp ./error/computeNormBis.cpp ./error/computeNormBis.hpp
./params/Params.hpp ./params/Params.cpp ./params/ibvFunction.hpp ./params/ibvFunction.cpp
./utils/utils.hpp ./utils/utils.cpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})
