SET(SRCS
main.cpp
./mesh/Mesh.cpp ./mesh/Mesh.hpp  ./mesh/displayMesh.cpp ./mesh/displayMesh.hpp
./matrices/buildM.cpp ./matrices/buildM.hpp ./matrices/buildS.cpp ./matrices/buildS.hpp ./matrices/buildMatrix.cpp ./matrices/buildMatrix.hpp ./matrices/matrix.hpp
./flux/buildFlux.cpp ./flux/buildFlux.hpp
./solver/timeInteg.cpp ./solver/timeInteg.hpp ./solver/field.hpp ./solver/RungeKutta.cpp ./solver/RungeKutta.hpp
./params/Params.hpp ./params/Params.cpp
./utils/utils.hpp ./utils/utils.cpp
./write/write.hpp ./write/write.cpp
./physics/shallowWater/flux.hpp ./physics/shallowWater/flux.cpp
./physics/shallowWater/phiPsi.hpp ./physics/shallowWater/phiPsi.cpp
./physics/shallowWater/boundaryCondition.hpp ./physics/shallowWater/boundaryCondition.cpp
./physics/shallowWater/writer.hpp ./physics/shallowWater/writer.cpp
./physics/shallowWater/source.hpp ./physics/shallowWater/source.cpp
./physics/linShallowWater/flux.hpp ./physics/linShallowWater/flux.cpp
./physics/linShallowWater/phiPsi.hpp ./physics/linShallowWater/phiPsi.cpp
./physics/linShallowWater/source.hpp ./physics/linShallowWater/source.cpp
./physics/linShallowWater/boundaryCondition.hpp ./physics/linShallowWater/boundaryCondition.cpp
./physics/linShallowWater/writer.hpp ./physics/linShallowWater/writer.cpp
./physics/transport/flux.hpp ./physics/transport/flux.cpp
./physics/transport/phiPsi.hpp ./physics/transport/phiPsi.cpp
./physics/transport/boundaryCondition.hpp ./physics/transport/boundaryCondition.cpp
./physics/transport/writer.hpp ./physics/transport/writer.cpp
./physics/fluxes.hpp ./physics/phiPsis.hpp ./physics/boundaryConditions.hpp ./physics/sources.hpp
./physics/ibcFunction.hpp ./physics/writers.hpp
./physics/meanPhiPsi.cpp ./physics/meanPhiPsi.hpp
./physics/commonBC.cpp ./physics/commonBC.hpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})
