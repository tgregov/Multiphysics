ADD_SUBDIRECTORY(dG/Mesh)

SET(SRCS
main.cpp)

ADD_EXECUTABLE(DGangoFEtt ${SRCS})
TARGET_LINK_LIBRARIES(DGangoFEtt PUBLIC dgMesh)
SET_TARGET_PROPERTIES(DGangoFEtt PROPERTIES COMPILE_FLAGS "-Wall -Wextra -pedantic-errors -Wold-style-cast -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wshadow")
