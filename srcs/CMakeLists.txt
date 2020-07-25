ADD_SUBDIRECTORY(dgMesh)

SET(SRCS
main.cpp)

ADD_EXECUTABLE(dgExe ${SRCS})
TARGET_LINK_LIBRARIES(dgExe PUBLIC dgMesh)
SET_TARGET_PROPERTIES(dgExe PROPERTIES COMPILE_FLAGS "-Wall -Wextra -pedantic-errors -Wold-style-cast -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wshadow")