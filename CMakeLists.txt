cmake_minimum_required(VERSION 3.10)
project(UnstableInvasion VERSION 0.1.0)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(MPI REQUIRED)

message(INFO "
To run MPI binaries:
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS
"
)

link_libraries(MPI::MPI_CXX)

add_subdirectory(3D/Flat/Bystander)
add_subdirectory(3D/Flat/Two-Species)
add_subdirectory(2D/Flat/Bystander)
add_subdirectory(2D/Flat/Two-Species)
