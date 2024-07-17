cmake_minimum_required(VERSION 3.10)
project(ParallelSolvers)

# Set the C and C++ compilers
set(CMAKE_CXX_COMPILER icpx)

# Set the include directories
include_directories("${MKLROOT}/include" "include")

# Set the compilation flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w -m64 -lm -qopenmp")

set(MKLROOT $ENV{MKLROOT})

# Set the link libraries
set(LINK_LIBS
    -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group
    iomp5
    pthread
    m
    dl
)

# Add executable targets
add_executable(direct_solver src/direct_solver.cpp src/direct_real.cpp src/direct_complex.cpp)
target_link_libraries(direct_solver ${LINK_LIBS})

add_executable(iterative_solver src/iterative_solver.cpp src/iterative_real.cpp src/iterative_complex.cpp)
target_link_libraries(iterative_solver ${LINK_LIBS})

add_executable(direct_solver_ilp64 src/direct_solver_ilp64.cpp src/direct_real_ilp64.cpp src/direct_complex_ilp64.cpp)
target_link_libraries(direct_solver_ilp64 ${LINK_LIBS})

add_executable(hybrid_solver src/hybrid_solver.cpp src/direct_real.cpp src/direct_complex.cpp src/iterative_real.cpp src/iterative_complex.cpp)
target_link_libraries(hybrid_solver ${LINK_LIBS})

# add_executable(hybrid_solver src/hybrid_solver.c)
# target_link_libraries(hybrid_solver ${LINK_LIBS})