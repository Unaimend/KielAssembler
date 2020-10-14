cmake_minimum_required(VERSION 3.17)
project(KielAssembler)

set(CMAKE_CXX_STANDARD 17)
add_executable(KielAssembler main.cpp deBruijn/DeBruijnGraph.cpp deBruijn/DeBruijnGraph.h lib/catch.hpp lib/bioio.hpp lib/cxxopts.h)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-z,--stack_size=10000000000")
#target_compile_options(KielAssembler PRIVATE -O3)
SET(CMAKE_CXX_FLAGS  "-Wall")
