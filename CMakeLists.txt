cmake_minimum_required(VERSION 3.17)
project(KielAssembler)

set(CMAKE_CXX_STANDARD 17)

add_executable(KielAssembler main.cpp deBruijn/DeBruijnGraph.cpp deBruijn/DeBruijnGraph.h)