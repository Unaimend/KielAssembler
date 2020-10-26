cmake_minimum_required(VERSION 3.16)
project(KielAssembler)

# compiler options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS  "-Wall -pthread"  )
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}" )

# prepare sources
set(SOURCE
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

set(HEADERS)

# add subdirectories
add_subdirectory(lib)
add_subdirectory(deBruijn)
add_subdirectory(dataStructures)

add_executable(${PROJECT_NAME} ${SOURCE} ${HEADERS})

target_compile_options(KielAssembler PRIVATE -O3)
