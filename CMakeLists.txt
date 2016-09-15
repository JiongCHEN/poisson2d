project(poisson2d)
cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_CXX_FLAGS "-std=c++0x")

include_directories(${PROJECT_SOURCE_DIR})

add_subdirectory(examples)
add_subdirectory(src)
