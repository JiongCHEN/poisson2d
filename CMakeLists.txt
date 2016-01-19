project(poisson-2d)
cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_CXX_FLAGS "-std=c++0x")

include_directories(${PROJECT_SOURCE_DIR})

add_executable(main util.h grid.h func.h problem.h main.cc)
