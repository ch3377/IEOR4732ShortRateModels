cmake_minimum_required(VERSION 3.13)
project(pyPigs)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(src)
add_subdirectory(pylib)

include_directories(pylib/include)
include_directories(src/diffusions)
#include_directories(src/fft)
include_directories(src/rates)


pybind11_add_module(pypigs pypigs.cpp)
target_link_libraries(pypigs PRIVATE diffusions rates)
#link fft diffusions