cmake_minimum_required(VERSION 3.10)
project(LibACE)

set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)

add_library(ace
    source/spherical_harmonic.cpp
)
target_include_directories(ace PUBLIC source)
