cmake_minimum_required(VERSION 3.10)
project(LibACE)

set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)

add_library(ace
    source/spherical_bessel.cpp
    source/spherical_harmonic.cpp
)
target_include_directories(ace PUBLIC source)

add_subdirectory(external/pybind11)
pybind11_add_module(ace-pybind MODULE source/pybind/ace-pybind.cpp)
target_link_libraries(ace-pybind PUBLIC ace)
set_target_properties(ace-pybind PROPERTIES OUTPUT_NAME ace)
