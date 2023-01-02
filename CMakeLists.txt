cmake_minimum_required(VERSION 3.10)
project(LibACE)

set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)

add_library(ace
    source/spherical_bessel.cpp
    source/spherical_harmonic.cpp
    source/acecore/simpleprodbasis.cpp
)
target_include_directories(ace PUBLIC source)

# PYBIND
add_subdirectory(external/pybind11)
pybind11_add_module(ace-pybind MODULE source/pybind/ace-pybind.cpp)
target_link_libraries(ace-pybind PUBLIC ace)
set_target_properties(ace-pybind PROPERTIES OUTPUT_NAME ace)

# CXXWRAP
execute_process(COMMAND julia -e "using CxxWrap; print(CxxWrap.prefix_path())"
                OUTPUT_VARIABLE CXXWRAP_PREFIX_PATH)
list(APPEND CMAKE_PREFIX_PATH ${CXXWRAP_PREFIX_PATH})
find_package(JlCxx)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")
message(STATUS "Found JlCxx at ${JlCxx_location}")
add_library(ace-cxxwrap SHARED source/cxxwrap/ace-cxxwrap.cpp)
target_link_libraries(ace-cxxwrap PRIVATE JlCxx::cxxwrap_julia)
target_link_libraries(ace-cxxwrap PRIVATE JlCxx::cxxwrap_julia_stl)
target_link_libraries(ace-cxxwrap PUBLIC ace)
