cmake_minimum_required(VERSION 3.0)
project(
    pybp
    LANGUAGES CXX
)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_BUILD_TYPE Release)
add_subdirectory(pybind11)
pybind11_add_module(pybp main.cc)
