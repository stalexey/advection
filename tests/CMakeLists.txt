cmake_minimum_required(VERSION 3.5)

add_compile_options(-std=c++14)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories("/opt/local/include/libomp")
add_definitions("-Xpreprocessor -fopenmp")

link_directories("/opt/local/lib/libomp")
add_link_options("-lomp")

add_executable(runAdvection runAdvection.cpp)
add_executable(runInterpolation runInterpolation.cpp)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plotAdvection.sh
    ${CMAKE_CURRENT_BINARY_DIR}/plotAdvection.sh COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/displayAdvection.gnu
    ${CMAKE_CURRENT_BINARY_DIR}/displayAdvection.gnu COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plotInterpolation.sh
    ${CMAKE_CURRENT_BINARY_DIR}/plotInterpolation.sh COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/displayInterpolation.gnu
    ${CMAKE_CURRENT_BINARY_DIR}/displayInterpolation.gnu COPYONLY)