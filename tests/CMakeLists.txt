add_executable(runAdvection runAdvection.cpp)
add_executable(runInterpolation runInterpolation.cpp)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plotAdvection.sh
    ${CMAKE_CURRENT_BINARY_DIR}/plotAdvection.sh COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/advection.gnu
    ${CMAKE_CURRENT_BINARY_DIR}/advection.gnu COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plotInterpolation.sh
    ${CMAKE_CURRENT_BINARY_DIR}/plotInterpolation.sh COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/interpolation.gnu
    ${CMAKE_CURRENT_BINARY_DIR}/interpolation.gnu COPYONLY)