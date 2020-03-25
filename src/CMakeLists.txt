add_executable(runAdvection runAdvection.cpp)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plotAdvection.sh
    ${CMAKE_CURRENT_BINARY_DIR}/plotAdvection.sh COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/advection.gnu
    ${CMAKE_CURRENT_BINARY_DIR}/advection.gnu COPYONLY)