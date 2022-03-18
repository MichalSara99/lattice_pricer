cmake_minimum_required(VERSION 3.0.0)
project(lattice_pricer VERSION 0.1.0)

set(SOURCES ${PROJECT_SOURCE_DIR}/main.cpp)

include(CTest)
enable_testing()

add_subdirectory(unit_tests)

include_directories(${PROJECT_SOURCE_DIR}/include)
add_executable(lattice_pricer main.cpp)

add_library(lattice_pricer_lib STATIC ${SOURCES})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
