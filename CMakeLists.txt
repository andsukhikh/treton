cmake_minimum_required(VERSION 3.10)

project(treton VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory("${CMAKE_SOURCE_DIR}/treton")

add_executable(${PROJECT_NAME} "${CMAKE_SOURCE_DIR}/treton/src/Main.cpp")
target_link_libraries(${PROJECT_NAME} object_treton)


if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(${PROJECT_NAME} PUBLIC -O3 -ffast-math -march=native -flto)
endif()


if(DEFINED FHI0)
    target_compile_definitions(${PROJECT_NAME} PUBLIC FHI0=${FHI0})
endif()


if(DEFINED TEST_CRD)
    target_compile_definitions(${PROJECT_NAME} PUBLIC TEST_CRD=${TEST_CRD})
endif()