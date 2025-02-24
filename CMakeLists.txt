cmake_minimum_required(VERSION 3.10)

# Название проекта
project(TRETON)

# Установить стандарт C++
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Указать исходные файлы (в корневой папке)
set(SOURCES
    Main.cpp
	Heat.cpp
	Hydro.cpp
	Matrix.cpp
	Thechyco.cpp
	ThechycoGlobalVar.cpp
	CoolantMaterials.cpp
	CRD_test.cpp
)

# Создание исполняемого файла
add_executable(${PROJECT_NAME} ${SOURCES})

# Указать директорию с заголовочными файлами
target_include_directories(${PROJECT_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/headers")

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(TRETON PRIVATE -O3 -ffast-math -march=native -flto)
endif()


if(DEFINED FHI0)
    target_compile_definitions(TRETON PRIVATE FHI0=${FHI0})
endif()


if(DEFINED TEST_CRD)
    target_compile_definitions(TRETON PRIVATE TEST_CRD=${TEST_CRD})
endif()