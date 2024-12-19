cmake_minimum_required(VERSION 3.10)

# Название проекта
project(TRETON)

# Установить стандарт C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Указать исходные файлы (в корневой папке)
set(SOURCES
    Main.cpp
    AMUV.cpp
	ANU.cpp
    CpSodium.cpp
	EntV.cpp
	Heat.cpp
	Hydro.cpp
	Matrix.cpp
	PRV.cpp
	SodiumProp.cpp
	Thechyco.cpp
	ThechycoGlobalVar.cpp
)

# Создание исполняемого файла
add_executable(${PROJECT_NAME} ${SOURCES})

# Указать директорию с заголовочными файлами
target_include_directories(${PROJECT_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/headers")