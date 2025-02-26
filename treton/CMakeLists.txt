cmake_minimum_required(VERSION 3.10)

project(treton VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(src_treton_dir "${CMAKE_CURRENT_SOURCE_DIR}/src")

set(src_obj "${src_treton_dir}/CoolantMaterials.cpp"
			"${src_treton_dir}/CRD_test.cpp"
			"${src_treton_dir}/Heat.cpp"
			"${src_treton_dir}/Hydro.cpp"
			"${src_treton_dir}/Matrix.cpp"
			"${src_treton_dir}/Thechyco.cpp"
			"${src_treton_dir}/GlobalVar.cpp"
)


add_library(object_treton OBJECT ${src_obj})
target_include_directories(object_treton PUBLIC "include")

