
cmake_minimum_required(VERSION 3.15...3.27)

project(${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

add_subdirectory(pybind11)
add_subdirectory(matplotplusplus)
add_subdirectory(AudioFile)


python_add_library(_core MODULE src/main.cpp WITH_SOABI)
target_link_libraries(_core PRIVATE pybind11::headers)
target_link_libraries(_core PUBLIC matplot)
target_link_libraries(_core PRIVATE AudioFile)


target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS _core DESTINATION scikit_build_example)
