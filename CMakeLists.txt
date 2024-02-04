cmake_minimum_required(VERSION 3.15...3.27)
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)
  
# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Find OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Include Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

python_add_library(_core MODULE scdemon/py_se.cpp WITH_SOABI)
target_link_libraries(_core PRIVATE pybind11::headers)
target_link_libraries(_core PRIVATE Eigen3::Eigen)

target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

# Include the C++ headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
install(TARGETS _core DESTINATION ${SKBUILD_PROJECT_NAME})
