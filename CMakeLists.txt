cmake_minimum_required(VERSION 3.9)
project(BCHCodes)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

set(SOURCE_FILES src/main.cpp src/BCHCodes.cpp include/BCHCodes.h include/GF.h src/GF.cpp include/GaussianElimination.h src/GaussianElimination.cpp include/ReedMullerCodes.h src/ReedMullerCodes.cpp)
add_executable(BCHCodes ${SOURCE_FILES})