###############################################
## This is general CMakeLists.txt file       ##
## which is needed to link your own program  ##
## with Pteros libraries                     ##
###############################################

cmake_minimum_required(VERSION 3.1)

project(ff_compute)
SET(CMAKE_BUILD_TYPE Release)

#Sets optmization in Release build
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
# Set specific options for GCC if used
if(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed -std=c++11")
endif()
# Strip libraries and executables to decrease size in Release builds
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -s")
    SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -s")
ENDIF()


find_package(pteros REQUIRED)

add_executable(ff_compute main.cpp FFMaps.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(ff_compute PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(ff_compute PRIVATE pteros::pteros pteros::pteros_analysis)
