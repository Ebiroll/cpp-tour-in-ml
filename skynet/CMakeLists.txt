cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(stdexecExample)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 20)

# Set the directory for the CPM.cmake file
set(CPM_PATH "${CMAKE_SOURCE_DIR}/cmake/get_cpm.cmake")




# Ensure CPM.cmake is downloaded if not already available
if(NOT EXISTS ${CPM_PATH})
    message(STATUS "Downloading get_cpm.cmake to ${CPM_PATH}")
    file(DOWNLOAD 
        https://github.com/cpm-cmake/CPM.cmake/releases/latest/download/get_cpm.cmake
        ${CPM_PATH}
        EXPECTED_HASH SHA256=8aa4ec46e3673e53387a8b9ea53b12e7d953958a5822c678a5c984cbd6ac1328
    )
endif()

# Include CPM
include(${CPM_PATH})


# Get CPM
# For more information on how to add CPM to your project, see: https://github.com/cpm-cmake/CPM.cmake#adding-cpm
include(cmake/CPM.cmake)
# wget -O cmake/CPM.cmake https://github.com/cpm-cmake/CPM.cmake/releases/latest/download/get_cpm.cmake


CPMAddPackage(
  NAME stdexec
  GITHUB_REPOSITORY NVIDIA/stdexec
  GIT_TAG main # This will always pull the latest code from the `main` branch. You may also use a specific release version or tag
)

include_directories(.)

add_executable(skynet skynet.cpp)

add_executable(nbody nbody_sim.cpp)


target_link_libraries(nbody STDEXEC::stdexec)

target_link_libraries(skynet STDEXEC::stdexec)
