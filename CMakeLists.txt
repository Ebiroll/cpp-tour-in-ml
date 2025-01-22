cmake_minimum_required(VERSION 3.10)
project(TF_Sine_Example CXX)

if(POLICY CMP0167)
cmake_policy(SET CMP0167 NEW)
endif(POLICY CMP0167)

# Build debug
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_USE_RELATIVE_PATHS ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Adjust this path to your TensorFlow C++ installation
# For example, if you built from source and have the library in /usr/local/lib
set(TENSORFLOW_INCLUDE_DIR "tensorflow/bazel-bin/tensorflow/include")
#set(TENSORFLOW_LIB "tensorflow/bazel-bin/tensorflow/libtensorflow_cc.so")
#For docker file
set(TENSORFLOW_LIB "/usr/local/lib/python3.11/dist-packages/tensorflow/libtensorflow_cc.so.2")

# Add library search dir, bazel-bin/tensorflow/
#link_directories(tensorflow/bazel-bin/tensorflow)
#link_directories(bazel-bin/tensorflow)
link_directories(/usr/local/lib/python3.11/dist-packages/tensorflow)

set(FRAMEWORK_LIB "/usr/local/lib/python3.11/dist-packages/tensorflow/libtensorflow_framework.so.2")

link_directories(/usr/local/lib/python3.11/dist-packages/tensorflow)

# packages.
find_package(Boost REQUIRED COMPONENTS
    iostreams system filesystem
)

add_library(gnuplot_iostream INTERFACE)
target_include_directories(gnuplot_iostream INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
    $<INSTALL_INTERFACE:>)

target_link_libraries(gnuplot_iostream INTERFACE
    Boost::iostreams
    Boost::system
    Boost::filesystem
)

set(GNUPLOT_IOSTREAM_LIBRARIES Boost::iostreams Boost::system Boost::filesystem)


#add include dir
include_directories(gnuplot-iostream)

include_directories(tensorflow)

include_directories(/usr/local/lib/python3.11/dist-packages/tensorflow/include/)

# For gnuplot-iostream:
find_package(PkgConfig REQUIRED)
#pkg_check_modules(GNUPLOT_IOSTREAM REQUIRED gnuplot-iostream)
# find gnuplot
find_package(Gnuplot REQUIRED)




include_directories(${TENSORFLOW_INCLUDE_DIR})
add_executable(tf_sine_example sin/main.cpp)
target_link_libraries(tf_sine_example PRIVATE ${TENSORFLOW_LIB} ${FRAMEWORK_LIB}  ${GNUPLOT_IOSTREAM_LIBRARIES})

add_executable(label_image label_image/main.cc)
target_link_libraries(label_image PRIVATE ${TENSORFLOW_LIB} ${FRAMEWORK_LIB}  ${GNUPLOT_IOSTREAM_LIBRARIES})

# Fetch Eigen with ExternalProject
##########################################

include(ExternalProject)

set(EIGEN_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-install/")

ExternalProject_Add(
    eigen
    URL  https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-src"
    BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-build"
    INSTALL_DIR "${EIGEN_INSTALL_DIR}"
    CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    -DCMAKE_BUILD_TYPE=Release
)

file(MAKE_DIRECTORY ${EIGEN_INSTALL_DIR}/include)  # avoid race condition

add_library(eigenlib INTERFACE IMPORTED GLOBAL)
add_dependencies(eigenlib eigen)
target_compile_features(eigenlib INTERFACE cxx_std_14)

set_target_properties(eigenlib PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES ${EIGEN_INSTALL_DIR}/include/eigen3
)
#####################################



#add_library(linalg INTERFACE)

#target_include_directories(linalg
#    INTERFACE 
#    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#    $<INSTALL_INTERFACE:include>)

add_executable(simple_nn eigen/simple_nn.cpp)


target_link_libraries(simple_nn INTERFACE eigenlib)





# If needed, explicitly link or include for gnuplot-iostream:
# e.g. target_include_directories(tf_sine_example PRIVATE /usr/include/gnuplot-iostream)


# Add this to a executable test-interactive gnuplot-iostream/example-interactive.cc
add_executable(test-interactive gnuplot-iostream/example-interactive.cc)
target_link_libraries(test-interactive PRIVATE ${GNUPLOT_IOSTREAM_LIBRARIES})
target_include_directories(test-interactive PRIVATE ${GNUPLOT_IOSTREAM_INCLUDE_DIRS})

# Add this to a executable test-interactive gnuplot-iostream/example-interactive.cc
#add_executable(test gnuplot-iostream/example.cc)
#target_link_libraries(test PRIVATE ${GNUPLOT_IOSTREAM_LIBRARIES})
#target_include_directories(test PRIVATE ${GNUPLOT_IOSTREAM_INCLUDE_DIRS})