cmake_minimum_required(VERSION 3.14)
project(simple_http_server)

# Enable C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

# Configure Beman execution library (new repository)
FetchContent_Declare(
    execution
    GIT_REPOSITORY https://github.com/bemanproject/execution
    GIT_TAG        main  # Use the appropriate branch or tag (e.g., main, v1.0.0, etc.)
)

# Configure Beman net (networking library)
FetchContent_Declare(
    net
    GIT_REPOSITORY https://github.com/bemanproject/net
    GIT_TAG        main  # Verify this tag matches the version you need
)

# Fetch both repositories
FetchContent_MakeAvailable(execution net)

include_directories(.)

add_executable(simple_http_server simple_http_server.cpp)


# Include the Beman libraries' headers
target_include_directories(simple_http_server PRIVATE
    ${execution_SOURCE_DIR}/include
    ${net_SOURCE_DIR}/include
    ${net_SOURCE_DIR}/examples
)



# Link threading library (if needed by the Beman libraries)
#target_link_libraries(simple_http_server  Threads::Threads)

