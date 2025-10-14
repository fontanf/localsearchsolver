cmake_minimum_required(VERSION 3.15.0)

project(LocalSearchSolver LANGUAGES CXX)

# Build options.
option(LOCALSEARCHSOLVER_BUILD_TEST "Build unit tests" ON)
option(LOCALSEARCHSOLVER_BUILD_EXAMPLES "Build examples" ON)

# Avoid FetchContent warning.
cmake_policy(SET CMP0135 NEW)

# Require C++14.
set(CMAKE_CXX_STANDARD 14)

# Enable output of compile commands during generation.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add sub-directories.
add_subdirectory(extern)
add_subdirectory(src)
if(LOCALSEARCHSOLVER_BUILD_TEST)
    add_subdirectory(test)
endif()
if(LOCALSEARCHSOLVER_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
