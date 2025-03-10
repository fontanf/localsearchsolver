cmake_minimum_required(VERSION 3.15.0)

project(LocalSearchSolver LANGUAGES CXX)

# Options.
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
add_subdirectory(test)
if(LOCALSEARCHSOLVER_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
