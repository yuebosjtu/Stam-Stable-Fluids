cmake_minimum_required(VERSION 3.16)
project(StamStableFluids)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use the included Eigen library
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Include directories
include_directories(include)
include_directories(visualize)

# Add main executable
add_executable(main 
    main.cpp
    src/fluid_simulator.cpp
    visualize/colorramp.cpp
)

# Set output directory to the source directory (same as main.cpp)
set_target_properties(main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_SOURCE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}"
)

# Set compiler flags for optimization
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(MSVC)
        target_compile_options(main PRIVATE /O2)
    else()
        target_compile_options(main PRIVATE -O3 -march=native)
    endif()
endif()

# Print build information
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Using included Eigen library")
message(STATUS "Target: main")
