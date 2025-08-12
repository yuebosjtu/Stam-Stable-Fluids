cmake_minimum_required(VERSION 3.16)
project(StamStableFluids)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use the included Eigen library
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Include directories
include_directories(include)

# Add main executable (original)
add_executable(main 
    main.cpp
    src/fluid_simulator.cpp
)

# Add circular source demo executable
add_executable(main_circular_source 
    main_circular_source.cpp
    src/fluid_simulator.cpp
)

# Set output directory to the source directory (same as main.cpp)
set_target_properties(main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_SOURCE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}"
)

set_target_properties(main_circular_source PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_SOURCE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}"
)

# Set compiler flags for optimization
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(MSVC)
        target_compile_options(main PRIVATE /O2)
        target_compile_options(main_circular_source PRIVATE /O2)
    else()
        target_compile_options(main PRIVATE -O3 -march=native)
        target_compile_options(main_circular_source PRIVATE -O3 -march=native)
    endif()
endif()

# Create output directory
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/fluid_output)

# Print build information
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Using included Eigen library")
message(STATUS "Targets: main, main_circular_source")

# Add custom target to run the circular source demo
add_custom_target(run_circular_source
    COMMAND main_circular_source
    DEPENDS main_circular_source
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running circular source fluid simulation demo"
)
