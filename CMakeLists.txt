cmake_minimum_required(VERSION 3.16)
project(DACPO VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Windows specific settings
if(WIN32)
    add_compile_definitions(WIN32 _CRT_SECURE_NO_WARNINGS)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_definitions(CURL_STATICLIB)
    endif()
endif()

# Find required packages
find_package(PCL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenMP REQUIRED)



# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/Src/include
    ${CMAKE_SOURCE_DIR}/Src/include/3dparty
    ${PCL_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
)

# Add compile definitions for PCL
add_definitions(${PCL_DEFINITIONS})

# Collect source files
file(GLOB_RECURSE SOURCES
    "Src/*.cpp"
)

# Remove any main.cpp duplicates and ensure we have the right main
list(FILTER SOURCES EXCLUDE REGEX ".*main\\.cpp$")
list(APPEND SOURCES "Src/main.cpp")

# Collect header files
file(GLOB_RECURSE HEADERS
    "Src/include/*.h"
    "Src/include/3dparty/*.h"
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Link libraries
target_link_libraries(${PROJECT_NAME}
    ${PCL_LIBRARIES}
    Eigen3::Eigen
    ${Boost_LIBRARIES}
    OpenMP::OpenMP_CXX
    ${OPENGL_LIBRARIES}
)

# Windows specific libraries
if(WIN32)
    target_link_libraries(${PROJECT_NAME}
        user32
        shell32
        gdi32
        opengl32
    )
endif()

# Compiler specific options
if(MSVC)
    # Add /bigobj flag for large object files
    target_compile_options(${PROJECT_NAME} PRIVATE /bigobj)
    
    # Set warning level
    target_compile_options(${PROJECT_NAME} PRIVATE /W3)
    
    # Enable OpenMP support
    target_compile_options(${PROJECT_NAME} PRIVATE /openmp)
    
    # Release optimizations
    target_compile_options(${PROJECT_NAME} PRIVATE 
        $<$<CONFIG:Release>:/O2 /Oi /GL>
    )
    
    # Release linker optimizations
    target_link_options(${PROJECT_NAME} PRIVATE 
        $<$<CONFIG:Release>:/LTCG /OPT:REF /OPT:ICF>
    )
endif()

# Set output directory
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/bin/Debug
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin/Release
)


# Print build information
message(STATUS "=== DACPO Build Configuration ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "PCL version: ${PCL_VERSION}")
message(STATUS "Eigen3 version: ${Eigen3_VERSION}")
message(STATUS "Boost version: ${Boost_VERSION}")
message(STATUS "OpenMP found: ${OpenMP_FOUND}")
message(STATUS "=================================") 