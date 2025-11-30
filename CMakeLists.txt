cmake_minimum_required(VERSION 3.10)
project(SatTrack VERSION 1.0.0 LANGUAGES CXX)

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Find required packages
find_package(PkgConfig REQUIRED)

# Find libcurl
find_package(CURL REQUIRED)

# Find GTest (optional, for testing)
find_package(GTest QUIET)

# RapidJSON (header-only library)
find_path(RAPIDJSON_INCLUDE_DIRS "rapidjson/rapidjson.h"
    PATHS
        /usr/include
        /usr/local/include
    PATH_SUFFIXES rapidjson
)

if(NOT RAPIDJSON_INCLUDE_DIRS)
    message(FATAL_ERROR "RapidJSON not found. Please install rapidjson-dev")
endif()

# curlpp
find_path(CURLPP_INCLUDE_DIRS "curlpp/cURLpp.hpp"
    PATHS
        /usr/include
        /usr/local/include
)

find_library(CURLPP_LIBRARIES
    NAMES curlpp
    PATHS
        /usr/lib
        /usr/local/lib
)

if(NOT CURLPP_INCLUDE_DIRS OR NOT CURLPP_LIBRARIES)
    message(FATAL_ERROR "curlpp not found. Please install libcurlpp-dev")
endif()

# CLI11 (header-only library)
find_path(CLI11_INCLUDE_DIR "CLI/CLI.hpp"
    PATHS
        /usr/include
        /usr/local/include
)

if(NOT CLI11_INCLUDE_DIR)
    message(FATAL_ERROR "CLI11 not found. Please install libcli11-dev")
endif()

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CURL_INCLUDE_DIRS}
    ${RAPIDJSON_INCLUDE_DIRS}
    ${CURLPP_INCLUDE_DIRS}
    ${CLI11_INCLUDE_DIR}
)

# Source files
set(SOURCES
    src/main.cpp
    src/n2yo.cpp
    src/config.cpp
    src/orbit.cpp
    # Add more source files here as you create them
)

# Create executable
add_executable(sattrack ${SOURCES})

# Link libraries
target_link_libraries(sattrack
    ${CURL_LIBRARIES}
    ${CURLPP_LIBRARIES}
)

# Install targets
install(TARGETS sattrack DESTINATION bin)

# Optional: Enable testing
if(GTest_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
    enable_testing()
    add_subdirectory(tests)
    message(STATUS "  GTest: Found - tests enabled")
else()
    if(NOT GTest_FOUND)
        message(STATUS "  GTest: Not found - tests disabled")
    endif()
endif()

# Print configuration summary
message(STATUS "")
message(STATUS "SatTrack Configuration Summary:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Install Prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
message(STATUS "Dependencies:")
message(STATUS "  libcurl: ${CURL_VERSION_STRING}")
message(STATUS "  RapidJSON: ${RAPIDJSON_INCLUDE_DIRS}")
message(STATUS "  curlpp: ${CURLPP_LIBRARIES}")
message(STATUS "  CLI11: ${CLI11_INCLUDE_DIR}")
message(STATUS "")
