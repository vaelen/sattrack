cmake_minimum_required(VERSION 3.10)
project(SatTrack VERSION 1.0.0 LANGUAGES CXX)

# Set C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Find required packages
find_package(PkgConfig REQUIRED)

# Find GTKmm
pkg_check_modules(GTKMM REQUIRED gtkmm-3.0)

# Find libcurl
find_package(CURL REQUIRED)

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

# TOML++ (header-only library)
find_path(TOMLPLUSPLUS_INCLUDE_DIR "toml++/toml.hpp"
    PATHS
        /usr/include
        /usr/local/include
)

if(NOT TOMLPLUSPLUS_INCLUDE_DIR)
    message(FATAL_ERROR "TOML++ not found. Please install libtomlplusplus-dev")
endif()

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${GTKMM_INCLUDE_DIRS}
    ${CURL_INCLUDE_DIRS}
    ${RAPIDJSON_INCLUDE_DIRS}
    ${CURLPP_INCLUDE_DIRS}
    ${CLI11_INCLUDE_DIR}
    ${TOMLPLUSPLUS_INCLUDE_DIR}
)

# Link directories
link_directories(
    ${GTKMM_LIBRARY_DIRS}
)

# Source files
set(SOURCES
    src/main.cpp
    src/gui.cpp
    src/n2yo.cpp
    # Add more source files here as you create them
)

# Create executable
add_executable(sattrack ${SOURCES})

# Link libraries
target_link_libraries(sattrack
    ${GTKMM_LIBRARIES}
    ${CURL_LIBRARIES}
    ${CURLPP_LIBRARIES}
)

# Compiler options
target_compile_options(sattrack PRIVATE
    ${GTKMM_CFLAGS_OTHER}
)

# Install targets
install(TARGETS sattrack DESTINATION bin)

# Optional: Enable testing
if(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
    enable_testing()
    add_subdirectory(tests)
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
message(STATUS "  GTKmm: ${GTKMM_VERSION}")
message(STATUS "  libcurl: ${CURL_VERSION_STRING}")
message(STATUS "  RapidJSON: ${RAPIDJSON_INCLUDE_DIRS}")
message(STATUS "  curlpp: ${CURLPP_LIBRARIES}")
message(STATUS "  CLI11: ${CLI11_INCLUDE_DIR}")
message(STATUS "  TOML++: ${TOMLPLUSPLUS_INCLUDE_DIR}")
message(STATUS "")
