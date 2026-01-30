cmake_minimum_required(VERSION 3.14)

project(SatTrack VERSION 1.0.0 LANGUAGES CXX C)

# =============================================================================
# Cross-compilation options
# =============================================================================
option(SATTRACK_STATIC_LINK "Prefer static linking for dependencies" OFF)

# When cross-compiling, default to static linking
if(CMAKE_CROSSCOMPILING)
    set(SATTRACK_STATIC_LINK ON CACHE BOOL "Prefer static linking" FORCE)
    message(STATUS "Cross-compilation detected - enabling static linking")
endif()

# =============================================================================
# C++ Standard
# =============================================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# =============================================================================
# FetchContent for header-only and easily-built libraries
# =============================================================================
include(FetchContent)

# RapidJSON (header-only) - only fetch, don't process its CMakeLists.txt
FetchContent_Declare(
    rapidjson
    GIT_REPOSITORY https://github.com/Tencent/rapidjson.git
    GIT_TAG        v1.1.0
    GIT_SHALLOW    TRUE
)
FetchContent_GetProperties(rapidjson)
if(NOT rapidjson_POPULATED)
    FetchContent_Populate(rapidjson)
endif()
set(RAPIDJSON_INCLUDE_DIRS ${rapidjson_SOURCE_DIR}/include)

# CLI11 (header-only)
FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG        v2.4.1
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(cli11)

# Asio (header-only, standalone version)
FetchContent_Declare(
    asio
    GIT_REPOSITORY https://github.com/chriskohlhoff/asio.git
    GIT_TAG        asio-1-30-2
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(asio)
add_library(asio INTERFACE)
target_include_directories(asio INTERFACE ${asio_SOURCE_DIR}/asio/include)
target_compile_definitions(asio INTERFACE ASIO_STANDALONE)

# i2clcd - LCD display control via I2C (Makefile project, compile manually)
FetchContent_Declare(
    i2clcd
    GIT_REPOSITORY https://github.com/vaelen/i2clcd.git
    GIT_TAG        main
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(i2clcd)
add_library(i2clcd STATIC ${i2clcd_SOURCE_DIR}/src/i2clcd.c)
target_include_directories(i2clcd PUBLIC
    ${i2clcd_SOURCE_DIR}/include
    ${i2clcd_SOURCE_DIR}/src
)

# spdlog - use FetchContent for static linking, system package otherwise
if(SATTRACK_STATIC_LINK)
    # Build spdlog with bundled fmt for static linking
    set(SPDLOG_FMT_EXTERNAL OFF CACHE BOOL "" FORCE)
    set(SPDLOG_BUILD_SHARED OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v1.13.0
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(spdlog)
else()
    # Use system packages for native builds
    find_package(spdlog REQUIRED)
    find_package(fmt REQUIRED)
endif()

# =============================================================================
# System Dependencies
# =============================================================================
find_package(PkgConfig REQUIRED)

# libcurl
if(SATTRACK_STATIC_LINK)
    set(CURL_USE_STATIC_LIBS ON)
endif()
find_package(CURL REQUIRED)

# curlpp (manual detection - no standard CMake config)
find_path(CURLPP_INCLUDE_DIRS "curlpp/cURLpp.hpp"
    PATHS
        ${CMAKE_SYSROOT}/usr/include
        /usr/include
        /usr/local/include
)

if(SATTRACK_STATIC_LINK)
    find_library(CURLPP_LIBRARIES
        NAMES libcurlpp.a curlpp
        PATHS
            ${CMAKE_SYSROOT}/usr/lib
            /usr/lib
            /usr/local/lib
        NO_DEFAULT_PATH
    )
    # Fallback to default paths if not found in sysroot
    if(NOT CURLPP_LIBRARIES)
        find_library(CURLPP_LIBRARIES
            NAMES libcurlpp.a curlpp
        )
    endif()
else()
    find_library(CURLPP_LIBRARIES
        NAMES curlpp
        PATHS
            /usr/lib
            /usr/local/lib
    )
endif()

if(NOT CURLPP_INCLUDE_DIRS OR NOT CURLPP_LIBRARIES)
    message(FATAL_ERROR "curlpp not found. Please install libcurlpp-dev or build the ARM sysroot")
endif()

# =============================================================================
# Handle static linking dependencies
# =============================================================================
if(SATTRACK_STATIC_LINK)
    # For fully static curl, we need to link its dependencies
    find_library(OPENSSL_SSL_LIBRARY
        NAMES libssl.a ssl
        PATHS ${CMAKE_SYSROOT}/usr/lib /usr/lib
        NO_DEFAULT_PATH
    )
    find_library(OPENSSL_CRYPTO_LIBRARY
        NAMES libcrypto.a crypto
        PATHS ${CMAKE_SYSROOT}/usr/lib /usr/lib
        NO_DEFAULT_PATH
    )
    find_library(ZLIB_LIBRARY_STATIC
        NAMES libz.a z
        PATHS ${CMAKE_SYSROOT}/usr/lib /usr/lib
        NO_DEFAULT_PATH
    )

    # Collect static dependencies for curl
    set(CURL_STATIC_DEPS "")
    if(OPENSSL_SSL_LIBRARY AND OPENSSL_CRYPTO_LIBRARY)
        list(APPEND CURL_STATIC_DEPS ${OPENSSL_SSL_LIBRARY} ${OPENSSL_CRYPTO_LIBRARY})
    endif()
    if(ZLIB_LIBRARY_STATIC)
        list(APPEND CURL_STATIC_DEPS ${ZLIB_LIBRARY_STATIC})
    endif()
    # pthread, dl are needed for OpenSSL; atomic is needed for 64-bit atomics on 32-bit ARM
    list(APPEND CURL_STATIC_DEPS pthread dl atomic)
endif()

# =============================================================================
# Include directories
# =============================================================================
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CURL_INCLUDE_DIRS}
    ${RAPIDJSON_INCLUDE_DIRS}
    ${CURLPP_INCLUDE_DIRS}
)

# =============================================================================
# Source files
# =============================================================================
set(SOURCES
    src/main.cpp
    src/celestrak.cpp
    src/config.cpp
    src/satellite.cpp
    src/sgp4.cpp
    src/gps.cpp
    src/daemon.cpp
    src/serialport.cpp
    src/rotator.cpp
    src/lcd.cpp
)

# =============================================================================
# Create executable
# =============================================================================
add_executable(sattrack ${SOURCES})

# Link libraries
if(SATTRACK_STATIC_LINK)
    target_link_libraries(sattrack
        ${CURLPP_LIBRARIES}
        ${CURL_LIBRARIES}
        ${CURL_STATIC_DEPS}
        spdlog::spdlog
        asio
        CLI11::CLI11
        i2clcd
    )
else()
    target_link_libraries(sattrack
        ${CURL_LIBRARIES}
        ${CURLPP_LIBRARIES}
        spdlog::spdlog
        fmt::fmt
        asio
        CLI11::CLI11
        i2clcd
    )
endif()

# =============================================================================
# Install targets
# =============================================================================
install(TARGETS sattrack DESTINATION bin)

# =============================================================================
# Testing (disabled for cross-compilation)
# =============================================================================
if(NOT CMAKE_CROSSCOMPILING)
    find_package(GTest QUIET)
    if(GTest_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
        enable_testing()
        add_subdirectory(tests)
        message(STATUS "  GTest: Found - tests enabled")
    else()
        if(NOT GTest_FOUND)
            message(STATUS "  GTest: Not found - tests disabled")
        endif()
    endif()
else()
    message(STATUS "  Tests: Disabled (cross-compiling)")
endif()

# =============================================================================
# Configuration summary
# =============================================================================
message(STATUS "")
message(STATUS "SatTrack Configuration Summary:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Cross-compiling: ${CMAKE_CROSSCOMPILING}")
message(STATUS "  Static linking: ${SATTRACK_STATIC_LINK}")
if(CMAKE_CROSSCOMPILING)
    message(STATUS "  Target system: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
    message(STATUS "  Sysroot: ${CMAKE_SYSROOT}")
endif()
message(STATUS "  Install Prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
message(STATUS "Dependencies:")
message(STATUS "  libcurl: ${CURL_VERSION_STRING}")
message(STATUS "  RapidJSON: FetchContent (v1.1.0)")
message(STATUS "  curlpp: ${CURLPP_LIBRARIES}")
message(STATUS "  CLI11: FetchContent (v2.4.1)")
if(SATTRACK_STATIC_LINK)
    message(STATUS "  spdlog: FetchContent (v1.13.0, static)")
else()
    message(STATUS "  spdlog: ${spdlog_VERSION}")
    message(STATUS "  fmt: ${fmt_VERSION}")
endif()
message(STATUS "  Asio: FetchContent (v1.30.2)")
message(STATUS "  i2clcd: FetchContent (main)")
message(STATUS "")
