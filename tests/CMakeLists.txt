# Test configuration
find_package(GTest REQUIRED)

# Include directories for tests
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${GTEST_INCLUDE_DIRS}
)

# Test source files (shared library sources, excluding main)
set(LIB_SOURCES
    ${CMAKE_SOURCE_DIR}/src/orbit.cpp
)

# Orbit tests
add_executable(orbit_test
    orbit_test.cpp
    ${LIB_SOURCES}
)

target_link_libraries(orbit_test
    GTest::gtest
    GTest::gtest_main
)

# Register tests with CTest
include(GoogleTest)
gtest_discover_tests(orbit_test)
