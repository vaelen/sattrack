# Test configuration
find_package(GTest REQUIRED)

# Include directories for tests
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${GTEST_INCLUDE_DIRS}
)

# Test source files (shared library sources, excluding main)
set(LIB_SOURCES
    ${CMAKE_SOURCE_DIR}/src/satellite.cpp
    ${CMAKE_SOURCE_DIR}/src/sgp4.cpp
)

# Satellite tests
add_executable(satellite_test
    satellite_test.cpp
    ${LIB_SOURCES}
)

target_link_libraries(satellite_test
    GTest::gtest
    GTest::gtest_main
)

# Celestrak tests
add_executable(celestrak_test
    celestrak_test.cpp
    ${CMAKE_SOURCE_DIR}/src/celestrak.cpp
)

target_link_libraries(celestrak_test
    GTest::gtest
    GTest::gtest_main
    ${CURL_LIBRARIES}
    ${CURLPP_LIBRARIES}
)

# GPS tests
add_executable(gps_test
    gps_test.cpp
    ${CMAKE_SOURCE_DIR}/src/gps.cpp
    ${LIB_SOURCES}
)

target_link_libraries(gps_test
    GTest::gtest
    GTest::gtest_main
    spdlog::spdlog
)

# Register tests with CTest
include(GoogleTest)
gtest_discover_tests(satellite_test)
gtest_discover_tests(celestrak_test)
gtest_discover_tests(gps_test)
