/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <gtest/gtest.h>
#include <sattrack/orbit.hpp>

#include <chrono>
#include <cmath>
#include <string>

namespace sattrack {
namespace {

// ISS TLE data from Celestrak (real example)
constexpr const char* ISS_TLE =
    "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
    "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850";

// NOAA 19 TLE from Celestrak
constexpr const char* NOAA19_TLE =
    "1 33591U 09005A   25333.78204194  .00000054  00000+0  52635-4 0  9999\n"
    "2 33591  98.9785  39.2910 0013037 231.6546 128.3455 14.13431889866318";

// TLE with name line (three-line format)
constexpr const char* TLE_WITH_NAME =
    "ISS (ZARYA)\n"
    "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
    "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850";

class OrbitTest : public ::testing::Test {
protected:
    Orbit orbit;
};

// Test default construction
TEST_F(OrbitTest, DefaultConstruction) {
    Orbit defaultOrbit;
    // Default constructed orbit has uninitialized numeric members
    // but string should be empty (std::string default constructor)
    EXPECT_TRUE(defaultOrbit.getDesignator().empty());
}

// Test parsing ISS TLE - Line 1 fields
TEST_F(OrbitTest, ParseISSLine1_NoradID) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_EQ(orbit.getNoradID(), 25544);
}

TEST_F(OrbitTest, ParseISSLine1_Classification) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_EQ(orbit.getClassification(), 'U');
}

TEST_F(OrbitTest, ParseISSLine1_Designator) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_EQ(orbit.getDesignator(), "98067A  ");
}

TEST_F(OrbitTest, ParseISSLine1_FirstDerivativeMeanMotion) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_NEAR(orbit.getFirstDerivativeMeanMotion(), 0.00008010, 1e-8);
}

TEST_F(OrbitTest, ParseISSLine1_SecondDerivativeMeanMotion) {
    orbit.updateFromTLE(ISS_TLE);
    // "00000+0" should parse to 0.0
    EXPECT_NEAR(orbit.getSecondDerivativeMeanMotion(), 0.0, 1e-10);
}

TEST_F(OrbitTest, ParseISSLine1_BstarDragTerm) {
    orbit.updateFromTLE(ISS_TLE);
    // "15237-3" -> 0.00015237
    EXPECT_NEAR(orbit.getBstarDragTerm(), 0.00015237, 1e-8);
}

TEST_F(OrbitTest, ParseISSLine1_ElementSetNumber) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_EQ(orbit.getElementSetNumber(), 999);
}

// Test parsing ISS TLE - Line 2 fields
TEST_F(OrbitTest, ParseISSLine2_Inclination) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_DOUBLE_EQ(orbit.getInclination(), 51.6312);
}

TEST_F(OrbitTest, ParseISSLine2_RAAN) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_DOUBLE_EQ(orbit.getRAAN(), 206.3646);
}

TEST_F(OrbitTest, ParseISSLine2_Eccentricity) {
    orbit.updateFromTLE(ISS_TLE);
    // "0003723" -> 0.0003723
    EXPECT_NEAR(orbit.getEccentricity(), 0.0003723, 1e-8);
}

TEST_F(OrbitTest, ParseISSLine2_ArgumentOfPerigee) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_DOUBLE_EQ(orbit.getArgumentOfPerigee(), 184.1118);
}

TEST_F(OrbitTest, ParseISSLine2_MeanAnomaly) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_DOUBLE_EQ(orbit.getMeanAnomaly(), 175.9840);
}

TEST_F(OrbitTest, ParseISSLine2_MeanMotion) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_NEAR(orbit.getMeanMotion(), 15.49193835, 1e-6);
}

TEST_F(OrbitTest, ParseISSLine2_RevolutionNumber) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_EQ(orbit.getRevolutionNumberAtEpoch(), 54085);
}

// Test parsing NOAA-19 TLE for different values
TEST_F(OrbitTest, ParseNOAA19_NoradID) {
    orbit.updateFromTLE(NOAA19_TLE);
    EXPECT_EQ(orbit.getNoradID(), 33591);
}

TEST_F(OrbitTest, ParseNOAA19_Inclination) {
    orbit.updateFromTLE(NOAA19_TLE);
    EXPECT_DOUBLE_EQ(orbit.getInclination(), 98.9785);
}

TEST_F(OrbitTest, ParseNOAA19_BstarDragTerm) {
    orbit.updateFromTLE(NOAA19_TLE);
    // "52635-4" -> 0.000052635
    EXPECT_NEAR(orbit.getBstarDragTerm(), 0.000052635, 1e-9);
}

TEST_F(OrbitTest, ParseNOAA19_MeanMotion) {
    orbit.updateFromTLE(NOAA19_TLE);
    EXPECT_NEAR(orbit.getMeanMotion(), 14.13431889, 1e-6);
}

// Test three-line format (with satellite name)
TEST_F(OrbitTest, ParseThreeLineFormat) {
    orbit.updateFromTLE(TLE_WITH_NAME);
    EXPECT_EQ(orbit.getNoradID(), 25544);
    EXPECT_DOUBLE_EQ(orbit.getInclination(), 51.6312);
}

// Test epoch parsing
TEST_F(OrbitTest, EpochParsing) {
    orbit.updateFromTLE(ISS_TLE);
    auto epoch = orbit.getEpoch();

    // Convert to time_t for easier verification
    auto time_t_epoch = std::chrono::system_clock::to_time_t(epoch);
    std::tm* tm = std::gmtime(&time_t_epoch);

    // Epoch "25333.83453771" = Day 333 of 2025 = November 29, 2025
    EXPECT_EQ(tm->tm_year + 1900, 2025);
    EXPECT_EQ(tm->tm_mon + 1, 11);  // November
    EXPECT_EQ(tm->tm_mday, 29);
}

// Test updating orbit with new TLE
TEST_F(OrbitTest, UpdateWithNewTLE) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_EQ(orbit.getNoradID(), 25544);

    orbit.updateFromTLE(NOAA19_TLE);
    EXPECT_EQ(orbit.getNoradID(), 33591);
}

// Test orbital element ranges (sanity checks)
TEST_F(OrbitTest, InclinationRange) {
    orbit.updateFromTLE(ISS_TLE);
    double incl = orbit.getInclination();
    EXPECT_GE(incl, 0.0);
    EXPECT_LE(incl, 180.0);
}

TEST_F(OrbitTest, EccentricityRange) {
    orbit.updateFromTLE(ISS_TLE);
    double ecc = orbit.getEccentricity();
    EXPECT_GE(ecc, 0.0);
    EXPECT_LT(ecc, 1.0);  // Must be less than 1 for elliptical orbit
}

TEST_F(OrbitTest, RAANRange) {
    orbit.updateFromTLE(ISS_TLE);
    double raan = orbit.getRAAN();
    EXPECT_GE(raan, 0.0);
    EXPECT_LT(raan, 360.0);
}

TEST_F(OrbitTest, ArgumentOfPerigeeRange) {
    orbit.updateFromTLE(ISS_TLE);
    double aop = orbit.getArgumentOfPerigee();
    EXPECT_GE(aop, 0.0);
    EXPECT_LT(aop, 360.0);
}

TEST_F(OrbitTest, MeanAnomalyRange) {
    orbit.updateFromTLE(ISS_TLE);
    double ma = orbit.getMeanAnomaly();
    EXPECT_GE(ma, 0.0);
    EXPECT_LT(ma, 360.0);
}

TEST_F(OrbitTest, MeanMotionPositive) {
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_GT(orbit.getMeanMotion(), 0.0);
}

// Test physical interpretation of orbital elements
TEST_F(OrbitTest, ISSIsLowEarthOrbit) {
    orbit.updateFromTLE(ISS_TLE);
    // ISS has mean motion ~15.5 rev/day, indicating LEO
    // LEO satellites typically have mean motion > 11 rev/day
    EXPECT_GT(orbit.getMeanMotion(), 11.0);
}

TEST_F(OrbitTest, ISSHasLowEccentricity) {
    orbit.updateFromTLE(ISS_TLE);
    // ISS has nearly circular orbit
    EXPECT_LT(orbit.getEccentricity(), 0.01);
}

TEST_F(OrbitTest, NOAA19IsSunSynchronous) {
    orbit.updateFromTLE(NOAA19_TLE);
    // Sun-synchronous orbits have inclination ~98-99 degrees
    double incl = orbit.getInclination();
    EXPECT_GT(incl, 97.0);
    EXPECT_LT(incl, 100.0);
}

}  // namespace
}  // namespace sattrack
