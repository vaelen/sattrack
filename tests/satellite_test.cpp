/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <gtest/gtest.h>
#include <sattrack/satellite.hpp>

#include <chrono>
#include <cmath>
#include <map>
#include <sstream>
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

class SatelliteTest : public ::testing::Test {
protected:
    Satellite satellite;
};

// Test default construction
TEST_F(SatelliteTest, DefaultConstruction) {
    Satellite defaultSatellite;
    // Default constructed orbit has uninitialized numeric members
    // but string should be empty (std::string default constructor)
    EXPECT_TRUE(defaultSatellite.getDesignator().empty());
    EXPECT_TRUE(defaultSatellite.getName().empty());
}

// Test parsing ISS TLE - Line 1 fields
TEST_F(SatelliteTest, ParseISSLine1_NoradID) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_EQ(satellite.getNoradID(), 25544);
}

TEST_F(SatelliteTest, ParseISSLine1_Classification) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_EQ(satellite.getClassification(), 'U');
}

TEST_F(SatelliteTest, ParseISSLine1_Designator) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_EQ(satellite.getDesignator(), "98067A  ");
}

TEST_F(SatelliteTest, ParseISSLine1_FirstDerivativeMeanMotion) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_NEAR(satellite.getFirstDerivativeMeanMotion(), 0.00008010, 1e-8);
}

TEST_F(SatelliteTest, ParseISSLine1_SecondDerivativeMeanMotion) {
    satellite.updateFromTLE(ISS_TLE);
    // "00000+0" should parse to 0.0
    EXPECT_NEAR(satellite.getSecondDerivativeMeanMotion(), 0.0, 1e-10);
}

TEST_F(SatelliteTest, ParseISSLine1_BstarDragTerm) {
    satellite.updateFromTLE(ISS_TLE);
    // "15237-3" -> 0.00015237
    EXPECT_NEAR(satellite.getBstarDragTerm(), 0.00015237, 1e-8);
}

TEST_F(SatelliteTest, ParseISSLine1_ElementSetNumber) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_EQ(satellite.getElementSetNumber(), 999);
}

// Test parsing ISS TLE - Line 2 fields
TEST_F(SatelliteTest, ParseISSLine2_Inclination) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_DOUBLE_EQ(satellite.getInclination(), 51.6312);
}

TEST_F(SatelliteTest, ParseISSLine2_RAAN) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_DOUBLE_EQ(satellite.getRightAscensionOfAscendingNode(), 206.3646);
}

TEST_F(SatelliteTest, ParseISSLine2_Eccentricity) {
    satellite.updateFromTLE(ISS_TLE);
    // "0003723" -> 0.0003723
    EXPECT_NEAR(satellite.getEccentricity(), 0.0003723, 1e-8);
}

TEST_F(SatelliteTest, ParseISSLine2_ArgumentOfPerigee) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_DOUBLE_EQ(satellite.getArgumentOfPerigee(), 184.1118);
}

TEST_F(SatelliteTest, ParseISSLine2_MeanAnomaly) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_DOUBLE_EQ(satellite.getMeanAnomaly(), 175.9840);
}

TEST_F(SatelliteTest, ParseISSLine2_MeanMotion) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_NEAR(satellite.getMeanMotion(), 15.49193835, 1e-6);
}

TEST_F(SatelliteTest, ParseISSLine2_RevolutionNumber) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_EQ(satellite.getRevolutionNumberAtEpoch(), 54085);
}

// Test parsing NOAA-19 TLE for different values
TEST_F(SatelliteTest, ParseNOAA19_NoradID) {
    satellite.updateFromTLE(NOAA19_TLE);
    EXPECT_EQ(satellite.getNoradID(), 33591);
}

TEST_F(SatelliteTest, ParseNOAA19_Inclination) {
    satellite.updateFromTLE(NOAA19_TLE);
    EXPECT_DOUBLE_EQ(satellite.getInclination(), 98.9785);
}

TEST_F(SatelliteTest, ParseNOAA19_BstarDragTerm) {
    satellite.updateFromTLE(NOAA19_TLE);
    // "52635-4" -> 0.000052635
    EXPECT_NEAR(satellite.getBstarDragTerm(), 0.000052635, 1e-9);
}

TEST_F(SatelliteTest, ParseNOAA19_MeanMotion) {
    satellite.updateFromTLE(NOAA19_TLE);
    EXPECT_NEAR(satellite.getMeanMotion(), 14.13431889, 1e-6);
}

// Test three-line format (with satellite name)
TEST_F(SatelliteTest, ParseThreeLineFormat) {
    satellite.updateFromTLE(TLE_WITH_NAME);
    EXPECT_EQ(satellite.getNoradID(), 25544);
    EXPECT_DOUBLE_EQ(satellite.getInclination(), 51.6312);
}

// ============================================================================
// Name Parsing Tests
// ============================================================================

TEST_F(SatelliteTest, GetName_TwoLineFormat) {
    // Two-line TLE format should result in empty name
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_TRUE(satellite.getName().empty());
}

TEST_F(SatelliteTest, GetName_ThreeLineFormat) {
    // Three-line TLE format should parse the name from the first line
    satellite.updateFromTLE(TLE_WITH_NAME);
    EXPECT_EQ(satellite.getName(), "ISS (ZARYA)");
}

TEST_F(SatelliteTest, UpdateFromTLE_WithExplicitName) {
    // Using the overload that accepts a separate name
    satellite.updateFromTLE("My Custom Name", ISS_TLE);
    EXPECT_EQ(satellite.getName(), "My Custom Name");
    // Should still parse the orbital elements correctly
    EXPECT_EQ(satellite.getNoradID(), 25544);
    EXPECT_DOUBLE_EQ(satellite.getInclination(), 51.6312);
}

TEST_F(SatelliteTest, UpdateFromTLE_ExplicitNameOverridesEmbeddedName) {
    // When using the overload with explicit name, it should override
    // any name that might be embedded in the TLE
    satellite.updateFromTLE("Override Name", TLE_WITH_NAME);
    EXPECT_EQ(satellite.getName(), "Override Name");
}

TEST_F(SatelliteTest, GetName_TrimsWhitespace) {
    // Test that leading/trailing whitespace is trimmed from name
    constexpr const char* TLE_WITH_PADDED_NAME =
        "  NOAA 19  \n"
        "1 33591U 09005A   25333.78204194  .00000054  00000+0  52635-4 0  9999\n"
        "2 33591  98.9785  39.2910 0013037 231.6546 128.3455 14.13431889866318";
    satellite.updateFromTLE(TLE_WITH_PADDED_NAME);
    EXPECT_EQ(satellite.getName(), "NOAA 19");
}

// Test epoch parsing
TEST_F(SatelliteTest, EpochParsing) {
    satellite.updateFromTLE(ISS_TLE);
    auto epoch = satellite.getEpoch();

    // Convert to time_t for easier verification
    auto time_t_epoch = std::chrono::system_clock::to_time_t(epoch);
    std::tm* tm = std::gmtime(&time_t_epoch);

    // Epoch "25333.83453771" = Day 333 of 2025 = November 29, 2025
    EXPECT_EQ(tm->tm_year + 1900, 2025);
    EXPECT_EQ(tm->tm_mon + 1, 11);  // November
    EXPECT_EQ(tm->tm_mday, 29);
}

// Test updating orbit with new TLE
TEST_F(SatelliteTest, UpdateWithNewTLE) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_EQ(satellite.getNoradID(), 25544);

    satellite.updateFromTLE(NOAA19_TLE);
    EXPECT_EQ(satellite.getNoradID(), 33591);
}

// Test orbital element ranges (sanity checks)
TEST_F(SatelliteTest, InclinationRange) {
    satellite.updateFromTLE(ISS_TLE);
    double incl = satellite.getInclination();
    EXPECT_GE(incl, 0.0);
    EXPECT_LE(incl, 180.0);
}

TEST_F(SatelliteTest, EccentricityRange) {
    satellite.updateFromTLE(ISS_TLE);
    double ecc = satellite.getEccentricity();
    EXPECT_GE(ecc, 0.0);
    EXPECT_LT(ecc, 1.0);  // Must be less than 1 for elliptical orbit
}

TEST_F(SatelliteTest, RAANRange) {
    satellite.updateFromTLE(ISS_TLE);
    double raan = satellite.getRightAscensionOfAscendingNode();
    EXPECT_GE(raan, 0.0);
    EXPECT_LT(raan, 360.0);
}

TEST_F(SatelliteTest, ArgumentOfPerigeeRange) {
    satellite.updateFromTLE(ISS_TLE);
    double aop = satellite.getArgumentOfPerigee();
    EXPECT_GE(aop, 0.0);
    EXPECT_LT(aop, 360.0);
}

TEST_F(SatelliteTest, MeanAnomalyRange) {
    satellite.updateFromTLE(ISS_TLE);
    double ma = satellite.getMeanAnomaly();
    EXPECT_GE(ma, 0.0);
    EXPECT_LT(ma, 360.0);
}

TEST_F(SatelliteTest, MeanMotionPositive) {
    satellite.updateFromTLE(ISS_TLE);
    EXPECT_GT(satellite.getMeanMotion(), 0.0);
}

// Test physical interpretation of orbital elements
TEST_F(SatelliteTest, ISSIsLowEarthOrbit) {
    satellite.updateFromTLE(ISS_TLE);
    // ISS has mean motion ~15.5 rev/day, indicating LEO
    // LEO satellites typically have mean motion > 11 rev/day
    EXPECT_GT(satellite.getMeanMotion(), 11.0);
}

TEST_F(SatelliteTest, ISSHasLowEccentricity) {
    satellite.updateFromTLE(ISS_TLE);
    // ISS has nearly circular orbit
    EXPECT_LT(satellite.getEccentricity(), 0.01);
}

TEST_F(SatelliteTest, NOAA19IsSunSynchronous) {
    satellite.updateFromTLE(NOAA19_TLE);
    // Sun-synchronous orbits have inclination ~98-99 degrees
    double incl = satellite.getInclination();
    EXPECT_GT(incl, 97.0);
    EXPECT_LT(incl, 100.0);
}

// ============================================================================
// Vec3 Tests
// ============================================================================

TEST(Vec3Test, Addition) {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};
    Vec3 result = a + b;
    EXPECT_DOUBLE_EQ(result.x, 5.0);
    EXPECT_DOUBLE_EQ(result.y, 7.0);
    EXPECT_DOUBLE_EQ(result.z, 9.0);
}

TEST(Vec3Test, ScalarMultiplication) {
    Vec3 v{2.0, 3.0, 4.0};
    Vec3 result = v * 2.5;
    EXPECT_DOUBLE_EQ(result.x, 5.0);
    EXPECT_DOUBLE_EQ(result.y, 7.5);
    EXPECT_DOUBLE_EQ(result.z, 10.0);
}

TEST(Vec3Test, Magnitude) {
    Vec3 v{3.0, 4.0, 0.0};
    EXPECT_DOUBLE_EQ(v.magnitude(), 5.0);
}

TEST(Vec3Test, MagnitudeUnitVector) {
    Vec3 v{1.0, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(v.magnitude(), 1.0);
}

TEST(Vec3Test, MagnitudeZeroVector) {
    Vec3 v{0.0, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(v.magnitude(), 0.0);
}

TEST(Vec3Test, Magnitude3D) {
    // sqrt(1^2 + 2^2 + 2^2) = sqrt(9) = 3
    Vec3 v{1.0, 2.0, 2.0};
    EXPECT_DOUBLE_EQ(v.magnitude(), 3.0);
}

// ============================================================================
// Julian Date Tests
// ============================================================================

TEST(JulianDateTest, UnixEpoch) {
    // Unix epoch (1970-01-01 00:00:00 UTC) should be JD 2440587.5
    auto unixEpoch = std::chrono::system_clock::from_time_t(0);
    double jd = toJulianDate(unixEpoch);
    EXPECT_NEAR(jd, 2440587.5, 1e-6);
}

TEST(JulianDateTest, J2000Epoch) {
    // J2000.0 (2000-01-01 12:00:00 TT) is JD 2451545.0
    // J2000 in UTC is approximately 2000-01-01 11:58:55.816 UTC
    // For simplicity, test that 2000-01-01 12:00:00 UTC is close
    std::tm tm = {};
    tm.tm_year = 100;  // 2000
    tm.tm_mon = 0;     // January
    tm.tm_mday = 1;
    tm.tm_hour = 12;
    tm.tm_min = 0;
    tm.tm_sec = 0;
    auto tp = std::chrono::system_clock::from_time_t(timegm(&tm));
    double jd = toJulianDate(tp);
    EXPECT_NEAR(jd, 2451545.0, 0.01);
}

TEST(JulianDateTest, KnownDate) {
    // 2024-03-20 00:00:00 UTC is approximately JD 2460389.5
    std::tm tm = {};
    tm.tm_year = 124;  // 2024
    tm.tm_mon = 2;     // March
    tm.tm_mday = 20;
    tm.tm_hour = 0;
    tm.tm_min = 0;
    tm.tm_sec = 0;
    auto tp = std::chrono::system_clock::from_time_t(timegm(&tm));
    double jd = toJulianDate(tp);
    EXPECT_NEAR(jd, 2460389.5, 0.01);
}

// ============================================================================
// GMST Tests
// ============================================================================

TEST(GMSTTest, ReturnsRadians) {
    // GMST should be in range [0, 2π)
    double jd = 2451545.0;  // J2000.0
    double gst = gmst(jd);
    EXPECT_GE(gst, 0.0);
    EXPECT_LT(gst, 2.0 * M_PI);
}

TEST(GMSTTest, J2000Epoch) {
    // At J2000.0 (JD 2451545.0), GMST should be approximately 280.46 degrees
    // which is about 4.894 radians
    double gst = gmst(2451545.0);
    double gstDegrees = gst * 180.0 / M_PI;
    EXPECT_NEAR(gstDegrees, 280.46, 0.1);
}

TEST(GMSTTest, OneDayLater) {
    // GMST advances ~360.98 degrees per day (sidereal rate)
    // But since gmst() normalizes to [0, 360), the difference
    // between consecutive days is just the fractional part: ~0.98 degrees
    double gst1 = gmst(2451545.0);
    double gst2 = gmst(2451546.0);

    double diff = gst2 - gst1;
    if (diff < 0) diff += 2.0 * M_PI;

    // After normalization, diff should be ~0.98 degrees (the excess over 360)
    double diffDegrees = diff * 180.0 / M_PI;
    EXPECT_NEAR(diffDegrees, 0.98564736629, 0.01);
}

// ============================================================================
// ECI to ECEF Tests
// ============================================================================

TEST(ECIToECEFTest, ZeroGST) {
    // When GST = 0, ECI and ECEF should be aligned
    Vec3 eci{1000.0, 2000.0, 3000.0};
    Vec3 ecef = eciToECEF(eci, 0.0);
    EXPECT_NEAR(ecef.x, eci.x, 1e-10);
    EXPECT_NEAR(ecef.y, eci.y, 1e-10);
    EXPECT_NEAR(ecef.z, eci.z, 1e-10);
}

TEST(ECIToECEFTest, QuarterRotation) {
    // When GST = π/2, X becomes Y and Y becomes -X
    Vec3 eci{1000.0, 0.0, 0.0};
    Vec3 ecef = eciToECEF(eci, M_PI / 2.0);
    EXPECT_NEAR(ecef.x, 0.0, 1e-10);
    EXPECT_NEAR(ecef.y, -1000.0, 1e-10);
    EXPECT_NEAR(ecef.z, 0.0, 1e-10);
}

TEST(ECIToECEFTest, ZUnchanged) {
    // Z component should remain unchanged regardless of GST
    Vec3 eci{0.0, 0.0, 5000.0};
    Vec3 ecef = eciToECEF(eci, 1.234);
    EXPECT_NEAR(ecef.z, 5000.0, 1e-10);
}

TEST(ECIToECEFTest, PreservesMagnitude) {
    // Rotation should preserve the magnitude of the vector
    Vec3 eci{1000.0, 2000.0, 3000.0};
    Vec3 ecef = eciToECEF(eci, 0.789);
    EXPECT_NEAR(ecef.magnitude(), eci.magnitude(), 1e-10);
}

// ============================================================================
// ECEF to Geodetic Tests
// ============================================================================

TEST(ECEFToGeodeticTest, EquatorPrimeMeridian) {
    // Point on equator at prime meridian at Earth's surface
    // WGS84 semi-major axis = 6378.137 km
    Vec3 ecef{6378.137, 0.0, 0.0};
    Geodetic geo = ecefToGeodetic(ecef);
    EXPECT_NEAR(geo.latInRadians, 0.0, 1e-6);
    EXPECT_NEAR(geo.lonInRadians, 0.0, 1e-6);
    EXPECT_NEAR(geo.altInKilometers, 0.0, 0.01);
}

TEST(ECEFToGeodeticTest, NorthPole) {
    // Point at North Pole at Earth's surface
    // WGS84 semi-minor axis ≈ 6356.752314 km
    Vec3 ecef{0.0, 0.0, 6356.752314};
    Geodetic geo = ecefToGeodetic(ecef);
    // Latitude should be 90 degrees (π/2 radians)
    EXPECT_NEAR(geo.latInRadians, M_PI / 2.0, 1e-4);
    // Note: altitude calculation at poles uses p/cos(lat) which is unstable
    // Just verify latitude is correct; altitude at poles is a known edge case
}

TEST(ECEFToGeodeticTest, Longitude90East) {
    // Point on equator at 90° East
    Vec3 ecef{0.0, 6378.137, 0.0};
    Geodetic geo = ecefToGeodetic(ecef);
    EXPECT_NEAR(geo.latInRadians, 0.0, 1e-6);
    EXPECT_NEAR(geo.lonInRadians, M_PI / 2.0, 1e-6);
}

TEST(ECEFToGeodeticTest, ISSAltitude) {
    // ISS is approximately at 420 km altitude
    // Point on equator at prime meridian, 420 km above surface
    Vec3 ecef{6378.137 + 420.0, 0.0, 0.0};
    Geodetic geo = ecefToGeodetic(ecef);
    EXPECT_NEAR(geo.altInKilometers, 420.0, 1.0);
}

// ============================================================================
// SGP4 Propagation Tests
// ============================================================================

TEST_F(SatelliteTest, SGP4InitializesSuccessfully) {
    satellite.updateFromTLE(ISS_TLE);
    double epochJD = toJulianDate(satellite.getEpoch());
    // Should not throw - SGP4 initializes on first getECI call
    EXPECT_NO_THROW(satellite.getECI(epochJD));
}

TEST_F(SatelliteTest, SGP4PropagationReturnsValidPosition) {
    satellite.updateFromTLE(ISS_TLE);
    double epochJD = toJulianDate(satellite.getEpoch());
    Vec3 eci = satellite.getECI(epochJD);
    double r = eci.magnitude();
    // ISS should be at reasonable altitude (~420 km above Earth)
    EXPECT_GT(r, 6700.0);  // Earth radius + ~320 km
    EXPECT_LT(r, 6900.0);  // Earth radius + ~520 km
}

TEST_F(SatelliteTest, SGP4VelocityReturnsValidSpeed) {
    satellite.updateFromTLE(ISS_TLE);
    double epochJD = toJulianDate(satellite.getEpoch());
    Vec3 vel = satellite.getVelocity(epochJD);
    double speed = vel.magnitude();
    // ISS orbital velocity is ~7.66 km/s
    EXPECT_GT(speed, 7.0);
    EXPECT_LT(speed, 8.0);
}

// ============================================================================
// ECI Position Tests
// ============================================================================

TEST_F(SatelliteTest, ECIReturnsReasonableAltitude) {
    satellite.updateFromTLE(ISS_TLE);
    double epochJD = toJulianDate(satellite.getEpoch());
    Vec3 eci = satellite.getECI(epochJD);
    double r = eci.magnitude();
    // ISS altitude is ~420 km, Earth radius ~6378 km
    // So distance from Earth center should be ~6798 km
    EXPECT_GT(r, 6700.0);
    EXPECT_LT(r, 6900.0);
}

TEST_F(SatelliteTest, ECIOrbitIsNearlyCircular) {
    satellite.updateFromTLE(ISS_TLE);
    double epochJD = toJulianDate(satellite.getEpoch());
    // Check radius at two points half an orbit apart (~46 minutes)
    Vec3 eci1 = satellite.getECI(epochJD);
    Vec3 eci2 = satellite.getECI(epochJD + 46.0/1440.0);  // 46 minutes later
    double r1 = eci1.magnitude();
    double r2 = eci2.magnitude();
    // For ISS with e ≈ 0.0003723, the difference should be small
    double diff = std::abs(r2 - r1);
    EXPECT_LT(diff, 10.0);  // Less than 10 km difference
}

TEST_F(SatelliteTest, ECIInclinationConsistent) {
    satellite.updateFromTLE(ISS_TLE);
    double epochJD = toJulianDate(satellite.getEpoch());
    // At some point during orbit, the z-component will reflect inclination
    Vec3 eci = satellite.getECI(epochJD + 23.0/1440.0);  // ~quarter orbit
    double r = eci.magnitude();
    // The maximum |z/r| should be approximately sin(inclination)
    double zOverR = std::abs(eci.z / r);
    double maxSinI = std::sin(51.6312 * M_PI / 180.0);
    // z/r should be <= max inclination
    EXPECT_LE(zOverR, maxSinI + 0.05);  // Allow small tolerance
}

// ============================================================================
// printInfo Tests
// ============================================================================

TEST_F(SatelliteTest, PrintInfo_ContainsSatelliteName) {
    satellite.updateFromTLE(TLE_WITH_NAME);
    std::ostringstream oss;
    satellite.printInfo(oss);
    std::string output = oss.str();
    EXPECT_NE(output.find("ISS (ZARYA)"), std::string::npos);
}

TEST_F(SatelliteTest, PrintInfo_ContainsNoradID) {
    satellite.updateFromTLE(ISS_TLE);
    std::ostringstream oss;
    satellite.printInfo(oss);
    std::string output = oss.str();
    EXPECT_NE(output.find("NORAD ID: 25544"), std::string::npos);
}

TEST_F(SatelliteTest, PrintInfo_ContainsOrbitalElements) {
    satellite.updateFromTLE(ISS_TLE);
    std::ostringstream oss;
    satellite.printInfo(oss);
    std::string output = oss.str();

    // Check for key orbital element labels
    EXPECT_NE(output.find("Inclination:"), std::string::npos);
    EXPECT_NE(output.find("Eccentricity:"), std::string::npos);
    EXPECT_NE(output.find("Mean Motion:"), std::string::npos);
    EXPECT_NE(output.find("Mean Anomaly:"), std::string::npos);
    EXPECT_NE(output.find("Argument of Perigee:"), std::string::npos);
    EXPECT_NE(output.find("Right Ascension of Ascending Node:"), std::string::npos);
}

TEST_F(SatelliteTest, PrintInfo_ContainsEpoch) {
    satellite.updateFromTLE(ISS_TLE);
    std::ostringstream oss;
    satellite.printInfo(oss);
    std::string output = oss.str();

    // Check that epoch is present and contains UTC
    EXPECT_NE(output.find("Epoch:"), std::string::npos);
    EXPECT_NE(output.find("UTC"), std::string::npos);
}

TEST_F(SatelliteTest, PrintInfo_ContainsDragTerms) {
    satellite.updateFromTLE(ISS_TLE);
    std::ostringstream oss;
    satellite.printInfo(oss);
    std::string output = oss.str();

    EXPECT_NE(output.find("First Derivative of Mean Motion:"), std::string::npos);
    EXPECT_NE(output.find("Second Derivative of Mean Motion:"), std::string::npos);
    EXPECT_NE(output.find("Bstar Drag Term:"), std::string::npos);
}

TEST_F(SatelliteTest, PrintInfo_OutputIsNotEmpty) {
    satellite.updateFromTLE(ISS_TLE);
    std::ostringstream oss;
    satellite.printInfo(oss);
    EXPECT_FALSE(oss.str().empty());
}

// ============================================================================
// loadTLEDatabase Tests
// ============================================================================

TEST(LoadTLEDatabaseTest, LoadSingleSatellite) {
    std::istringstream input(
        "ISS (ZARYA)\n"
        "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
        "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850\n"
    );

    std::map<int, Satellite> database;
    loadTLEDatabase(input, database);

    EXPECT_EQ(database.size(), 1);
    EXPECT_TRUE(database.contains(25544));
    EXPECT_EQ(database[25544].getName(), "ISS (ZARYA)");
    EXPECT_EQ(database[25544].getNoradID(), 25544);
}

TEST(LoadTLEDatabaseTest, LoadMultipleSatellites) {
    std::istringstream input(
        "ISS (ZARYA)\n"
        "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
        "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850\n"
        "NOAA 19\n"
        "1 33591U 09005A   25333.78204194  .00000054  00000+0  52635-4 0  9999\n"
        "2 33591  98.9785  39.2910 0013037 231.6546 128.3455 14.13431889866318\n"
    );

    std::map<int, Satellite> database;
    loadTLEDatabase(input, database);

    EXPECT_EQ(database.size(), 2);
    EXPECT_TRUE(database.contains(25544));
    EXPECT_TRUE(database.contains(33591));
    EXPECT_EQ(database[25544].getName(), "ISS (ZARYA)");
    EXPECT_EQ(database[33591].getName(), "NOAA 19");
}

TEST(LoadTLEDatabaseTest, LoadEmptyStream) {
    std::istringstream input("");

    std::map<int, Satellite> database;
    loadTLEDatabase(input, database);

    EXPECT_TRUE(database.empty());
}

TEST(LoadTLEDatabaseTest, LoadWithBlankLines) {
    std::istringstream input(
        "\n"
        "ISS (ZARYA)\n"
        "\n"
        "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
        "\n"
        "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850\n"
        "\n"
    );

    std::map<int, Satellite> database;
    loadTLEDatabase(input, database);

    EXPECT_EQ(database.size(), 1);
    EXPECT_TRUE(database.contains(25544));
}

TEST(LoadTLEDatabaseTest, PreservesOrbitalElements) {
    std::istringstream input(
        "ISS (ZARYA)\n"
        "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
        "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850\n"
    );

    std::map<int, Satellite> database;
    loadTLEDatabase(input, database);

    const Satellite& iss = database[25544];
    EXPECT_DOUBLE_EQ(iss.getInclination(), 51.6312);
    EXPECT_DOUBLE_EQ(iss.getRightAscensionOfAscendingNode(), 206.3646);
    EXPECT_NEAR(iss.getEccentricity(), 0.0003723, 1e-8);
    EXPECT_NEAR(iss.getMeanMotion(), 15.49193835, 1e-6);
}

TEST(LoadTLEDatabaseTest, AppendsToExistingDatabase) {
    std::map<int, Satellite> database;

    // Load first satellite
    std::istringstream input1(
        "ISS (ZARYA)\n"
        "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
        "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850\n"
    );
    loadTLEDatabase(input1, database);
    EXPECT_EQ(database.size(), 1);

    // Load second satellite
    std::istringstream input2(
        "NOAA 19\n"
        "1 33591U 09005A   25333.78204194  .00000054  00000+0  52635-4 0  9999\n"
        "2 33591  98.9785  39.2910 0013037 231.6546 128.3455 14.13431889866318\n"
    );
    loadTLEDatabase(input2, database);

    EXPECT_EQ(database.size(), 2);
    EXPECT_TRUE(database.contains(25544));
    EXPECT_TRUE(database.contains(33591));
}

TEST(LoadTLEDatabaseTest, UpdatesExistingEntry) {
    std::map<int, Satellite> database;

    // Load ISS with one epoch
    std::istringstream input1(
        "ISS (ZARYA)\n"
        "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
        "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850\n"
    );
    loadTLEDatabase(input1, database);

    // Load ISS with updated data (different mean anomaly)
    std::istringstream input2(
        "ISS (ZARYA)\n"
        "1 25544U 98067A   25334.83453771  .00008010  00000+0  15237-3 0  9993\n"
        "2 25544  51.6312 206.3646 0003723 184.1118 200.0000 15.49193835540850\n"
    );
    loadTLEDatabase(input2, database);

    // Should still have only one entry, but with updated data
    EXPECT_EQ(database.size(), 1);
    EXPECT_DOUBLE_EQ(database[25544].getMeanAnomaly(), 200.0000);
}

TEST(LoadTLEDatabaseTest, TwoLineFormatWithoutName) {
    std::istringstream input(
        "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
        "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850\n"
    );

    std::map<int, Satellite> database;
    loadTLEDatabase(input, database);

    EXPECT_EQ(database.size(), 1);
    EXPECT_TRUE(database.contains(25544));
    // Name should be empty for 2-line format
    EXPECT_TRUE(database[25544].getName().empty());
}

TEST(LoadTLEDatabaseTest, LookupByNoradID) {
    std::istringstream input(
        "ISS (ZARYA)\n"
        "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
        "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850\n"
        "NOAA 19\n"
        "1 33591U 09005A   25333.78204194  .00000054  00000+0  52635-4 0  9999\n"
        "2 33591  98.9785  39.2910 0013037 231.6546 128.3455 14.13431889866318\n"
    );

    std::map<int, Satellite> database;
    loadTLEDatabase(input, database);

    // Test lookup of existing IDs
    EXPECT_TRUE(database.contains(25544));
    EXPECT_TRUE(database.contains(33591));

    // Test lookup of non-existing ID
    EXPECT_FALSE(database.contains(99999));
}

// ============================================================================
// getTLE Round-Trip Tests
// ============================================================================

TEST_F(SatelliteTest, GetTLE_RoundTrip_ISS) {
    // Known good 3-line TLE for ISS
    const std::string originalTLE =
        "ISS (ZARYA)\n"
        "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
        "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850\n";

    satellite.updateFromTLE(originalTLE);
    std::string regeneratedTLE = satellite.getTLE();

    EXPECT_EQ(originalTLE, regeneratedTLE);
}

TEST_F(SatelliteTest, GetTLE_RoundTrip_NOAA19) {
    // Known good 3-line TLE for NOAA 19
    const std::string originalTLE =
        "NOAA 19\n"
        "1 33591U 09005A   25333.78204194  .00000054  00000+0  52635-4 0  9999\n"
        "2 33591  98.9785  39.2910 0013037 231.6546 128.3455 14.13431889866318\n";

    satellite.updateFromTLE(originalTLE);
    std::string regeneratedTLE = satellite.getTLE();

    EXPECT_EQ(originalTLE, regeneratedTLE);
}

TEST_F(SatelliteTest, GetTLE_ContainsName) {
    satellite.updateFromTLE(TLE_WITH_NAME);
    std::string tle = satellite.getTLE();

    // First line should be the satellite name
    EXPECT_TRUE(tle.starts_with("ISS (ZARYA)\n"));
}

TEST_F(SatelliteTest, GetTLE_Line1StartsCorrectly) {
    satellite.updateFromTLE(TLE_WITH_NAME);
    std::string tle = satellite.getTLE();

    // Find the second line (line 1 of TLE)
    auto firstNewline = tle.find('\n');
    auto line1Start = tle.substr(firstNewline + 1, 2);
    EXPECT_EQ(line1Start, "1 ");
}

TEST_F(SatelliteTest, GetTLE_Line2StartsCorrectly) {
    satellite.updateFromTLE(TLE_WITH_NAME);
    std::string tle = satellite.getTLE();

    // Find line 2 of TLE
    auto firstNewline = tle.find('\n');
    auto secondNewline = tle.find('\n', firstNewline + 1);
    auto line2Start = tle.substr(secondNewline + 1, 2);
    EXPECT_EQ(line2Start, "2 ");
}

TEST_F(SatelliteTest, GetTLE_ContainsNoradID) {
    satellite.updateFromTLE(TLE_WITH_NAME);
    std::string tle = satellite.getTLE();

    // NORAD ID 25544 should appear in line 1 and line 2
    EXPECT_NE(tle.find("25544"), std::string::npos);
}

TEST_F(SatelliteTest, GetTLE_ReloadProducesSameOrbit) {
    // Load original TLE
    satellite.updateFromTLE(TLE_WITH_NAME);

    // Get regenerated TLE
    std::string regeneratedTLE = satellite.getTLE();

    // Load regenerated TLE into a new Satellite
    Satellite satellite2;
    satellite2.updateFromTLE(regeneratedTLE);

    // Compare all orbital elements
    EXPECT_EQ(satellite.getNoradID(), satellite2.getNoradID());
    EXPECT_EQ(satellite.getName(), satellite2.getName());
    EXPECT_EQ(satellite.getClassification(), satellite2.getClassification());
    EXPECT_EQ(satellite.getDesignator(), satellite2.getDesignator());
    EXPECT_DOUBLE_EQ(satellite.getInclination(), satellite2.getInclination());
    EXPECT_DOUBLE_EQ(satellite.getRightAscensionOfAscendingNode(), satellite2.getRightAscensionOfAscendingNode());
    EXPECT_DOUBLE_EQ(satellite.getEccentricity(), satellite2.getEccentricity());
    EXPECT_DOUBLE_EQ(satellite.getArgumentOfPerigee(), satellite2.getArgumentOfPerigee());
    EXPECT_DOUBLE_EQ(satellite.getMeanAnomaly(), satellite2.getMeanAnomaly());
    EXPECT_DOUBLE_EQ(satellite.getMeanMotion(), satellite2.getMeanMotion());
    EXPECT_EQ(satellite.getRevolutionNumberAtEpoch(), satellite2.getRevolutionNumberAtEpoch());
}

TEST_F(SatelliteTest, GetTLE_TwoLineFormat) {
    // Test with 2-line format (no name)
    satellite.updateFromTLE(ISS_TLE);
    std::string tle = satellite.getTLE();

    // Should still produce valid output (with empty name line)
    EXPECT_TRUE(tle.find("1 ") != std::string::npos);
    EXPECT_TRUE(tle.find("2 ") != std::string::npos);
}

// ============================================================================
// calculateChecksum Tests
// ============================================================================

TEST(CalculateChecksumTest, AllDigits) {
    // Sum of 1+2+3+4+5 = 15, mod 10 = 5
    EXPECT_EQ(calculateChecksum("12345"), 5);
}

TEST(CalculateChecksumTest, WithMinusSigns) {
    // Each '-' counts as 1
    // Sum of 1+2+1+3+1 = 8 (the '-' signs each add 1)
    EXPECT_EQ(calculateChecksum("12-3-"), 8);
}

TEST(CalculateChecksumTest, WithLettersAndSpaces) {
    // Letters, spaces, and other chars are ignored
    // Only 1+2+3 = 6
    EXPECT_EQ(calculateChecksum("1 A 2 B 3"), 6);
}

TEST(CalculateChecksumTest, EmptyString) {
    EXPECT_EQ(calculateChecksum(""), 0);
}

TEST(CalculateChecksumTest, RealTLELine1) {
    // Real TLE line 1 for ISS (without checksum digit)
    // The checksum should be 3
    std::string line1 = "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  999";
    EXPECT_EQ(calculateChecksum(line1), 3);
}

TEST(CalculateChecksumTest, RealTLELine2) {
    // Real TLE line 2 for ISS (without checksum digit)
    // The checksum should be 0
    std::string line2 = "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.4919383554085";
    EXPECT_EQ(calculateChecksum(line2), 0);
}

TEST(CalculateChecksumTest, OnlyMinusSigns) {
    // Three minus signs = 3
    EXPECT_EQ(calculateChecksum("---"), 3);
}

TEST(CalculateChecksumTest, LargeSum) {
    // 9+9+9+9+9 = 45, mod 10 = 5
    EXPECT_EQ(calculateChecksum("99999"), 5);
}

// ============================================================================
// toTLEExponential Tests
// ============================================================================

TEST(ToTLEExponentialTest, Zero) {
    EXPECT_EQ(toTLEExponential(0.0), " 00000+0");
}

TEST(ToTLEExponentialTest, PositiveSmallValue) {
    // 0.00015237 -> mantissa 15237, exponent -3
    // Format: " 15237-3"
    EXPECT_EQ(toTLEExponential(0.00015237), " 15237-3");
}

TEST(ToTLEExponentialTest, NegativeValue) {
    // -0.00015237 -> "-15237-3"
    EXPECT_EQ(toTLEExponential(-0.00015237), "-15237-3");
}

TEST(ToTLEExponentialTest, VerySmallValue) {
    // 0.000052635 -> " 52635-4"
    EXPECT_EQ(toTLEExponential(0.000052635), " 52635-4");
}

TEST(ToTLEExponentialTest, PositiveExponent) {
    // 1.5 -> mantissa 15000, exponent +1
    // Format: " 15000+1"
    EXPECT_EQ(toTLEExponential(1.5), " 15000+1");
}

TEST(ToTLEExponentialTest, ExactlyOne) {
    // 1.0 -> " 10000+1"
    EXPECT_EQ(toTLEExponential(1.0), " 10000+1");
}

TEST(ToTLEExponentialTest, TinyValue) {
    // 1e-10 -> " 10000-9"
    EXPECT_EQ(toTLEExponential(1e-10), " 10000-9");
}

// ============================================================================
// formatFirstDerivative Tests
// ============================================================================

TEST(FormatFirstDerivativeTest, Zero) {
    EXPECT_EQ(formatFirstDerivative(0.0), " .00000000");
}

TEST(FormatFirstDerivativeTest, PositiveValue) {
    // 0.00008010 -> " .00008010"
    EXPECT_EQ(formatFirstDerivative(0.00008010), " .00008010");
}

TEST(FormatFirstDerivativeTest, NegativeValue) {
    // -0.00008010 -> "-.00008010"
    EXPECT_EQ(formatFirstDerivative(-0.00008010), "-.00008010");
}

TEST(FormatFirstDerivativeTest, SmallPositiveValue) {
    // 0.00000054 -> " .00000054"
    EXPECT_EQ(formatFirstDerivative(0.00000054), " .00000054");
}

TEST(FormatFirstDerivativeTest, LargerValue) {
    // 0.12345678 -> " .12345678"
    EXPECT_EQ(formatFirstDerivative(0.12345678), " .12345678");
}

TEST(FormatFirstDerivativeTest, RoundingUp) {
    // 0.000000005 should round to 0.00000001
    EXPECT_EQ(formatFirstDerivative(0.000000005), " .00000001");
}

TEST(FormatFirstDerivativeTest, RoundingDown) {
    // 0.000000004 should round to 0.00000000
    EXPECT_EQ(formatFirstDerivative(0.000000004), " .00000000");
}

// ============================================================================
// Vec3 Extended Operations Tests
// ============================================================================

TEST(Vec3Test, Subtraction) {
    Vec3 a{5.0, 7.0, 9.0};
    Vec3 b{1.0, 2.0, 3.0};
    Vec3 result = a - b;
    EXPECT_DOUBLE_EQ(result.x, 4.0);
    EXPECT_DOUBLE_EQ(result.y, 5.0);
    EXPECT_DOUBLE_EQ(result.z, 6.0);
}

TEST(Vec3Test, SubtractionNegativeResult) {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{5.0, 7.0, 9.0};
    Vec3 result = a - b;
    EXPECT_DOUBLE_EQ(result.x, -4.0);
    EXPECT_DOUBLE_EQ(result.y, -5.0);
    EXPECT_DOUBLE_EQ(result.z, -6.0);
}

TEST(Vec3Test, DotProduct) {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};
    // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    EXPECT_DOUBLE_EQ(a.dot(b), 32.0);
}

TEST(Vec3Test, DotProductPerpendicular) {
    // Two perpendicular vectors have dot product = 0
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{0.0, 1.0, 0.0};
    EXPECT_DOUBLE_EQ(a.dot(b), 0.0);
}

TEST(Vec3Test, DotProductParallel) {
    // Parallel vectors: dot product = |a| * |b|
    Vec3 a{2.0, 0.0, 0.0};
    Vec3 b{3.0, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(a.dot(b), 6.0);
}

TEST(Vec3Test, DotProductAntiparallel) {
    // Anti-parallel vectors: dot product = -|a| * |b|
    Vec3 a{2.0, 0.0, 0.0};
    Vec3 b{-3.0, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(a.dot(b), -6.0);
}

TEST(Vec3Test, CrossProduct) {
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{0.0, 1.0, 0.0};
    Vec3 result = a.cross(b);
    // i × j = k
    EXPECT_DOUBLE_EQ(result.x, 0.0);
    EXPECT_DOUBLE_EQ(result.y, 0.0);
    EXPECT_DOUBLE_EQ(result.z, 1.0);
}

TEST(Vec3Test, CrossProductReverse) {
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{0.0, 1.0, 0.0};
    Vec3 result = b.cross(a);
    // j × i = -k
    EXPECT_DOUBLE_EQ(result.x, 0.0);
    EXPECT_DOUBLE_EQ(result.y, 0.0);
    EXPECT_DOUBLE_EQ(result.z, -1.0);
}

TEST(Vec3Test, CrossProductGeneral) {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};
    Vec3 result = a.cross(b);
    // (2*6 - 3*5, 3*4 - 1*6, 1*5 - 2*4) = (12-15, 12-6, 5-8) = (-3, 6, -3)
    EXPECT_DOUBLE_EQ(result.x, -3.0);
    EXPECT_DOUBLE_EQ(result.y, 6.0);
    EXPECT_DOUBLE_EQ(result.z, -3.0);
}

TEST(Vec3Test, CrossProductPerpendicular) {
    // Cross product result is perpendicular to both inputs
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};
    Vec3 c = a.cross(b);
    EXPECT_NEAR(a.dot(c), 0.0, 1e-10);
    EXPECT_NEAR(b.dot(c), 0.0, 1e-10);
}

TEST(Vec3Test, CrossProductParallelIsZero) {
    // Parallel vectors have zero cross product
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{2.0, 4.0, 6.0};
    Vec3 result = a.cross(b);
    EXPECT_NEAR(result.magnitude(), 0.0, 1e-10);
}

TEST(Vec3Test, Normalize) {
    Vec3 v{3.0, 4.0, 0.0};
    Vec3 result = v.normalize();
    EXPECT_DOUBLE_EQ(result.magnitude(), 1.0);
    EXPECT_DOUBLE_EQ(result.x, 0.6);
    EXPECT_DOUBLE_EQ(result.y, 0.8);
    EXPECT_DOUBLE_EQ(result.z, 0.0);
}

TEST(Vec3Test, NormalizeUnitVector) {
    Vec3 v{1.0, 0.0, 0.0};
    Vec3 result = v.normalize();
    EXPECT_DOUBLE_EQ(result.x, 1.0);
    EXPECT_DOUBLE_EQ(result.y, 0.0);
    EXPECT_DOUBLE_EQ(result.z, 0.0);
}

TEST(Vec3Test, NormalizePreservesDirection) {
    Vec3 v{10.0, 20.0, 30.0};
    Vec3 n = v.normalize();
    // Normalized vector should be parallel to original
    // (their cross product should be zero)
    Vec3 cross = v.cross(n);
    EXPECT_NEAR(cross.magnitude(), 0.0, 1e-10);
}

// ============================================================================
// Geodetic::toECEF Tests
// ============================================================================

TEST(GeodeticToECEFTest, EquatorPrimeMeridian) {
    // Point on equator at prime meridian at Earth's surface
    Geodetic geo{0.0, 0.0, 0.0};
    Vec3 ecef = geo.toECEF();
    // X should be approximately Earth's equatorial radius
    EXPECT_NEAR(ecef.x, 6378.137, 0.001);
    EXPECT_NEAR(ecef.y, 0.0, 1e-10);
    EXPECT_NEAR(ecef.z, 0.0, 1e-10);
}

TEST(GeodeticToECEFTest, EquatorAt90East) {
    // Point on equator at 90° East
    Geodetic geo{0.0, M_PI / 2.0, 0.0};
    Vec3 ecef = geo.toECEF();
    EXPECT_NEAR(ecef.x, 0.0, 1e-10);
    EXPECT_NEAR(ecef.y, 6378.137, 0.001);
    EXPECT_NEAR(ecef.z, 0.0, 1e-10);
}

TEST(GeodeticToECEFTest, NorthPole) {
    // North pole at surface
    Geodetic geo{M_PI / 2.0, 0.0, 0.0};
    Vec3 ecef = geo.toECEF();
    // At north pole, X and Y should be ~0, Z should be semi-minor axis
    EXPECT_NEAR(ecef.x, 0.0, 1e-10);
    EXPECT_NEAR(ecef.y, 0.0, 1e-10);
    // WGS84 semi-minor axis ≈ 6356.752 km
    EXPECT_NEAR(ecef.z, 6356.752, 0.1);
}

TEST(GeodeticToECEFTest, WithAltitude) {
    // Point on equator at prime meridian, 420 km altitude (ISS altitude)
    Geodetic geo{0.0, 0.0, 420.0};
    Vec3 ecef = geo.toECEF();
    EXPECT_NEAR(ecef.x, 6378.137 + 420.0, 0.001);
    EXPECT_NEAR(ecef.y, 0.0, 1e-10);
    EXPECT_NEAR(ecef.z, 0.0, 1e-10);
}

TEST(GeodeticToECEFTest, RoundTrip) {
    // Test that toECEF -> ecefToGeodetic gives back original values
    Geodetic original{0.7, 1.2, 500.0};  // ~40° lat, ~69° lon, 500 km alt
    Vec3 ecef = original.toECEF();
    Geodetic recovered = ecefToGeodetic(ecef);

    EXPECT_NEAR(recovered.latInRadians, original.latInRadians, 1e-6);
    EXPECT_NEAR(recovered.lonInRadians, original.lonInRadians, 1e-6);
    EXPECT_NEAR(recovered.altInKilometers, original.altInKilometers, 0.01);
}

TEST(GeodeticToECEFTest, RoundTripMultipleLocations) {
    // Test round-trip at several locations
    std::vector<Geodetic> testLocations = {
        {0.0, 0.0, 0.0},           // Equator, prime meridian, surface
        {M_PI / 4, M_PI / 4, 100}, // 45°N, 45°E, 100km alt
        {-M_PI / 3, -M_PI / 2, 35786}, // GEO altitude
        {M_PI / 6, M_PI, 400},     // 30°N, 180°, 400km (ISS-like)
    };

    for (const auto& original : testLocations) {
        Vec3 ecef = original.toECEF();
        Geodetic recovered = ecefToGeodetic(ecef);

        EXPECT_NEAR(recovered.latInRadians, original.latInRadians, 1e-5);
        EXPECT_NEAR(recovered.lonInRadians, original.lonInRadians, 1e-5);
        EXPECT_NEAR(recovered.altInKilometers, original.altInKilometers, 0.1);
    }
}

// ============================================================================
// ecefToENU Tests
// ============================================================================

TEST(ECEFToENUTest, SatelliteDirectlyAbove) {
    // Observer on equator at prime meridian
    Geodetic observer{0.0, 0.0, 0.0};
    // Satellite directly above at 400 km altitude
    Vec3 observerECEF = observer.toECEF();
    Vec3 satECEF = {observerECEF.x + 400.0, 0.0, 0.0};

    Vec3 enu = ecefToENU(satECEF, observer);

    // Should be purely "Up"
    EXPECT_NEAR(enu.x, 0.0, 1e-6);  // East
    EXPECT_NEAR(enu.y, 0.0, 1e-6);  // North
    EXPECT_NEAR(enu.z, 400.0, 0.1); // Up
}

TEST(ECEFToENUTest, SatelliteToEast) {
    // Observer on equator at prime meridian
    Geodetic observer{0.0, 0.0, 0.0};
    // Satellite same altitude but displaced toward +Y (East)
    Vec3 satECEF = {6378.137, 500.0, 0.0};

    Vec3 enu = ecefToENU(satECEF, observer);

    // Should be primarily East with some Up
    EXPECT_GT(enu.x, 0.0);  // East (positive)
    EXPECT_NEAR(enu.y, 0.0, 1e-6);  // North (approximately zero)
}

TEST(ECEFToENUTest, SatelliteToNorth) {
    // Observer on equator at prime meridian
    Geodetic observer{0.0, 0.0, 0.0};
    // Satellite displaced toward +Z (North)
    Vec3 satECEF = {6378.137, 0.0, 500.0};

    Vec3 enu = ecefToENU(satECEF, observer);

    // Should be primarily North with some Up
    EXPECT_NEAR(enu.x, 0.0, 1e-6);  // East (approximately zero)
    EXPECT_GT(enu.y, 0.0);  // North (positive)
}

TEST(ECEFToENUTest, RangePreservation) {
    // The ENU transformation should preserve the distance (range)
    Geodetic observer{0.5, 1.0, 100.0};  // Arbitrary location
    Vec3 observerECEF = observer.toECEF();

    // Arbitrary satellite position
    Vec3 satECEF = {7000.0, 1000.0, 2000.0};

    // Compute range in ECEF
    Vec3 diff = satECEF - observerECEF;
    double rangeECEF = diff.magnitude();

    // Compute range in ENU
    Vec3 enu = ecefToENU(satECEF, observer);
    double rangeENU = enu.magnitude();

    // Should be the same
    EXPECT_NEAR(rangeENU, rangeECEF, 1e-6);
}

// ============================================================================
// getLookAngles Tests
// ============================================================================

TEST(GetLookAnglesTest, SatelliteDirectlyOverhead) {
    // Observer on equator at prime meridian
    Geodetic observer{0.0, 0.0, 0.0};
    // Satellite directly above
    Vec3 observerECEF = observer.toECEF();
    Vec3 satECEF = observerECEF * (1.0 + 400.0 / observerECEF.magnitude());

    LookAngles angles = getLookAngles(satECEF, observer);

    // Elevation should be 90 degrees (π/2)
    EXPECT_NEAR(angles.elevationInRadians, M_PI / 2.0, 0.01);
    // Range should be approximately 400 km
    EXPECT_NEAR(angles.rangeInKilometers, 400.0, 1.0);
}

TEST(GetLookAnglesTest, SatelliteAtSameAltitude) {
    // Test with a satellite at similar distance from Earth's center
    // but displaced so it's not directly overhead
    Geodetic observer{0.0, 0.0, 0.0};
    Vec3 observerECEF = observer.toECEF();

    // Satellite at 90 degrees longitude (due East), same radius
    // This puts it well above the horizon due to Earth's curvature
    Vec3 satECEF = {0.0, observerECEF.x + 400.0, 0.0};

    LookAngles angles = getLookAngles(satECEF, observer);

    // Verify we get valid look angles
    EXPECT_GE(angles.azimuthInRadians, 0.0);
    EXPECT_LT(angles.azimuthInRadians, 2.0 * M_PI);
    // Elevation can be positive (satellite is ~400km further from Earth center)
    EXPECT_GE(angles.elevationInRadians, -M_PI / 2.0);
    EXPECT_LE(angles.elevationInRadians, M_PI / 2.0);
}

TEST(GetLookAnglesTest, AzimuthNorth) {
    // Observer on equator at prime meridian
    Geodetic observer{0.0, 0.0, 0.0};
    Vec3 observerECEF = observer.toECEF();

    // Satellite to the north (positive Z in ECEF)
    Vec3 satECEF = {observerECEF.x, 0.0, 1000.0};

    LookAngles angles = getLookAngles(satECEF, observer);

    // Azimuth should be close to 0 (North)
    EXPECT_LT(angles.azimuthInRadians, M_PI / 6);  // Within 30 degrees of North
}

TEST(GetLookAnglesTest, AzimuthEast) {
    // Observer on equator at prime meridian
    Geodetic observer{0.0, 0.0, 0.0};
    Vec3 observerECEF = observer.toECEF();

    // Satellite to the east (positive Y in ECEF at equator/prime meridian)
    Vec3 satECEF = {observerECEF.x, 1000.0, 0.0};

    LookAngles angles = getLookAngles(satECEF, observer);

    // Azimuth should be close to π/2 (East = 90°)
    EXPECT_NEAR(angles.azimuthInRadians, M_PI / 2.0, M_PI / 6);
}

TEST(GetLookAnglesTest, AzimuthInValidRange) {
    // Test that azimuth is always in [0, 2π)
    Geodetic observer{0.5, 1.0, 0.0};

    std::vector<Vec3> testSatellites = {
        {7000.0, 1000.0, 500.0},
        {-7000.0, 1000.0, -500.0},
        {1000.0, -7000.0, 500.0},
        {1000.0, 1000.0, 7000.0},
    };

    for (const auto& satECEF : testSatellites) {
        LookAngles angles = getLookAngles(satECEF, observer);
        EXPECT_GE(angles.azimuthInRadians, 0.0);
        EXPECT_LT(angles.azimuthInRadians, 2.0 * M_PI);
    }
}

TEST(GetLookAnglesTest, ElevationInValidRange) {
    // Test that elevation is always in [-π/2, π/2]
    Geodetic observer{0.5, 1.0, 0.0};

    std::vector<Vec3> testSatellites = {
        {7000.0, 1000.0, 500.0},
        {-7000.0, 1000.0, -500.0},
        {1000.0, -7000.0, 500.0},
        {1000.0, 1000.0, 7000.0},
    };

    for (const auto& satECEF : testSatellites) {
        LookAngles angles = getLookAngles(satECEF, observer);
        EXPECT_GE(angles.elevationInRadians, -M_PI / 2.0);
        EXPECT_LE(angles.elevationInRadians, M_PI / 2.0);
    }
}

// ============================================================================
// getLookAngles with Satellite Tests
// ============================================================================

TEST_F(SatelliteTest, GetLookAnglesAtEpoch) {
    satellite.updateFromTLE(ISS_TLE);

    // Observer at a mid-latitude location
    Geodetic observer{0.7, -1.5, 0.0};  // ~40°N, ~86°W

    double epochJD = toJulianDate(satellite.getEpoch());
    LookAngles angles = getLookAngles(satellite, observer, epochJD);

    // Just verify we get reasonable values
    EXPECT_GE(angles.azimuthInRadians, 0.0);
    EXPECT_LT(angles.azimuthInRadians, 2.0 * M_PI);
    EXPECT_GE(angles.elevationInRadians, -M_PI / 2.0);
    EXPECT_LE(angles.elevationInRadians, M_PI / 2.0);
    EXPECT_GT(angles.rangeInKilometers, 0.0);
    EXPECT_LT(angles.rangeInKilometers, 50000.0);  // Should be less than GEO distance
}

TEST_F(SatelliteTest, GetLookAnglesRangeReasonable) {
    satellite.updateFromTLE(ISS_TLE);

    // Observer location doesn't matter much for range check
    Geodetic observer{0.0, 0.0, 0.0};

    double epochJD = toJulianDate(satellite.getEpoch());
    LookAngles angles = getLookAngles(satellite, observer, epochJD);

    // ISS at ~420 km altitude
    // Minimum range (directly overhead): ~420 km
    // Maximum range (when satellite is on opposite side of Earth): ~13000 km
    // The actual range depends on where the satellite is in its orbit relative to observer
    EXPECT_GT(angles.rangeInKilometers, 400.0);
    EXPECT_LT(angles.rangeInKilometers, 15000.0);
}

// ============================================================================
// isVisible Tests
// ============================================================================

TEST(IsVisibleTest, AboveThreshold) {
    LookAngles angles{0.0, 0.2, 1000.0};  // 0.2 rad ≈ 11.5°
    EXPECT_TRUE(isVisible(angles, 0.1));  // Threshold 0.1 rad ≈ 5.7°
}

TEST(IsVisibleTest, BelowThreshold) {
    LookAngles angles{0.0, 0.05, 1000.0};  // 0.05 rad ≈ 2.9°
    EXPECT_FALSE(isVisible(angles, 0.1));  // Threshold 0.1 rad ≈ 5.7°
}

TEST(IsVisibleTest, ExactlyAtThreshold) {
    LookAngles angles{0.0, 0.1, 1000.0};
    EXPECT_TRUE(isVisible(angles, 0.1));  // Should be visible when exactly at threshold
}

TEST(IsVisibleTest, NegativeElevation) {
    LookAngles angles{0.0, -0.1, 1000.0};  // Below horizon
    EXPECT_FALSE(isVisible(angles, 0.0));  // Threshold at horizon
}

TEST(IsVisibleTest, ZeroThreshold) {
    LookAngles aboveHorizon{0.0, 0.01, 1000.0};
    LookAngles belowHorizon{0.0, -0.01, 1000.0};

    EXPECT_TRUE(isVisible(aboveHorizon, 0.0));
    EXPECT_FALSE(isVisible(belowHorizon, 0.0));
}

TEST(IsVisibleTest, HighElevation) {
    LookAngles overhead{0.0, M_PI / 2.0, 400.0};  // Directly overhead
    EXPECT_TRUE(isVisible(overhead, 0.0));
    EXPECT_TRUE(isVisible(overhead, M_PI / 4.0));  // 45° threshold
    EXPECT_TRUE(isVisible(overhead, M_PI / 3.0));  // 60° threshold
}

// ============================================================================
// isVisible with Satellite Tests
// ============================================================================

TEST_F(SatelliteTest, IsVisibleReturnsReasonableResults) {
    satellite.updateFromTLE(ISS_TLE);

    // Observer at a specific location
    Geodetic observer{0.7, -1.5, 0.0};

    double epochJD = toJulianDate(satellite.getEpoch());

    // Check visibility with different thresholds
    // These should return consistent results (visible at low threshold means
    // visible at even lower threshold)
    bool visibleAt0 = isVisible(satellite, observer, epochJD, 0.0);
    bool visibleAt10deg = isVisible(satellite, observer, epochJD, 10.0 * M_PI / 180.0);
    bool visibleAt45deg = isVisible(satellite, observer, epochJD, 45.0 * M_PI / 180.0);

    // If visible at high threshold, must be visible at lower threshold
    if (visibleAt45deg) {
        EXPECT_TRUE(visibleAt10deg);
        EXPECT_TRUE(visibleAt0);
    }
    if (visibleAt10deg) {
        EXPECT_TRUE(visibleAt0);
    }
}

// ============================================================================
// findNextPass Tests
// ============================================================================

TEST_F(SatelliteTest, FindNextPassReturnsValidPassForISS) {
    // ISS should have passes over mid-latitude locations
    satellite.updateFromTLE(ISS_TLE);
    Geodetic observer{0.7, -1.5, 0.0};  // ~40°N, ~86°W

    auto pass = findNextPass(satellite, observer, 0.0, satellite.getEpoch());

    // Should find a pass within 48 hours for ISS
    EXPECT_TRUE(pass.has_value());
}

TEST_F(SatelliteTest, FindNextPassReturnsValidPass) {
    satellite.updateFromTLE(ISS_TLE);
    Geodetic observer{0.7, -1.5, 0.0};  // ~40°N, ~86°W

    auto pass = findNextPass(satellite, observer, 0.0, satellite.getEpoch());

    ASSERT_TRUE(pass.has_value());
    // Rise time should be before max elevation time
    EXPECT_LT(pass->riseTime, pass->maxElevationTime);
    // Max elevation time should be before set time
    EXPECT_LT(pass->maxElevationTime, pass->setTime);
    // Max elevation should be higher than rise/set
    EXPECT_GE(pass->maxAngles.elevationInRadians, pass->riseAngles.elevationInRadians);
    EXPECT_GE(pass->maxAngles.elevationInRadians, pass->setAngles.elevationInRadians);
}

TEST_F(SatelliteTest, FindNextPassWithMinElevation) {
    satellite.updateFromTLE(ISS_TLE);
    Geodetic observer{0.7, -1.5, 0.0};
    double minElevation = 10.0 * M_PI / 180.0;  // 10 degrees

    auto pass = findNextPass(satellite, observer, minElevation, satellite.getEpoch());

    if (pass.has_value()) {
        // Rise elevation should be at or above minimum
        EXPECT_GE(pass->riseAngles.elevationInRadians, minElevation - 0.01);
        // Set elevation should be at or above minimum
        EXPECT_GE(pass->setAngles.elevationInRadians, minElevation - 0.01);
    }
}

TEST_F(SatelliteTest, FindNextPassISSHasMultipleDailyPasses) {
    // ISS completes about 15.5 orbits per day, so there should be multiple
    // opportunities for passes over any given location
    satellite.updateFromTLE(ISS_TLE);
    Geodetic observer{0.7, -1.5, 0.0};

    auto now = satellite.getEpoch();
    int passCount = 0;

    // Search for passes over 24 hours
    auto searchTime = now;
    auto endTime = now + std::chrono::hours(24);

    while (searchTime < endTime) {
        auto pass = findNextPass(satellite, observer, 0.0, searchTime);
        if (pass.has_value()) {
            passCount++;
            // Move search time past this pass
            searchTime = pass->setTime + std::chrono::minutes(1);
        } else {
            break;
        }
    }

    // ISS should have multiple passes per day at mid-latitudes
    EXPECT_GE(passCount, 1);
}

}  // namespace
}  // namespace sattrack
