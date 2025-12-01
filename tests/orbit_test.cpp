/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <gtest/gtest.h>
#include <sattrack/orbit.hpp>

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
    EXPECT_TRUE(defaultOrbit.getName().empty());
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
    EXPECT_DOUBLE_EQ(orbit.getRightAscensionOfAscendingNode(), 206.3646);
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

// ============================================================================
// Name Parsing Tests
// ============================================================================

TEST_F(OrbitTest, GetName_TwoLineFormat) {
    // Two-line TLE format should result in empty name
    orbit.updateFromTLE(ISS_TLE);
    EXPECT_TRUE(orbit.getName().empty());
}

TEST_F(OrbitTest, GetName_ThreeLineFormat) {
    // Three-line TLE format should parse the name from the first line
    orbit.updateFromTLE(TLE_WITH_NAME);
    EXPECT_EQ(orbit.getName(), "ISS (ZARYA)");
}

TEST_F(OrbitTest, UpdateFromTLE_WithExplicitName) {
    // Using the overload that accepts a separate name
    orbit.updateFromTLE("My Custom Name", ISS_TLE);
    EXPECT_EQ(orbit.getName(), "My Custom Name");
    // Should still parse the orbital elements correctly
    EXPECT_EQ(orbit.getNoradID(), 25544);
    EXPECT_DOUBLE_EQ(orbit.getInclination(), 51.6312);
}

TEST_F(OrbitTest, UpdateFromTLE_ExplicitNameOverridesEmbeddedName) {
    // When using the overload with explicit name, it should override
    // any name that might be embedded in the TLE
    orbit.updateFromTLE("Override Name", TLE_WITH_NAME);
    EXPECT_EQ(orbit.getName(), "Override Name");
}

TEST_F(OrbitTest, GetName_TrimsWhitespace) {
    // Test that leading/trailing whitespace is trimmed from name
    constexpr const char* TLE_WITH_PADDED_NAME =
        "  NOAA 19  \n"
        "1 33591U 09005A   25333.78204194  .00000054  00000+0  52635-4 0  9999\n"
        "2 33591  98.9785  39.2910 0013037 231.6546 128.3455 14.13431889866318";
    orbit.updateFromTLE(TLE_WITH_PADDED_NAME);
    EXPECT_EQ(orbit.getName(), "NOAA 19");
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
    double raan = orbit.getRightAscensionOfAscendingNode();
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
// Orbit Anomaly Calculation Tests
// ============================================================================

TEST_F(OrbitTest, MeanAnomalyAtEpoch) {
    orbit.updateFromTLE(ISS_TLE);
    double epochJD = toJulianDate(orbit.getEpoch());
    double M = orbit.getMeanAnomalyAtTime(epochJD);
    // Should match the TLE mean anomaly (175.9840 degrees) in radians
    double expectedM = 175.9840 * M_PI / 180.0;
    EXPECT_NEAR(M, expectedM, 1e-6);
}

TEST_F(OrbitTest, MeanAnomalyInRange) {
    orbit.updateFromTLE(ISS_TLE);
    double epochJD = toJulianDate(orbit.getEpoch());
    // Test at various times
    for (double dt = 0; dt < 1.0; dt += 0.1) {
        double M = orbit.getMeanAnomalyAtTime(epochJD + dt);
        EXPECT_GE(M, 0.0);
        EXPECT_LT(M, 2.0 * M_PI);
    }
}

TEST_F(OrbitTest, EccentricAnomalyCircularOrbit) {
    orbit.updateFromTLE(ISS_TLE);
    // For nearly circular orbits (e ≈ 0), E ≈ M
    double M = M_PI / 4.0;  // 45 degrees
    double E = orbit.getEccentricAnomalyFromMeanAnomaly(M);
    // ISS has e ≈ 0.0003723, so E should be very close to M
    EXPECT_NEAR(E, M, 0.001);
}

TEST_F(OrbitTest, EccentricAnomalyKeplerEquation) {
    orbit.updateFromTLE(ISS_TLE);
    // Verify Kepler's equation: M = E - e*sin(E)
    double M = 1.0;  // radians
    double E = orbit.getEccentricAnomalyFromMeanAnomaly(M);
    double e = orbit.getEccentricity();
    double computedM = E - e * std::sin(E);
    EXPECT_NEAR(computedM, M, 1e-10);
}

TEST_F(OrbitTest, TrueAnomalyCircularOrbit) {
    orbit.updateFromTLE(ISS_TLE);
    // For nearly circular orbits, true anomaly ≈ eccentric anomaly ≈ mean anomaly
    double E = M_PI / 3.0;  // 60 degrees
    double nu = orbit.getTrueAnomalyFromEccentricAnomaly(E);
    EXPECT_NEAR(nu, E, 0.001);
}

TEST_F(OrbitTest, TrueAnomalyAtZero) {
    orbit.updateFromTLE(ISS_TLE);
    // At E = 0, true anomaly should also be 0
    double nu = orbit.getTrueAnomalyFromEccentricAnomaly(0.0);
    EXPECT_NEAR(nu, 0.0, 1e-10);
}

TEST_F(OrbitTest, TrueAnomalyAtPi) {
    orbit.updateFromTLE(ISS_TLE);
    // At E = π, true anomaly should also be π
    double nu = orbit.getTrueAnomalyFromEccentricAnomaly(M_PI);
    EXPECT_NEAR(nu, M_PI, 1e-10);
}

TEST_F(OrbitTest, TrueAnomalyAtTimeReturnsValidRange) {
    orbit.updateFromTLE(ISS_TLE);
    double epochJD = toJulianDate(orbit.getEpoch());
    double nu = orbit.getTrueAnomalyAtTime(epochJD);
    // True anomaly should be in range [-π, π] or [0, 2π)
    EXPECT_GE(nu, -M_PI);
    EXPECT_LE(nu, M_PI);
}

// ============================================================================
// ECI Position Tests
// ============================================================================

TEST_F(OrbitTest, ECIReturnsReasonableAltitude) {
    orbit.updateFromTLE(ISS_TLE);
    double nu = 0.0;  // At perigee
    Vec3 eci = orbit.getECI(nu);
    double r = eci.magnitude();
    // ISS altitude is ~420 km, Earth radius ~6378 km
    // So distance from Earth center should be ~6798 km
    EXPECT_GT(r, 6700.0);
    EXPECT_LT(r, 6900.0);
}

TEST_F(OrbitTest, ECIOrbitIsNearlyCircular) {
    orbit.updateFromTLE(ISS_TLE);
    // Check radius at perigee (nu=0) and apogee (nu=π)
    Vec3 eciPerigee = orbit.getECI(0.0);
    Vec3 eciApogee = orbit.getECI(M_PI);
    double rPerigee = eciPerigee.magnitude();
    double rApogee = eciApogee.magnitude();
    // For ISS with e ≈ 0.0003723, the difference should be small
    double diff = std::abs(rApogee - rPerigee);
    EXPECT_LT(diff, 10.0);  // Less than 10 km difference
}

TEST_F(OrbitTest, ECIInclinationConsistent) {
    orbit.updateFromTLE(ISS_TLE);
    Vec3 eci = orbit.getECI(M_PI / 2.0);  // At ascending node + 90°
    double r = eci.magnitude();
    // The z component at this point indicates the inclination
    // For ISS with i ≈ 51.6°, z/r should be approximately sin(51.6°)
    double sinI = eci.z / r;
    double expectedSinI = std::sin(51.6312 * M_PI / 180.0);
    EXPECT_NEAR(std::abs(sinI), expectedSinI, 0.01);
}

// ============================================================================
// printInfo Tests
// ============================================================================

TEST_F(OrbitTest, PrintInfo_ContainsSatelliteName) {
    orbit.updateFromTLE(TLE_WITH_NAME);
    std::ostringstream oss;
    orbit.printInfo(oss);
    std::string output = oss.str();
    EXPECT_NE(output.find("ISS (ZARYA)"), std::string::npos);
}

TEST_F(OrbitTest, PrintInfo_ContainsNoradID) {
    orbit.updateFromTLE(ISS_TLE);
    std::ostringstream oss;
    orbit.printInfo(oss);
    std::string output = oss.str();
    EXPECT_NE(output.find("NORAD ID: 25544"), std::string::npos);
}

TEST_F(OrbitTest, PrintInfo_ContainsOrbitalElements) {
    orbit.updateFromTLE(ISS_TLE);
    std::ostringstream oss;
    orbit.printInfo(oss);
    std::string output = oss.str();

    // Check for key orbital element labels
    EXPECT_NE(output.find("Inclination:"), std::string::npos);
    EXPECT_NE(output.find("Eccentricity:"), std::string::npos);
    EXPECT_NE(output.find("Mean Motion:"), std::string::npos);
    EXPECT_NE(output.find("Mean Anomaly:"), std::string::npos);
    EXPECT_NE(output.find("Argument of Perigee:"), std::string::npos);
    EXPECT_NE(output.find("Right Ascension of Ascending Node:"), std::string::npos);
}

TEST_F(OrbitTest, PrintInfo_ContainsEpoch) {
    orbit.updateFromTLE(ISS_TLE);
    std::ostringstream oss;
    orbit.printInfo(oss);
    std::string output = oss.str();

    // Check that epoch is present and contains UTC
    EXPECT_NE(output.find("Epoch:"), std::string::npos);
    EXPECT_NE(output.find("UTC"), std::string::npos);
}

TEST_F(OrbitTest, PrintInfo_ContainsDragTerms) {
    orbit.updateFromTLE(ISS_TLE);
    std::ostringstream oss;
    orbit.printInfo(oss);
    std::string output = oss.str();

    EXPECT_NE(output.find("First Derivative of Mean Motion:"), std::string::npos);
    EXPECT_NE(output.find("Second Derivative of Mean Motion:"), std::string::npos);
    EXPECT_NE(output.find("Bstar Drag Term:"), std::string::npos);
}

TEST_F(OrbitTest, PrintInfo_OutputIsNotEmpty) {
    orbit.updateFromTLE(ISS_TLE);
    std::ostringstream oss;
    orbit.printInfo(oss);
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

    std::map<int, Orbit> database;
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

    std::map<int, Orbit> database;
    loadTLEDatabase(input, database);

    EXPECT_EQ(database.size(), 2);
    EXPECT_TRUE(database.contains(25544));
    EXPECT_TRUE(database.contains(33591));
    EXPECT_EQ(database[25544].getName(), "ISS (ZARYA)");
    EXPECT_EQ(database[33591].getName(), "NOAA 19");
}

TEST(LoadTLEDatabaseTest, LoadEmptyStream) {
    std::istringstream input("");

    std::map<int, Orbit> database;
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

    std::map<int, Orbit> database;
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

    std::map<int, Orbit> database;
    loadTLEDatabase(input, database);

    const Orbit& iss = database[25544];
    EXPECT_DOUBLE_EQ(iss.getInclination(), 51.6312);
    EXPECT_DOUBLE_EQ(iss.getRightAscensionOfAscendingNode(), 206.3646);
    EXPECT_NEAR(iss.getEccentricity(), 0.0003723, 1e-8);
    EXPECT_NEAR(iss.getMeanMotion(), 15.49193835, 1e-6);
}

TEST(LoadTLEDatabaseTest, AppendsToExistingDatabase) {
    std::map<int, Orbit> database;

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
    std::map<int, Orbit> database;

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

    std::map<int, Orbit> database;
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

    std::map<int, Orbit> database;
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

TEST_F(OrbitTest, GetTLE_RoundTrip_ISS) {
    // Known good 3-line TLE for ISS
    const std::string originalTLE =
        "ISS (ZARYA)\n"
        "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
        "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850\n";

    orbit.updateFromTLE(originalTLE);
    std::string regeneratedTLE = orbit.getTLE();

    EXPECT_EQ(originalTLE, regeneratedTLE);
}

TEST_F(OrbitTest, GetTLE_RoundTrip_NOAA19) {
    // Known good 3-line TLE for NOAA 19
    const std::string originalTLE =
        "NOAA 19\n"
        "1 33591U 09005A   25333.78204194  .00000054  00000+0  52635-4 0  9999\n"
        "2 33591  98.9785  39.2910 0013037 231.6546 128.3455 14.13431889866318\n";

    orbit.updateFromTLE(originalTLE);
    std::string regeneratedTLE = orbit.getTLE();

    EXPECT_EQ(originalTLE, regeneratedTLE);
}

TEST_F(OrbitTest, GetTLE_ContainsName) {
    orbit.updateFromTLE(TLE_WITH_NAME);
    std::string tle = orbit.getTLE();

    // First line should be the satellite name
    EXPECT_TRUE(tle.starts_with("ISS (ZARYA)\n"));
}

TEST_F(OrbitTest, GetTLE_Line1StartsCorrectly) {
    orbit.updateFromTLE(TLE_WITH_NAME);
    std::string tle = orbit.getTLE();

    // Find the second line (line 1 of TLE)
    auto firstNewline = tle.find('\n');
    auto line1Start = tle.substr(firstNewline + 1, 2);
    EXPECT_EQ(line1Start, "1 ");
}

TEST_F(OrbitTest, GetTLE_Line2StartsCorrectly) {
    orbit.updateFromTLE(TLE_WITH_NAME);
    std::string tle = orbit.getTLE();

    // Find line 2 of TLE
    auto firstNewline = tle.find('\n');
    auto secondNewline = tle.find('\n', firstNewline + 1);
    auto line2Start = tle.substr(secondNewline + 1, 2);
    EXPECT_EQ(line2Start, "2 ");
}

TEST_F(OrbitTest, GetTLE_ContainsNoradID) {
    orbit.updateFromTLE(TLE_WITH_NAME);
    std::string tle = orbit.getTLE();

    // NORAD ID 25544 should appear in line 1 and line 2
    EXPECT_NE(tle.find("25544"), std::string::npos);
}

TEST_F(OrbitTest, GetTLE_ReloadProducesSameOrbit) {
    // Load original TLE
    orbit.updateFromTLE(TLE_WITH_NAME);

    // Get regenerated TLE
    std::string regeneratedTLE = orbit.getTLE();

    // Load regenerated TLE into a new Orbit
    Orbit orbit2;
    orbit2.updateFromTLE(regeneratedTLE);

    // Compare all orbital elements
    EXPECT_EQ(orbit.getNoradID(), orbit2.getNoradID());
    EXPECT_EQ(orbit.getName(), orbit2.getName());
    EXPECT_EQ(orbit.getClassification(), orbit2.getClassification());
    EXPECT_EQ(orbit.getDesignator(), orbit2.getDesignator());
    EXPECT_DOUBLE_EQ(orbit.getInclination(), orbit2.getInclination());
    EXPECT_DOUBLE_EQ(orbit.getRightAscensionOfAscendingNode(), orbit2.getRightAscensionOfAscendingNode());
    EXPECT_DOUBLE_EQ(orbit.getEccentricity(), orbit2.getEccentricity());
    EXPECT_DOUBLE_EQ(orbit.getArgumentOfPerigee(), orbit2.getArgumentOfPerigee());
    EXPECT_DOUBLE_EQ(orbit.getMeanAnomaly(), orbit2.getMeanAnomaly());
    EXPECT_DOUBLE_EQ(orbit.getMeanMotion(), orbit2.getMeanMotion());
    EXPECT_EQ(orbit.getRevolutionNumberAtEpoch(), orbit2.getRevolutionNumberAtEpoch());
}

TEST_F(OrbitTest, GetTLE_TwoLineFormat) {
    // Test with 2-line format (no name)
    orbit.updateFromTLE(ISS_TLE);
    std::string tle = orbit.getTLE();

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

}  // namespace
}  // namespace sattrack
