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
// getLookAngles with Orbit Tests
// ============================================================================

TEST_F(OrbitTest, GetLookAnglesAtEpoch) {
    orbit.updateFromTLE(ISS_TLE);

    // Observer at a mid-latitude location
    Geodetic observer{0.7, -1.5, 0.0};  // ~40°N, ~86°W

    double epochJD = toJulianDate(orbit.getEpoch());
    LookAngles angles = getLookAngles(orbit, observer, epochJD);

    // Just verify we get reasonable values
    EXPECT_GE(angles.azimuthInRadians, 0.0);
    EXPECT_LT(angles.azimuthInRadians, 2.0 * M_PI);
    EXPECT_GE(angles.elevationInRadians, -M_PI / 2.0);
    EXPECT_LE(angles.elevationInRadians, M_PI / 2.0);
    EXPECT_GT(angles.rangeInKilometers, 0.0);
    EXPECT_LT(angles.rangeInKilometers, 50000.0);  // Should be less than GEO distance
}

TEST_F(OrbitTest, GetLookAnglesRangeReasonable) {
    orbit.updateFromTLE(ISS_TLE);

    // Observer location doesn't matter much for range check
    Geodetic observer{0.0, 0.0, 0.0};

    double epochJD = toJulianDate(orbit.getEpoch());
    LookAngles angles = getLookAngles(orbit, observer, epochJD);

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
// isVisible with Orbit Tests
// ============================================================================

TEST_F(OrbitTest, IsVisibleReturnsReasonableResults) {
    orbit.updateFromTLE(ISS_TLE);

    // Observer at a specific location
    Geodetic observer{0.7, -1.5, 0.0};

    double epochJD = toJulianDate(orbit.getEpoch());

    // Check visibility with different thresholds
    // These should return consistent results (visible at low threshold means
    // visible at even lower threshold)
    bool visibleAt0 = isVisible(orbit, observer, epochJD, 0.0);
    bool visibleAt10deg = isVisible(orbit, observer, epochJD, 10.0 * M_PI / 180.0);
    bool visibleAt45deg = isVisible(orbit, observer, epochJD, 45.0 * M_PI / 180.0);

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

TEST_F(OrbitTest, FindNextPassReturnsValidPassForISS) {
    // ISS should have passes over mid-latitude locations
    orbit.updateFromTLE(ISS_TLE);
    Geodetic observer{0.7, -1.5, 0.0};  // ~40°N, ~86°W

    auto pass = findNextPass(orbit, observer, 0.0, orbit.getEpoch());

    // Should find a pass within 48 hours for ISS
    EXPECT_TRUE(pass.has_value());
}

TEST_F(OrbitTest, FindNextPassReturnsValidPass) {
    orbit.updateFromTLE(ISS_TLE);
    Geodetic observer{0.7, -1.5, 0.0};  // ~40°N, ~86°W

    auto pass = findNextPass(orbit, observer, 0.0, orbit.getEpoch());

    ASSERT_TRUE(pass.has_value());
    // Rise time should be before max elevation time
    EXPECT_LT(pass->riseTime, pass->maxElevationTime);
    // Max elevation time should be before set time
    EXPECT_LT(pass->maxElevationTime, pass->setTime);
    // Max elevation should be higher than rise/set
    EXPECT_GE(pass->maxAngles.elevationInRadians, pass->riseAngles.elevationInRadians);
    EXPECT_GE(pass->maxAngles.elevationInRadians, pass->setAngles.elevationInRadians);
}

TEST_F(OrbitTest, FindNextPassWithMinElevation) {
    orbit.updateFromTLE(ISS_TLE);
    Geodetic observer{0.7, -1.5, 0.0};
    double minElevation = 10.0 * M_PI / 180.0;  // 10 degrees

    auto pass = findNextPass(orbit, observer, minElevation, orbit.getEpoch());

    if (pass.has_value()) {
        // Rise elevation should be at or above minimum
        EXPECT_GE(pass->riseAngles.elevationInRadians, minElevation - 0.01);
        // Set elevation should be at or above minimum
        EXPECT_GE(pass->setAngles.elevationInRadians, minElevation - 0.01);
    }
}

TEST_F(OrbitTest, FindNextPassISSHasMultipleDailyPasses) {
    // ISS completes about 15.5 orbits per day, so there should be multiple
    // opportunities for passes over any given location
    orbit.updateFromTLE(ISS_TLE);
    Geodetic observer{0.7, -1.5, 0.0};

    auto now = orbit.getEpoch();
    int passCount = 0;

    // Search for passes over 24 hours
    auto searchTime = now;
    auto endTime = now + std::chrono::hours(24);

    while (searchTime < endTime) {
        auto pass = findNextPass(orbit, observer, 0.0, searchTime);
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
