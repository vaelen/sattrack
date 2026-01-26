/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <gtest/gtest.h>
#include <sattrack/gps.hpp>

#include <chrono>
#include <cmath>
#include <thread>
#include <vector>

namespace sattrack {
namespace {

// Helper to compute NMEA checksum
uint8_t computeChecksum(const std::string& data) {
    uint8_t checksum = 0;
    for (char c : data) {
        checksum ^= static_cast<uint8_t>(c);
    }
    return checksum;
}

// Helper to create a valid NMEA sentence with checksum
std::string makeNMEA(const std::string& data) {
    uint8_t cs = computeChecksum(data);
    char buf[8];
    snprintf(buf, sizeof(buf), "*%02X", cs);
    return "$" + data + buf;
}

// Valid test sentences
const std::string VALID_GGA = makeNMEA("GPGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4,M,47.0,M,,");
const std::string VALID_RMC = makeNMEA("GPRMC,123519,A,4807.038,N,01131.000,E,022.4,084.4,230394,003.1,W");
const std::string VALID_GSA = makeNMEA("GPGSA,A,3,04,05,,09,12,,,24,,,,,2.5,1.3,2.1");

// GGA with no fix
const std::string GGA_NO_FIX = makeNMEA("GPGGA,123519,,,,,0,00,,,,,,,");

// RMC with void status
const std::string RMC_VOID = makeNMEA("GPRMC,123519,V,,,,,,,230394,,,N");

// GSA with no fix
const std::string GSA_NO_FIX = makeNMEA("GPGSA,A,1,,,,,,,,,,,,,,,");

// GSA with 2D fix
const std::string GSA_2D_FIX = makeNMEA("GPGSA,A,2,04,05,,09,12,,,24,,,,,3.0,1.5,2.6");

// Different talker IDs
const std::string GN_GGA = makeNMEA("GNGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4,M,47.0,M,,");
const std::string GL_GGA = makeNMEA("GLGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4,M,47.0,M,,");
const std::string GA_GGA = makeNMEA("GAGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4,M,47.0,M,,");

// South/West position
const std::string GGA_SOUTH_WEST = makeNMEA("GPGGA,123519,3356.789,S,11812.345,W,1,08,0.9,100.0,M,47.0,M,,");

// Tolerance for floating point comparisons
constexpr double COORD_TOLERANCE = 1e-6;
constexpr double DOP_TOLERANCE = 1e-2;

class GPSTest : public ::testing::Test {
protected:
    GPS gps;
};

// ============================================================================
// Checksum Validation Tests
// ============================================================================

TEST_F(GPSTest, ValidChecksum_ReturnsTrue) {
    EXPECT_TRUE(gps.update(VALID_GGA));
}

TEST_F(GPSTest, InvalidChecksum_ReturnsFalse) {
    // Modify the checksum to be wrong
    std::string invalid = VALID_GGA;
    invalid[invalid.size() - 1] = '0';  // Change last digit of checksum
    EXPECT_FALSE(gps.update(invalid));
}

TEST_F(GPSTest, MissingChecksum_ReturnsFalse) {
    std::string noChecksum = "$GPGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4,M,47.0,M,,";
    EXPECT_FALSE(gps.update(noChecksum));
}

TEST_F(GPSTest, MissingAsterisk_ReturnsFalse) {
    std::string noAsterisk = "$GPGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4,M,47.0,M,,47";
    EXPECT_FALSE(gps.update(noAsterisk));
}

TEST_F(GPSTest, ChecksumWithLeadingDollar) {
    EXPECT_TRUE(gps.update(VALID_GGA));
}

TEST_F(GPSTest, ChecksumWithoutLeadingDollar) {
    // Remove the leading $
    std::string noDollar = VALID_GGA.substr(1);
    EXPECT_TRUE(gps.update(noDollar));
}

// ============================================================================
// GGA Sentence Parsing Tests
// ============================================================================

TEST_F(GPSTest, GGA_ValidFix_ParsesPosition) {
    ASSERT_TRUE(gps.update(VALID_GGA));
    ASSERT_TRUE(gps.hasValidPosition());

    auto pos = gps.getPosition();
    ASSERT_TRUE(pos.has_value());

    // 48°07.038' N = 48.1173° N = 0.8401 radians (approximately)
    double expectedLat = (48.0 + 7.038 / 60.0) * DEGREES_TO_RADIANS;
    EXPECT_NEAR(pos->latInRadians, expectedLat, COORD_TOLERANCE);

    // 011°31.000' E = 11.5167° E = 0.2010 radians (approximately)
    double expectedLon = (11.0 + 31.0 / 60.0) * DEGREES_TO_RADIANS;
    EXPECT_NEAR(pos->lonInRadians, expectedLon, COORD_TOLERANCE);

    // Altitude: 545.4 M = 0.5454 km
    EXPECT_NEAR(pos->altInKilometers, 0.5454, 0.0001);
}

TEST_F(GPSTest, GGA_ValidFix_ParsesFixQuality) {
    ASSERT_TRUE(gps.update(VALID_GGA));
    EXPECT_EQ(gps.getFixQuality(), GPSFixQuality::GPSFix);
}

TEST_F(GPSTest, GGA_ValidFix_ParsesSatelliteCount) {
    ASSERT_TRUE(gps.update(VALID_GGA));
    EXPECT_EQ(gps.getSatellitesInUse(), 8);
}

TEST_F(GPSTest, GGA_ValidFix_ParsesHDOP) {
    ASSERT_TRUE(gps.update(VALID_GGA));
    auto hdop = gps.getHDOP();
    ASSERT_TRUE(hdop.has_value());
    EXPECT_NEAR(*hdop, 0.9, DOP_TOLERANCE);
}

TEST_F(GPSTest, GGA_NoFix_PositionEmpty) {
    ASSERT_TRUE(gps.update(GGA_NO_FIX));
    EXPECT_FALSE(gps.hasValidPosition());
    EXPECT_EQ(gps.getFixQuality(), GPSFixQuality::Invalid);
}

TEST_F(GPSTest, GGA_NorthLatitude_PositiveRadians) {
    ASSERT_TRUE(gps.update(VALID_GGA));
    auto pos = gps.getPosition();
    ASSERT_TRUE(pos.has_value());
    EXPECT_GT(pos->latInRadians, 0.0);
}

TEST_F(GPSTest, GGA_SouthLatitude_NegativeRadians) {
    ASSERT_TRUE(gps.update(GGA_SOUTH_WEST));
    auto pos = gps.getPosition();
    ASSERT_TRUE(pos.has_value());
    EXPECT_LT(pos->latInRadians, 0.0);
}

TEST_F(GPSTest, GGA_EastLongitude_PositiveRadians) {
    ASSERT_TRUE(gps.update(VALID_GGA));
    auto pos = gps.getPosition();
    ASSERT_TRUE(pos.has_value());
    EXPECT_GT(pos->lonInRadians, 0.0);
}

TEST_F(GPSTest, GGA_WestLongitude_NegativeRadians) {
    ASSERT_TRUE(gps.update(GGA_SOUTH_WEST));
    auto pos = gps.getPosition();
    ASSERT_TRUE(pos.has_value());
    EXPECT_LT(pos->lonInRadians, 0.0);
}

TEST_F(GPSTest, GGA_AltitudeConversion_MetersToKm) {
    ASSERT_TRUE(gps.update(GGA_SOUTH_WEST));
    auto pos = gps.getPosition();
    ASSERT_TRUE(pos.has_value());
    // 100.0 M = 0.1 km
    EXPECT_NEAR(pos->altInKilometers, 0.1, 0.0001);
}

TEST_F(GPSTest, GGA_EmptyFields_GracefulHandling) {
    // GGA with empty optional fields should still parse
    std::string ggaEmpty = makeNMEA("GPGGA,123519,4807.038,N,01131.000,E,1,08,,,M,,M,,");
    ASSERT_TRUE(gps.update(ggaEmpty));
    EXPECT_TRUE(gps.hasValidPosition());
    EXPECT_FALSE(gps.getHDOP().has_value());
}

// ============================================================================
// RMC Sentence Parsing Tests
// ============================================================================

TEST_F(GPSTest, RMC_Active_ParsesTimeAndDate) {
    ASSERT_TRUE(gps.update(VALID_RMC));
    ASSERT_TRUE(gps.hasValidTime());

    auto utc = gps.getUTCTime();
    ASSERT_TRUE(utc.has_value());

    // Date: 23/03/1994, Time: 12:35:19 UTC
    // We can't easily verify the exact time_point, but we can check it's set
}

TEST_F(GPSTest, RMC_Active_ParsesSpeed) {
    ASSERT_TRUE(gps.update(VALID_RMC));
    auto speed = gps.getSpeedKnots();
    ASSERT_TRUE(speed.has_value());
    EXPECT_NEAR(*speed, 22.4, 0.1);
}

TEST_F(GPSTest, RMC_Active_ParsesCourse) {
    ASSERT_TRUE(gps.update(VALID_RMC));
    auto course = gps.getCourseDegrees();
    ASSERT_TRUE(course.has_value());
    EXPECT_NEAR(*course, 84.4, 0.1);
}

TEST_F(GPSTest, RMC_Void_NoTimeUpdate) {
    // First set valid time
    ASSERT_TRUE(gps.update(VALID_RMC));
    ASSERT_TRUE(gps.hasValidTime());

    // Reset and parse void RMC
    gps.reset();
    ASSERT_TRUE(gps.update(RMC_VOID));
    // Time should not be updated when status is void
    EXPECT_FALSE(gps.hasValidTime());
}

TEST_F(GPSTest, RMC_DataStatus_Active) {
    ASSERT_TRUE(gps.update(VALID_RMC));
    EXPECT_EQ(gps.getDataStatus(), GPSDataStatus::Active);
}

TEST_F(GPSTest, RMC_DataStatus_Void) {
    ASSERT_TRUE(gps.update(RMC_VOID));
    EXPECT_EQ(gps.getDataStatus(), GPSDataStatus::Void);
}

// ============================================================================
// GSA Sentence Parsing Tests
// ============================================================================

TEST_F(GPSTest, GSA_NoFix_FixTypeIsNoFix) {
    ASSERT_TRUE(gps.update(GSA_NO_FIX));
    EXPECT_EQ(gps.getFixType(), GPSFixType::NoFix);
    EXPECT_FALSE(gps.hasFix());
}

TEST_F(GPSTest, GSA_2DFix_FixTypeIs2D) {
    ASSERT_TRUE(gps.update(GSA_2D_FIX));
    EXPECT_EQ(gps.getFixType(), GPSFixType::Fix2D);
    EXPECT_TRUE(gps.hasFix());
    EXPECT_FALSE(gps.has3DFix());
}

TEST_F(GPSTest, GSA_3DFix_FixTypeIs3D) {
    ASSERT_TRUE(gps.update(VALID_GSA));
    EXPECT_EQ(gps.getFixType(), GPSFixType::Fix3D);
    EXPECT_TRUE(gps.hasFix());
    EXPECT_TRUE(gps.has3DFix());
}

TEST_F(GPSTest, GSA_ParsesPDOP) {
    ASSERT_TRUE(gps.update(VALID_GSA));
    auto pdop = gps.getPDOP();
    ASSERT_TRUE(pdop.has_value());
    EXPECT_NEAR(*pdop, 2.5, DOP_TOLERANCE);
}

TEST_F(GPSTest, GSA_ParsesHDOP) {
    ASSERT_TRUE(gps.update(VALID_GSA));
    auto hdop = gps.getHDOP();
    ASSERT_TRUE(hdop.has_value());
    EXPECT_NEAR(*hdop, 1.3, DOP_TOLERANCE);
}

TEST_F(GPSTest, GSA_ParsesVDOP) {
    ASSERT_TRUE(gps.update(VALID_GSA));
    auto vdop = gps.getVDOP();
    ASSERT_TRUE(vdop.has_value());
    EXPECT_NEAR(*vdop, 2.1, DOP_TOLERANCE);
}

// ============================================================================
// State Query Tests
// ============================================================================

TEST_F(GPSTest, HasValidPosition_NoData_ReturnsFalse) {
    EXPECT_FALSE(gps.hasValidPosition());
}

TEST_F(GPSTest, HasValidPosition_AfterGGA_ReturnsTrue) {
    gps.update(VALID_GGA);
    EXPECT_TRUE(gps.hasValidPosition());
}

TEST_F(GPSTest, HasValidTime_NoData_ReturnsFalse) {
    EXPECT_FALSE(gps.hasValidTime());
}

TEST_F(GPSTest, HasValidTime_AfterRMC_ReturnsTrue) {
    gps.update(VALID_RMC);
    EXPECT_TRUE(gps.hasValidTime());
}

TEST_F(GPSTest, HasFix_NoFix_ReturnsFalse) {
    gps.update(GSA_NO_FIX);
    EXPECT_FALSE(gps.hasFix());
}

TEST_F(GPSTest, HasFix_2DFix_ReturnsTrue) {
    gps.update(GSA_2D_FIX);
    EXPECT_TRUE(gps.hasFix());
}

TEST_F(GPSTest, Has3DFix_2DFix_ReturnsFalse) {
    gps.update(GSA_2D_FIX);
    EXPECT_FALSE(gps.has3DFix());
}

TEST_F(GPSTest, Has3DFix_3DFix_ReturnsTrue) {
    gps.update(VALID_GSA);
    EXPECT_TRUE(gps.has3DFix());
}

// ============================================================================
// Coordinate Conversion Tests
// ============================================================================

TEST_F(GPSTest, LatitudeConversion_DDMM_ToRadians) {
    ASSERT_TRUE(gps.update(VALID_GGA));
    auto pos = gps.getPosition();
    ASSERT_TRUE(pos.has_value());

    // 4807.038 = 48 degrees, 7.038 minutes = 48.1173 degrees
    double expectedDegrees = 48.0 + 7.038 / 60.0;
    double expectedRadians = expectedDegrees * DEGREES_TO_RADIANS;
    EXPECT_NEAR(pos->latInRadians, expectedRadians, COORD_TOLERANCE);
}

TEST_F(GPSTest, LongitudeConversion_DDDMM_ToRadians) {
    ASSERT_TRUE(gps.update(VALID_GGA));
    auto pos = gps.getPosition();
    ASSERT_TRUE(pos.has_value());

    // 01131.000 = 11 degrees, 31.000 minutes = 11.5167 degrees
    double expectedDegrees = 11.0 + 31.0 / 60.0;
    double expectedRadians = expectedDegrees * DEGREES_TO_RADIANS;
    EXPECT_NEAR(pos->lonInRadians, expectedRadians, COORD_TOLERANCE);
}

TEST_F(GPSTest, CoordinateAccuracy_SubMeterPrecision) {
    // High precision coordinates
    std::string preciseGGA = makeNMEA("GPGGA,123519,4807.03800,N,01131.00000,E,1,08,0.9,545.4,M,47.0,M,,");
    ASSERT_TRUE(gps.update(preciseGGA));
    auto pos = gps.getPosition();
    ASSERT_TRUE(pos.has_value());

    // Verify we maintain precision
    double expectedLat = (48.0 + 7.038 / 60.0) * DEGREES_TO_RADIANS;
    EXPECT_NEAR(pos->latInRadians, expectedLat, 1e-8);
}

// ============================================================================
// Talker ID Tests
// ============================================================================

TEST_F(GPSTest, GP_TalkerID_Accepted) {
    EXPECT_TRUE(gps.update(VALID_GGA));
    EXPECT_TRUE(gps.hasValidPosition());
}

TEST_F(GPSTest, GN_TalkerID_Accepted) {
    EXPECT_TRUE(gps.update(GN_GGA));
    EXPECT_TRUE(gps.hasValidPosition());
}

TEST_F(GPSTest, GL_TalkerID_Accepted) {
    EXPECT_TRUE(gps.update(GL_GGA));
    EXPECT_TRUE(gps.hasValidPosition());
}

TEST_F(GPSTest, GA_TalkerID_Accepted) {
    EXPECT_TRUE(gps.update(GA_GGA));
    EXPECT_TRUE(gps.hasValidPosition());
}

TEST_F(GPSTest, UnknownSentence_Ignored) {
    std::string vtg = makeNMEA("GPVTG,054.7,T,034.4,M,005.5,N,010.2,K");
    // Should return false (not processed) but not crash
    EXPECT_FALSE(gps.update(vtg));
}

// ============================================================================
// Edge Case Tests
// ============================================================================

TEST_F(GPSTest, EmptySentence_ReturnsFalse) {
    EXPECT_FALSE(gps.update(""));
}

TEST_F(GPSTest, WhitespaceOnly_ReturnsFalse) {
    EXPECT_FALSE(gps.update("   "));
}

TEST_F(GPSTest, TruncatedSentence_ReturnsFalse) {
    EXPECT_FALSE(gps.update("$GPGGA,123"));
}

TEST_F(GPSTest, ExtraCommas_GracefulHandling) {
    // Extra trailing commas should be handled
    std::string extraCommas = makeNMEA("GPGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4,M,47.0,M,,,,,");
    // May or may not parse depending on field count, but shouldn't crash
    gps.update(extraCommas);
}

TEST_F(GPSTest, VeryLongSentence_NoOverflow) {
    // Create an excessively long sentence
    std::string longData = "GPGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4,M,47.0,M,,";
    for (int i = 0; i < 100; i++) {
        longData += ",extra";
    }
    std::string longSentence = makeNMEA(longData);
    // Should handle gracefully without crashing
    gps.update(longSentence);
}

// ============================================================================
// Thread Safety Tests
// ============================================================================

TEST_F(GPSTest, ConcurrentUpdates_NoDataRace) {
    std::vector<std::thread> threads;

    for (int i = 0; i < 10; i++) {
        threads.emplace_back([this]() {
            for (int j = 0; j < 100; j++) {
                gps.update(VALID_GGA);
                gps.update(VALID_RMC);
                gps.update(VALID_GSA);
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // If we get here without crashing, the test passes
    EXPECT_TRUE(gps.hasValidPosition());
}

TEST_F(GPSTest, ConcurrentReads_NoDataRace) {
    // Set up initial state
    gps.update(VALID_GGA);
    gps.update(VALID_RMC);
    gps.update(VALID_GSA);

    std::vector<std::thread> threads;

    for (int i = 0; i < 10; i++) {
        threads.emplace_back([this]() {
            for (int j = 0; j < 100; j++) {
                gps.hasValidPosition();
                gps.getPosition();
                gps.hasValidTime();
                gps.getUTCTime();
                gps.getFixQuality();
                gps.getFixType();
                gps.hasFix();
                gps.has3DFix();
                gps.getSatellitesInUse();
                gps.getHDOP();
                gps.getVDOP();
                gps.getPDOP();
                gps.getSpeedKnots();
                gps.getCourseDegrees();
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // If we get here without crashing, the test passes
    SUCCEED();
}

// ============================================================================
// Reset Tests
// ============================================================================

TEST_F(GPSTest, Reset_ClearsPosition) {
    gps.update(VALID_GGA);
    ASSERT_TRUE(gps.hasValidPosition());

    gps.reset();
    EXPECT_FALSE(gps.hasValidPosition());
}

TEST_F(GPSTest, Reset_ClearsTime) {
    gps.update(VALID_RMC);
    ASSERT_TRUE(gps.hasValidTime());

    gps.reset();
    EXPECT_FALSE(gps.hasValidTime());
}

TEST_F(GPSTest, Reset_ResetsFixQuality) {
    gps.update(VALID_GGA);
    ASSERT_EQ(gps.getFixQuality(), GPSFixQuality::GPSFix);

    gps.reset();
    EXPECT_EQ(gps.getFixQuality(), GPSFixQuality::Invalid);
}

TEST_F(GPSTest, Reset_ResetsFixType) {
    gps.update(VALID_GSA);
    ASSERT_EQ(gps.getFixType(), GPSFixType::Fix3D);

    gps.reset();
    EXPECT_EQ(gps.getFixType(), GPSFixType::NoFix);
}

// ============================================================================
// Integration Tests
// ============================================================================

TEST_F(GPSTest, FullSequence_AllDataPopulated) {
    // Simulate receiving a full sequence of NMEA sentences
    ASSERT_TRUE(gps.update(VALID_GGA));
    ASSERT_TRUE(gps.update(VALID_RMC));
    ASSERT_TRUE(gps.update(VALID_GSA));

    // Verify all data is populated
    EXPECT_TRUE(gps.hasValidPosition());
    EXPECT_TRUE(gps.hasValidTime());
    EXPECT_TRUE(gps.has3DFix());
    EXPECT_EQ(gps.getDataStatus(), GPSDataStatus::Active);
    EXPECT_TRUE(gps.getHDOP().has_value());
    EXPECT_TRUE(gps.getVDOP().has_value());
    EXPECT_TRUE(gps.getPDOP().has_value());
    EXPECT_TRUE(gps.getSpeedKnots().has_value());
    EXPECT_TRUE(gps.getCourseDegrees().has_value());
}

} // namespace
} // namespace sattrack
