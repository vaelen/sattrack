/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <gtest/gtest.h>
#include <sattrack/celestrak.hpp>

#include <string>
#include <vector>

namespace celestrak {
namespace {

// Sample 3-line TLE strings for testing
const std::string ISS_TLE =
    "ISS (ZARYA)\n"
    "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993\n"
    "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850\n";

const std::string NOAA19_TLE =
    "NOAA 19\n"
    "1 33591U 09005A   25333.78204194  .00000054  00000+0  52635-4 0  9999\n"
    "2 33591  98.9785  39.2910 0013037 231.6546 128.3455 14.13431889866318\n";

// ============================================================================
// TLEResponse Tests
// ============================================================================

TEST(TLEResponseTest, DefaultConstruction) {
    TLEResponse response;
    EXPECT_TRUE(response.group.empty());
    EXPECT_TRUE(response.entries.empty());
}

TEST(TLEResponseTest, InitializerConstruction) {
    TLEResponse response{
        .group = "stations",
        .entries = { ISS_TLE }
    };
    EXPECT_EQ(response.group, "stations");
    EXPECT_EQ(response.entries.size(), 1);
    EXPECT_TRUE(response.entries[0].find("ISS (ZARYA)") != std::string::npos);
}

TEST(TLEResponseTest, MultipleEntries) {
    TLEResponse response{
        .group = "active",
        .entries = { ISS_TLE, NOAA19_TLE }
    };
    EXPECT_EQ(response.entries.size(), 2);
    EXPECT_TRUE(response.entries[0].find("ISS") != std::string::npos);
    EXPECT_TRUE(response.entries[1].find("NOAA 19") != std::string::npos);
}

// ============================================================================
// Group Name Tests
// ============================================================================

TEST(TLEResponseTest, GroupNameStations) {
    TLEResponse response{.group = "stations", .entries = {}};
    EXPECT_EQ(response.group, "stations");
}

TEST(TLEResponseTest, GroupNameActive) {
    TLEResponse response{.group = "active", .entries = {}};
    EXPECT_EQ(response.group, "active");
}

TEST(TLEResponseTest, GroupNameWeather) {
    TLEResponse response{.group = "weather", .entries = {}};
    EXPECT_EQ(response.group, "weather");
}

TEST(TLEResponseTest, GroupNameLastDays) {
    TLEResponse response{.group = "last-30-days", .entries = {}};
    EXPECT_EQ(response.group, "last-30-days");
}

// ============================================================================
// Entry Access Tests
// ============================================================================

TEST(TLEResponseTest, AccessEntryByIndex) {
    TLEResponse response{
        .group = "stations",
        .entries = { ISS_TLE, NOAA19_TLE }
    };
    EXPECT_TRUE(response.entries[0].find("ISS") != std::string::npos);
    EXPECT_TRUE(response.entries[1].find("NOAA") != std::string::npos);
}

TEST(TLEResponseTest, IterateEntries) {
    TLEResponse response{
        .group = "stations",
        .entries = { ISS_TLE, NOAA19_TLE, ISS_TLE }
    };

    int count = 0;
    for (const auto& entry : response.entries) {
        EXPECT_FALSE(entry.empty());
        count++;
    }
    EXPECT_EQ(count, 3);
}

TEST(TLEResponseTest, EmptyResponse) {
    TLEResponse response{.group = "empty-group", .entries = {}};
    EXPECT_TRUE(response.entries.empty());
    EXPECT_EQ(response.entries.size(), 0);
}

// ============================================================================
// TLE Entry Format Tests
// ============================================================================

TEST(TLEEntryFormatTest, EntryContainsThreeLines) {
    // Each entry should be a 3-line TLE string
    std::string entry = ISS_TLE;

    // Count newlines
    int newlineCount = 0;
    for (char c : entry) {
        if (c == '\n') newlineCount++;
    }
    EXPECT_EQ(newlineCount, 3);
}

TEST(TLEEntryFormatTest, EntryStartsWithName) {
    std::string entry = ISS_TLE;
    EXPECT_TRUE(entry.find("ISS (ZARYA)") == 0);
}

TEST(TLEEntryFormatTest, EntryContainsLine1) {
    std::string entry = ISS_TLE;
    EXPECT_TRUE(entry.find("\n1 ") != std::string::npos);
}

TEST(TLEEntryFormatTest, EntryContainsLine2) {
    std::string entry = ISS_TLE;
    EXPECT_TRUE(entry.find("\n2 ") != std::string::npos);
}

// ============================================================================
// NORAD ID Extraction Tests (from TLE entry string)
// ============================================================================

TEST(TLEEntryFormatTest, ExtractNoradIDFromEntry) {
    std::string entry = ISS_TLE;

    // Find line 1 (starts after first newline with "1 ")
    auto line1Start = entry.find("\n1 ");
    ASSERT_NE(line1Start, std::string::npos);

    // NORAD ID is at columns 3-7 of line 1 (0-indexed: 2-6 from line start)
    std::string noradIdStr = entry.substr(line1Start + 3, 5);
    int noradId = std::stoi(noradIdStr);
    EXPECT_EQ(noradId, 25544);
}

TEST(TLEEntryFormatTest, ExtractNoradIDFromLine2) {
    std::string entry = ISS_TLE;

    // Find line 2 (starts after "1 " line with "2 ")
    auto line2Start = entry.find("\n2 ");
    ASSERT_NE(line2Start, std::string::npos);

    // NORAD ID is at columns 3-7 of line 2
    std::string noradIdStr = entry.substr(line2Start + 3, 5);
    int noradId = std::stoi(noradIdStr);
    EXPECT_EQ(noradId, 25544);
}

TEST(TLEEntryFormatTest, NoradIDConsistentBetweenLines) {
    std::string entry = ISS_TLE;

    auto line1Start = entry.find("\n1 ");
    auto line2Start = entry.find("\n2 ");
    ASSERT_NE(line1Start, std::string::npos);
    ASSERT_NE(line2Start, std::string::npos);

    std::string id1 = entry.substr(line1Start + 3, 5);
    std::string id2 = entry.substr(line2Start + 3, 5);
    EXPECT_EQ(id1, id2);
}

}  // namespace
}  // namespace celestrak
