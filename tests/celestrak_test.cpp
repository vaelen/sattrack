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

// ============================================================================
// TLEEntry Tests
// ============================================================================

TEST(TLEEntryTest, DefaultConstruction) {
    TLEEntry entry;
    EXPECT_TRUE(entry.name.empty());
    EXPECT_TRUE(entry.line1.empty());
    EXPECT_TRUE(entry.line2.empty());
}

TEST(TLEEntryTest, InitializerConstruction) {
    TLEEntry entry{
        .name = "ISS (ZARYA)",
        .line1 = "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993",
        .line2 = "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850"
    };
    EXPECT_EQ(entry.name, "ISS (ZARYA)");
    EXPECT_EQ(entry.line1[0], '1');
    EXPECT_EQ(entry.line2[0], '2');
}

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
        .entries = {
            {
                .name = "ISS (ZARYA)",
                .line1 = "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993",
                .line2 = "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850"
            }
        }
    };
    EXPECT_EQ(response.group, "stations");
    EXPECT_EQ(response.entries.size(), 1);
    EXPECT_EQ(response.entries[0].name, "ISS (ZARYA)");
}

TEST(TLEResponseTest, MultipleEntries) {
    TLEResponse response{
        .group = "active",
        .entries = {
            {
                .name = "ISS (ZARYA)",
                .line1 = "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993",
                .line2 = "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850"
            },
            {
                .name = "NOAA 19",
                .line1 = "1 33591U 09005A   25333.78204194  .00000054  00000+0  52635-4 0  9999",
                .line2 = "2 33591  98.9785  39.2910 0013037 231.6546 128.3455 14.13431889866318"
            }
        }
    };
    EXPECT_EQ(response.entries.size(), 2);
    EXPECT_EQ(response.entries[0].name, "ISS (ZARYA)");
    EXPECT_EQ(response.entries[1].name, "NOAA 19");
}

// ============================================================================
// TLE Line Validation Tests
// ============================================================================

TEST(TLEEntryTest, Line1StartsWithOne) {
    TLEEntry entry{
        .name = "ISS (ZARYA)",
        .line1 = "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993",
        .line2 = "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850"
    };
    EXPECT_TRUE(entry.line1.length() >= 2);
    EXPECT_EQ(entry.line1[0], '1');
    EXPECT_EQ(entry.line1[1], ' ');
}

TEST(TLEEntryTest, Line2StartsWithTwo) {
    TLEEntry entry{
        .name = "ISS (ZARYA)",
        .line1 = "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993",
        .line2 = "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850"
    };
    EXPECT_TRUE(entry.line2.length() >= 2);
    EXPECT_EQ(entry.line2[0], '2');
    EXPECT_EQ(entry.line2[1], ' ');
}

TEST(TLEEntryTest, TLELineLength) {
    // Standard TLE lines are 69 characters
    TLEEntry entry{
        .name = "ISS (ZARYA)",
        .line1 = "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993",
        .line2 = "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850"
    };
    EXPECT_EQ(entry.line1.length(), 69);
    EXPECT_EQ(entry.line2.length(), 69);
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
        .entries = {
            {.name = "ISS", .line1 = "1 ...", .line2 = "2 ..."},
            {.name = "CSS", .line1 = "1 ...", .line2 = "2 ..."}
        }
    };
    EXPECT_EQ(response.entries[0].name, "ISS");
    EXPECT_EQ(response.entries[1].name, "CSS");
}

TEST(TLEResponseTest, IterateEntries) {
    TLEResponse response{
        .group = "stations",
        .entries = {
            {.name = "SAT1", .line1 = "1 ...", .line2 = "2 ..."},
            {.name = "SAT2", .line1 = "1 ...", .line2 = "2 ..."},
            {.name = "SAT3", .line1 = "1 ...", .line2 = "2 ..."}
        }
    };

    int count = 0;
    for (const auto& entry : response.entries) {
        EXPECT_FALSE(entry.name.empty());
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
// NORAD ID Extraction Tests (from TLE line 1)
// ============================================================================

TEST(TLEEntryTest, ExtractNoradIDFromLine1) {
    TLEEntry entry{
        .name = "ISS (ZARYA)",
        .line1 = "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993",
        .line2 = "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850"
    };
    // NORAD ID is at columns 3-7 (0-indexed: 2-6)
    std::string noradIdStr = entry.line1.substr(2, 5);
    int noradId = std::stoi(noradIdStr);
    EXPECT_EQ(noradId, 25544);
}

TEST(TLEEntryTest, ExtractNoradIDFromLine2) {
    TLEEntry entry{
        .name = "ISS (ZARYA)",
        .line1 = "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993",
        .line2 = "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850"
    };
    // NORAD ID is also at columns 3-7 in line 2
    std::string noradIdStr = entry.line2.substr(2, 5);
    int noradId = std::stoi(noradIdStr);
    EXPECT_EQ(noradId, 25544);
}

TEST(TLEEntryTest, NoradIDConsistentBetweenLines) {
    TLEEntry entry{
        .name = "ISS (ZARYA)",
        .line1 = "1 25544U 98067A   25333.83453771  .00008010  00000+0  15237-3 0  9993",
        .line2 = "2 25544  51.6312 206.3646 0003723 184.1118 175.9840 15.49193835540850"
    };
    std::string id1 = entry.line1.substr(2, 5);
    std::string id2 = entry.line2.substr(2, 5);
    EXPECT_EQ(id1, id2);
}

}  // namespace
}  // namespace celestrak
