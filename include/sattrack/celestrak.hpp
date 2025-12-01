/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_CELESTRAK_HPP
#define __SATTRACK_CELESTRAK_HPP

#include <string>
#include <vector>

namespace celestrak {

/**
 * A single TLE entry containing the satellite name and two-line elements
 */
struct TLEEntry {
    std::string name;
    std::string line1;
    std::string line2;
};

/**
 * Response containing multiple TLE entries from a group query
 */
struct TLEResponse {
    std::string group;
    std::vector<TLEEntry> entries;
};

/**
 * Download TLE data for a satellite group from Celestrak
 * @param group The group name (e.g., "active", "stations", "weather", etc.)
 *              Defaults to "active" if empty
 * @param debug Print debug information if true
 * @return TLE data for all satellites in the group
 * @throws std::runtime_error if the download fails
 */
TLEResponse getTLE(const std::string& group = "active", bool debug = false);

} // namespace celestrak

#endif
