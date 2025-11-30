/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_N2YO_HPP
#define __SATTRACK_N2YO_HPP

#include <string>
#include <vector>

namespace n2yo {

struct SatelliteInfo {
    int id;
    std::string name;
};

struct TLEResponse {
    SatelliteInfo info;
    std::string tle;
};

struct Position {
    double latitude;
    double longitude;
    double altitude;
    double azimuth;
    double elevation;
    double ra;           // right ascension
    double dec;          // declination
    long timestamp;
};

struct PositionsResponse {
    SatelliteInfo info;
    std::vector<Position> positions;
};

struct RadioPass {
    double startAz;
    std::string startAzCompass;
    long startUTC;
    double maxAz;
    std::string maxAzCompass;
    double maxEl;
    long maxUTC;
    double endAz;
    std::string endAzCompass;
    long endUTC;
};

struct RadioPassesResponse {
    SatelliteInfo info;
    std::vector<RadioPass> passes;
};

/**
 * Get TLE data for a satellite from the N2YO API
 * @param apiKey The N2YO API key
 * @param noradId The NORAD catalog ID of the satellite
 * @return The TLE data as a string (two lines separated by \r\n)
 * @throws std::runtime_error if the API call fails
 */
TLEResponse getTLE(const std::string& apiKey, int noradId, bool debug = false);

/**
 * Get satellite positions
 * @param apiKey The N2YO API key
 * @param noradId The NORAD catalog ID of the satellite
 * @param observerLat Observer's latitude in decimal degrees
 * @param observerLng Observer's longitude in decimal degrees
 * @param observerAlt Observer's altitude above sea level in meters
 * @param seconds Number of future positions to retrieve (max 300)
 * @return Position data for the satellite
 * @throws std::runtime_error if the API call fails
 */
PositionsResponse getPositions(const std::string& apiKey, int noradId,
                               double observerLat, double observerLng,
                               double observerAlt, int seconds,
                               bool debug = false);

/**
 * Get radio passes for a satellite
 * @param apiKey The N2YO API key
 * @param noradId The NORAD catalog ID of the satellite
 * @param observerLat Observer's latitude in decimal degrees
 * @param observerLng Observer's longitude in decimal degrees
 * @param observerAlt Observer's altitude above sea level in meters
 * @param days Number of days to predict (max 10)
 * @param minElevation Minimum elevation threshold in degrees
 * @return Radio pass predictions
 * @throws std::runtime_error if the API call fails
 */
RadioPassesResponse getRadioPasses(const std::string& apiKey, int noradId,
                                   double observerLat, double observerLng,
                                   double observerAlt, int days, int minElevation,
                                   bool debug = false);

} // namespace n2yo

#endif