/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_GPS_HPP
#define __SATTRACK_GPS_HPP

#include <sattrack/satellite.hpp>

#include <chrono>
#include <cstdint>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace sattrack {

/**
 * GPS fix quality from GGA sentence (field 6).
 */
enum class GPSFixQuality : uint8_t {
    Invalid = 0,
    GPSFix = 1,
    DGPSFix = 2,
    PPSFix = 3,
    RTKFixed = 4,
    RTKFloat = 5,
    Estimated = 6,
    Manual = 7,
    Simulation = 8
};

/**
 * Fix type from GSA sentence (field 2).
 */
enum class GPSFixType : uint8_t {
    NoFix = 1,
    Fix2D = 2,
    Fix3D = 3
};

/**
 * GPS data status from RMC sentence (field 2).
 */
enum class GPSDataStatus : char {
    Active = 'A',
    Void = 'V'
};

/**
 * GPS NMEA parser class.
 *
 * Parses NMEA sentences (GGA, RMC, GSA) from a GPS receiver and maintains
 * current position, time, and fix state. Thread-safe for concurrent access.
 *
 * Usage:
 *   GPS gps;
 *   gps.update("$GPGGA,123456.00,4807.038,N,01131.000,E,1,08,0.9,545.4,M,47.0,M,,*47");
 *   if (gps.hasValidPosition()) {
 *       auto pos = gps.getPosition();
 *   }
 */
class GPS {
public:
    GPS() = default;
    ~GPS() = default;

    // Non-copyable, non-movable (due to mutex)
    GPS(const GPS&) = delete;
    GPS& operator=(const GPS&) = delete;
    GPS(GPS&&) = delete;
    GPS& operator=(GPS&&) = delete;

    /**
     * Process an NMEA sentence and update internal state.
     *
     * Handles GGA, RMC, and GSA sentences. Other sentence types are ignored.
     * Invalid sentences (bad checksum, malformed) are silently ignored.
     *
     * @param sentence Complete NMEA sentence (with or without leading $)
     * @return true if the sentence was valid and processed, false otherwise
     */
    bool update(std::string_view sentence);

    // === Position getters ===

    /** Check if a valid position has been received */
    bool hasValidPosition() const;

    /** Get current position (latitude, longitude, altitude) */
    std::optional<Geodetic> getPosition() const;

    // === Time getters ===

    /** Check if valid time has been received */
    bool hasValidTime() const;

    /** Get current UTC time from GPS */
    std::optional<time_point> getUTCTime() const;

    // === Fix status getters ===

    /** Get fix quality (from GGA sentence) */
    GPSFixQuality getFixQuality() const;

    /** Get fix type: NoFix, 2D, or 3D (from GSA sentence) */
    GPSFixType getFixType() const;

    /** Check if we have a valid fix (at least 2D) */
    bool hasFix() const;

    /** Check if we have a 3D fix */
    bool has3DFix() const;

    /** Get number of satellites in use */
    int getSatellitesInUse() const;

    /** Get data status (Active/Void from RMC) */
    GPSDataStatus getDataStatus() const;

    // === Dilution of Precision getters ===

    std::optional<double> getHDOP() const;
    std::optional<double> getVDOP() const;
    std::optional<double> getPDOP() const;

    // === Speed and course getters ===

    /** Get speed over ground in knots */
    std::optional<double> getSpeedKnots() const;

    /** Get course over ground in degrees (true north) */
    std::optional<double> getCourseDegrees() const;

    // === Utility methods ===

    /** Reset all state to defaults */
    void reset();

private:
    mutable std::mutex mutex_;

    // Position data (from GGA)
    std::optional<Geodetic> position_;

    // Time (from RMC)
    std::optional<time_point> utcTime_;

    // Fix information (from GGA)
    GPSFixQuality fixQuality_ = GPSFixQuality::Invalid;
    int satellitesInUse_ = 0;
    std::optional<double> hdop_;

    // Fix type (from GSA)
    GPSFixType fixType_ = GPSFixType::NoFix;
    std::optional<double> pdop_;
    std::optional<double> vdop_;

    // Speed and course (from RMC)
    std::optional<double> speedOverGroundKnots_;
    std::optional<double> courseOverGroundDegrees_;

    // Data status (from RMC)
    GPSDataStatus dataStatus_ = GPSDataStatus::Void;

    // NMEA parsing helpers
    bool validateChecksum(std::string_view sentence) const;
    static uint8_t calculateChecksum(std::string_view data);
    static std::vector<std::string_view> splitFields(std::string_view sentence);

    // Sentence parsers (return true on success)
    bool parseGGA(const std::vector<std::string_view>& fields);
    bool parseRMC(const std::vector<std::string_view>& fields);
    bool parseGSA(const std::vector<std::string_view>& fields);

    // Field parsing helpers
    static std::optional<double> parseLatitude(std::string_view field, char hemisphere);
    static std::optional<double> parseLongitude(std::string_view field, char hemisphere);
    static std::optional<double> parseDouble(std::string_view field);
    static std::optional<int> parseInt(std::string_view field);
    static std::optional<time_point> parseDateTime(std::string_view timeField, std::string_view dateField);
};

} // namespace sattrack

#endif
