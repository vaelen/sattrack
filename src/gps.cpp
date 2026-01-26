/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack/gps.hpp>
#include <spdlog/spdlog.h>

#include <charconv>
#include <cmath>
#include <date/date.h>

using spdlog::debug;

namespace sattrack {

// === Public methods ===

bool GPS::update(std::string_view sentence) {
    // Remove leading $ if present
    if (!sentence.empty() && sentence[0] == '$') {
        sentence.remove_prefix(1);
    }

    // Validate checksum
    if (!validateChecksum(sentence)) {
        debug("GPS: Invalid checksum");
        return false;
    }

    // Find the asterisk to get just the data portion
    auto asteriskPos = sentence.find('*');
    if (asteriskPos == std::string_view::npos) {
        return false;
    }
    sentence = sentence.substr(0, asteriskPos);

    // Split into fields
    auto fields = splitFields(sentence);
    if (fields.empty()) {
        return false;
    }

    // Determine sentence type (skip talker ID, check positions 2-5)
    std::string_view sentenceId = fields[0];
    if (sentenceId.size() < 5) {
        return false;
    }

    // Extract the sentence type (last 3 characters of the identifier)
    std::string_view sentenceType = sentenceId.substr(sentenceId.size() - 3);

    std::lock_guard<std::mutex> lock(mutex_);

    if (sentenceType == "GGA") {
        return parseGGA(fields);
    } else if (sentenceType == "RMC") {
        return parseRMC(fields);
    } else if (sentenceType == "GSA") {
        return parseGSA(fields);
    }

    // Unknown sentence type - not an error, just ignore
    return false;
}

bool GPS::hasValidPosition() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return position_.has_value();
}

std::optional<Geodetic> GPS::getPosition() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return position_;
}

bool GPS::hasValidTime() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return utcTime_.has_value();
}

std::optional<time_point> GPS::getUTCTime() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return utcTime_;
}

GPSFixQuality GPS::getFixQuality() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return fixQuality_;
}

GPSFixType GPS::getFixType() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return fixType_;
}

bool GPS::hasFix() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return fixType_ == GPSFixType::Fix2D || fixType_ == GPSFixType::Fix3D;
}

bool GPS::has3DFix() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return fixType_ == GPSFixType::Fix3D;
}

int GPS::getSatellitesInUse() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return satellitesInUse_;
}

GPSDataStatus GPS::getDataStatus() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return dataStatus_;
}

std::optional<double> GPS::getHDOP() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return hdop_;
}

std::optional<double> GPS::getVDOP() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return vdop_;
}

std::optional<double> GPS::getPDOP() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return pdop_;
}

std::optional<double> GPS::getSpeedKnots() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return speedOverGroundKnots_;
}

std::optional<double> GPS::getCourseDegrees() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return courseOverGroundDegrees_;
}

void GPS::reset() {
    std::lock_guard<std::mutex> lock(mutex_);
    position_ = std::nullopt;
    utcTime_ = std::nullopt;
    fixQuality_ = GPSFixQuality::Invalid;
    satellitesInUse_ = 0;
    hdop_ = std::nullopt;
    fixType_ = GPSFixType::NoFix;
    pdop_ = std::nullopt;
    vdop_ = std::nullopt;
    speedOverGroundKnots_ = std::nullopt;
    courseOverGroundDegrees_ = std::nullopt;
    dataStatus_ = GPSDataStatus::Void;
}

// === Private methods ===

bool GPS::validateChecksum(std::string_view sentence) const {
    // Find the asterisk
    auto asteriskPos = sentence.find('*');
    if (asteriskPos == std::string_view::npos || asteriskPos + 2 >= sentence.size()) {
        return false;
    }

    // Calculate checksum of data between start and *
    auto data = sentence.substr(0, asteriskPos);
    uint8_t calculated = calculateChecksum(data);

    // Parse the provided checksum (two hex digits after *)
    auto checksumStr = sentence.substr(asteriskPos + 1, 2);
    unsigned int provided = 0;
    auto result = std::from_chars(checksumStr.data(), checksumStr.data() + 2, provided, 16);
    if (result.ec != std::errc{}) {
        return false;
    }

    return calculated == static_cast<uint8_t>(provided);
}

uint8_t GPS::calculateChecksum(std::string_view data) {
    uint8_t checksum = 0;
    for (char c : data) {
        checksum ^= static_cast<uint8_t>(c);
    }
    return checksum;
}

std::vector<std::string_view> GPS::splitFields(std::string_view sentence) {
    std::vector<std::string_view> fields;
    size_t start = 0;

    while (start <= sentence.size()) {
        auto pos = sentence.find(',', start);
        if (pos == std::string_view::npos) {
            fields.push_back(sentence.substr(start));
            break;
        }
        fields.push_back(sentence.substr(start, pos - start));
        start = pos + 1;
    }

    return fields;
}

// GGA format: $GPGGA,hhmmss.ss,llll.ll,a,yyyyy.yy,a,x,xx,x.x,x.x,M,x.x,M,x.x,xxxx*hh
// Fields:
//  0: Sentence ID
//  1: UTC time (hhmmss.ss)
//  2: Latitude (ddmm.mmmm)
//  3: N/S indicator
//  4: Longitude (dddmm.mmmm)
//  5: E/W indicator
//  6: GPS quality (0-8)
//  7: Number of satellites in use
//  8: HDOP
//  9: Altitude above MSL
// 10: Altitude units (M)
// 11: Geoid separation
// 12: Geoid units (M)
// 13: Age of differential GPS data
// 14: Differential reference station ID
bool GPS::parseGGA(const std::vector<std::string_view>& fields) {
    if (fields.size() < 10) {
        return false;
    }

    // Parse fix quality (field 6)
    auto quality = parseInt(fields[6]);
    if (quality) {
        fixQuality_ = static_cast<GPSFixQuality>(*quality);
    } else {
        fixQuality_ = GPSFixQuality::Invalid;
    }

    // Only parse position if we have a fix
    if (fixQuality_ != GPSFixQuality::Invalid) {
        // Parse latitude (fields 2, 3)
        char latHemi = fields[3].empty() ? '\0' : fields[3][0];
        auto lat = parseLatitude(fields[2], latHemi);

        // Parse longitude (fields 4, 5)
        char lonHemi = fields[5].empty() ? '\0' : fields[5][0];
        auto lon = parseLongitude(fields[4], lonHemi);

        // Parse altitude (field 9)
        auto alt = parseDouble(fields[9]);

        if (lat && lon) {
            Geodetic pos;
            pos.latInRadians = *lat;
            pos.lonInRadians = *lon;
            pos.altInKilometers = alt.value_or(0.0) / 1000.0; // Convert m to km
            position_ = pos;
        }
    } else {
        position_ = std::nullopt;
    }

    // Parse satellites in use (field 7)
    if (auto sats = parseInt(fields[7])) {
        satellitesInUse_ = *sats;
    }

    // Parse HDOP (field 8)
    hdop_ = parseDouble(fields[8]);

    return true;
}

// RMC format: $GPRMC,hhmmss.ss,A,llll.ll,a,yyyyy.yy,a,x.x,x.x,ddmmyy,x.x,a*hh
// Fields:
//  0: Sentence ID
//  1: UTC time (hhmmss.ss)
//  2: Status (A=active, V=void)
//  3: Latitude (ddmm.mmmm)
//  4: N/S indicator
//  5: Longitude (dddmm.mmmm)
//  6: E/W indicator
//  7: Speed over ground (knots)
//  8: Course over ground (degrees)
//  9: Date (ddmmyy)
// 10: Magnetic variation
// 11: E/W indicator
bool GPS::parseRMC(const std::vector<std::string_view>& fields) {
    if (fields.size() < 10) {
        return false;
    }

    // Parse status (field 2)
    if (!fields[2].empty()) {
        dataStatus_ = (fields[2][0] == 'A') ?
            GPSDataStatus::Active : GPSDataStatus::Void;
    }

    // Only parse data if status is active
    if (dataStatus_ == GPSDataStatus::Active) {
        // Parse time and date (fields 1, 9)
        utcTime_ = parseDateTime(fields[1], fields[9]);

        // Parse speed over ground (field 7)
        speedOverGroundKnots_ = parseDouble(fields[7]);

        // Parse course over ground (field 8)
        courseOverGroundDegrees_ = parseDouble(fields[8]);
    }

    return true;
}

// GSA format: $GPGSA,A,x,xx,xx,...,x.x,x.x,x.x*hh
// Fields:
//  0: Sentence ID
//  1: Mode (A=auto, M=manual)
//  2: Fix type (1=no fix, 2=2D, 3=3D)
//  3-14: PRN of satellites used (12 fields)
// 15: PDOP
// 16: HDOP
// 17: VDOP
bool GPS::parseGSA(const std::vector<std::string_view>& fields) {
    if (fields.size() < 18) {
        return false;
    }

    // Parse fix type (field 2)
    if (auto ft = parseInt(fields[2])) {
        fixType_ = static_cast<GPSFixType>(*ft);
    }

    // Parse PDOP (field 15)
    pdop_ = parseDouble(fields[15]);

    // Parse HDOP (field 16)
    auto hdop = parseDouble(fields[16]);
    if (hdop) {
        hdop_ = hdop;
    }

    // Parse VDOP (field 17)
    vdop_ = parseDouble(fields[17]);

    return true;
}

std::optional<double> GPS::parseLatitude(std::string_view field, char hemisphere) {
    if (field.empty() || (hemisphere != 'N' && hemisphere != 'S')) {
        return std::nullopt;
    }

    // Format: ddmm.mmmm (2 digits for degrees)
    if (field.size() < 4) {
        return std::nullopt;
    }

    // Parse degrees (first 2 digits)
    int degrees = 0;
    auto result = std::from_chars(field.data(), field.data() + 2, degrees);
    if (result.ec != std::errc{}) {
        return std::nullopt;
    }

    // Parse minutes (rest of the field)
    double minutes = 0.0;
    // std::from_chars for double may not be available on all platforms
    // Use a manual approach
    std::string minutesStr(field.substr(2));
    try {
        minutes = std::stod(minutesStr);
    } catch (...) {
        return std::nullopt;
    }

    // Convert to decimal degrees, then to radians
    double decimalDegrees = degrees + (minutes / 60.0);
    if (hemisphere == 'S') {
        decimalDegrees = -decimalDegrees;
    }

    return decimalDegrees * DEGREES_TO_RADIANS;
}

std::optional<double> GPS::parseLongitude(std::string_view field, char hemisphere) {
    if (field.empty() || (hemisphere != 'E' && hemisphere != 'W')) {
        return std::nullopt;
    }

    // Format: dddmm.mmmm (3 digits for degrees)
    if (field.size() < 5) {
        return std::nullopt;
    }

    // Parse degrees (first 3 digits)
    int degrees = 0;
    auto result = std::from_chars(field.data(), field.data() + 3, degrees);
    if (result.ec != std::errc{}) {
        return std::nullopt;
    }

    // Parse minutes (rest of the field)
    std::string minutesStr(field.substr(3));
    double minutes = 0.0;
    try {
        minutes = std::stod(minutesStr);
    } catch (...) {
        return std::nullopt;
    }

    // Convert to decimal degrees, then to radians
    double decimalDegrees = degrees + (minutes / 60.0);
    if (hemisphere == 'W') {
        decimalDegrees = -decimalDegrees;
    }

    return decimalDegrees * DEGREES_TO_RADIANS;
}

std::optional<double> GPS::parseDouble(std::string_view field) {
    if (field.empty()) {
        return std::nullopt;
    }

    std::string str(field);
    try {
        return std::stod(str);
    } catch (...) {
        return std::nullopt;
    }
}

std::optional<int> GPS::parseInt(std::string_view field) {
    if (field.empty()) {
        return std::nullopt;
    }

    int value = 0;
    auto result = std::from_chars(field.data(), field.data() + field.size(), value);
    if (result.ec != std::errc{}) {
        return std::nullopt;
    }

    return value;
}

std::optional<time_point> GPS::parseDateTime(std::string_view timeField, std::string_view dateField) {
    using namespace std::chrono;

    if (timeField.size() < 6 || dateField.size() < 6) {
        return std::nullopt;
    }

    // Parse time: hhmmss.ss
    int h = 0, m = 0, s = 0;
    auto result = std::from_chars(timeField.data(), timeField.data() + 2, h);
    if (result.ec != std::errc{}) return std::nullopt;

    result = std::from_chars(timeField.data() + 2, timeField.data() + 4, m);
    if (result.ec != std::errc{}) return std::nullopt;

    result = std::from_chars(timeField.data() + 4, timeField.data() + 6, s);
    if (result.ec != std::errc{}) return std::nullopt;

    // Parse fractional seconds if present
    int ms = 0;
    if (timeField.size() > 7 && timeField[6] == '.') {
        std::string fracStr(timeField.substr(7));
        // Pad or truncate to 3 digits for milliseconds
        while (fracStr.size() < 3) fracStr += '0';
        if (fracStr.size() > 3) fracStr = fracStr.substr(0, 3);
        try {
            ms = std::stoi(fracStr);
        } catch (...) {
            ms = 0;
        }
    }

    // Parse date: ddmmyy
    int dd = 0, mm = 0, yy = 0;
    result = std::from_chars(dateField.data(), dateField.data() + 2, dd);
    if (result.ec != std::errc{}) return std::nullopt;

    result = std::from_chars(dateField.data() + 2, dateField.data() + 4, mm);
    if (result.ec != std::errc{}) return std::nullopt;

    result = std::from_chars(dateField.data() + 4, dateField.data() + 6, yy);
    if (result.ec != std::errc{}) return std::nullopt;

    // Convert 2-digit year to 4-digit year
    // Following TLE convention: < 57 = 2000s, >= 57 = 1900s
    int y = (yy < 57) ? (2000 + yy) : (1900 + yy);

    // Use date library from the project - construct year_month_day
    auto ymd = year{y} / month{static_cast<unsigned>(mm)} / day{static_cast<unsigned>(dd)};
    if (!ymd.ok()) {
        return std::nullopt;
    }

    // Create time_point from date and time components
    auto tp = sys_days{ymd} + hours{h} + minutes{m} + seconds{s} + milliseconds{ms};

    return tp;
}

} // namespace sattrack
