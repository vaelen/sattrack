/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack/orbit.hpp>

#include <charconv>
#include <chrono>
#include <cmath>
#include <ranges>
#include <string>
#include <string_view>

namespace sattrack {

// Helper function to convert substring to numeric type
template <typename T>
inline T toNumber(const std::string_view &str) {
    T value;
    auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), value);
    if (ec != std::errc()) {
        throw std::invalid_argument("Couldn't convert value: " + std::string(str));
    }
    return value;
}

// Helper function to convert exponential substring to numeric type
inline double fromExponentialString(const std::string_view &str) {
    // Example input: "11606-4" -> 0.00011606

    auto pos = str.find('-');
    if (pos == std::string_view::npos) {
        throw std::invalid_argument("Invalid exponential format: " + std::string(str));
    }

    std::string_view baseView = str.substr(0, pos);
    std::string_view exponentView = str.substr(pos + 1);

    double base = toNumber<double>(baseView);;
    int exponent = toNumber<int>(exponentView);
    double value = base * std::pow(10.0, -exponent);
    return value;
}

// Helper function to parse epoch from TLE format (YYDDD.DDDDDDDD)
auto parseEpoch(const std::string_view &epochStr) {
    using namespace std::chrono;

    int y = toNumber<int>(epochStr.substr(0, 2));
    double dayOfYear = toNumber<double>(epochStr.substr(2));

    // Convert two-digit year to four-digit year
    if (y < 57) {
        y += 2000;
    } else {
        y += 1900;
    }

    int wholeDays = static_cast<int>(dayOfYear);
    double fracDays = dayOfYear - wholeDays;
    
    auto date = sys_days{year{y}/January/1} + days{wholeDays - 1};
    auto time = duration_cast<microseconds>(duration<double, std::ratio<86400>>{fracDays});
    
    return date + time;
}

// Update orbital elements from TLE data
void Orbit::updateFromTLE(const std::string_view &tle) {
    bool firstLineParsed = false;
    bool secondLineParsed = false;
    for (auto line : tle | std::views::split('\n')) {
        std::string_view lineView(line.begin(), line.end());
        if (lineView.starts_with("1 ")) {
            // NORAD ID is columns 3-7
            noradID = toNumber<int>(lineView.substr(2, 5));
            // Classification is column 8
            classification = lineView[7];
            // Designator is columns 10-17
            designator = std::string(lineView.substr(9, 8));
            // First Derivative of Mean Motion is columns 34-43
            firstDerivativeMeanMotion = toNumber<double>(lineView.substr(33, 10));
            // Second Derivative of Mean Motion is columns 45-52 (exponential format)
            secondDerivativeMeanMotion = fromExponentialString(lineView.substr(44, 8));
            // Bstar Drag Term is columns 54-61 (exponential format)
            bstarDragTerm = fromExponentialString(lineView.substr(53, 8));
            // Time since epoch is columns 19-32
            epoch = parseEpoch(lineView.substr(18, 14));
            // Element Set Number is columns 65-68
            elementSetNumber = toNumber<int>(lineView.substr(64, 4));
            firstLineParsed = true;
        } else if (lineView.starts_with("2 ")) {
            // Inclination is columns 9-16
            inclination = toNumber<double>(lineView.substr(8, 8));
            // RAAN is columns 18-25
            raan = toNumber<double>(lineView.substr(17, 8));
            // Eccentricity is columns 27-33 (decimal implied)
            eccentricity = toNumber<double>("0." + std::string(lineView.substr(26, 7)));
            // Argument of perigee is columns 35-42
            argumentOfPerigee = toNumber<double>(lineView.substr(34, 8));
            // Mean Anomaly is columns 44-51
            meanAnomaly = toNumber<double>(lineView.substr(43, 8));
            // Mean Motion is columns 53-63
            meanMotion = toNumber<double>(lineView.substr(52, 11));
            // Revolution number at epoch is columns 64-68
            revolutionNumberAtEpoch = toNumber<int>(lineView.substr(63, 5));
            secondLineParsed = true;
        }
        if (firstLineParsed && secondLineParsed) {
            break;
        }
    }
}

int Orbit::getNoradID() const {
    return noradID;
}

char Orbit::getClassification() const {
    return classification;
}

std::string Orbit::getDesignator() const {
    return designator;
}

time_point Orbit::getEpoch() const {
    return epoch;
}

double Orbit::getFirstDerivativeMeanMotion() const {
    return firstDerivativeMeanMotion;
}

double Orbit::getSecondDerivativeMeanMotion() const {
    return secondDerivativeMeanMotion;
}

double Orbit::getBstarDragTerm() const {
    return bstarDragTerm;
}

int Orbit::getElementSetNumber() const {
    return elementSetNumber;
}

double Orbit::getInclination() const {
    return inclination;
}

double Orbit::getRAAN() const {
    return raan;
}

double Orbit::getEccentricity() const {
    return eccentricity;
}

double Orbit::getArgumentOfPerigee() const {
    return argumentOfPerigee;
}

double Orbit::getMeanAnomaly() const {
    return meanAnomaly;
}

double Orbit::getMeanMotion() const {
    return meanMotion;
}

int Orbit::getRevolutionNumberAtEpoch() const {
    return revolutionNumberAtEpoch;
}

} // namespace sattrack