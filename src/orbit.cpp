/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack/orbit.hpp>

#include <charconv>
#include <chrono>
#include <cmath>
#include <numbers>
#include <ranges>
#include <string>
#include <string_view>

namespace sattrack {

// Astronomical constants
constexpr double J2000_JD = 2451545.0;                      // Julian Date of J2000.0 epoch
constexpr double DAYS_PER_JULIAN_CENTURY = 36525.0;         // Days in a Julian century
constexpr double GMST_AT_J2000 = 280.46061837;              // GMST at J2000.0 epoch (degrees)
constexpr double EARTH_SIDEREAL_RATE = 360.98564736629;     // Earth's rotation rate (deg/day)

// IAU polynomial correction coefficients for long-term variations in Earth's rotation
constexpr double GMST_T2_COEFF = 0.000387933;   // Quadratic correction for precession (T² term)
constexpr double GMST_T3_DIVISOR = 38710000.0;  // Cubic correction divisor (T³ term)

// Degree-radian conversion factors
constexpr double DEGREES_TO_RADIANS = std::numbers::pi / 180.0;
// constexpr double RADIANS_TO_DEGREES = 180.0 / std::numbers::pi;

// Helper function to trim leading spaces from a string_view
std::string_view trimLeft(const std::string_view &str) {
    auto pos = str.find_first_not_of(' ');
    return pos == std::string_view::npos ? "" : str.substr(pos);
}

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

    bool isNegative = true;
    auto pos = str.find('-');
    if (pos == std::string_view::npos) {
        isNegative = false;
        pos = str.find('+');
        if (pos == std::string_view::npos) {
            throw std::invalid_argument("Invalid exponential format: " + std::string(str));
        }
    }

    std::string baseStr = "0." + std::string(str.substr(0, pos));
    std::string_view exponentView = str.substr(pos + 1);

    double base = toNumber<double>(baseStr);
    int exponent = toNumber<int>(exponentView);
    double value = base * std::pow(10.0, isNegative ? -exponent : exponent);
    return value;
}

// Helper function to parse epoch from TLE format (YYDDD.DDDDDDDD)
auto parseEpoch(const std::string_view &epochStr) {
    using namespace std::chrono;

    int y = toNumber<int>(trimLeft(epochStr.substr(0, 2)));
    double dayOfYear = toNumber<double>(trimLeft(epochStr.substr(2)));

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

// Convert a time_point to Julian Date
double toJulianDate(time_point tp) {
    using namespace std::chrono;
    
    auto daysSinceEpoch = duration_cast<duration<double, days::period>>(
        tp.time_since_epoch()
    ).count();
    
    // Unix epoch (1970-01-01) in Julian Date is 2440587.5
    return 2440587.5 + daysSinceEpoch;
}

// Greenwich Mean Sidereal Time in radians
double gmst(double julianDate) {
    // Julian centuries since J2000.0
    double T = (julianDate - J2000_JD) / DAYS_PER_JULIAN_CENTURY;

    // GMST in degrees at 0h UT, then add rotation for time of day
    double gmstInDegrees = GMST_AT_J2000
                    + EARTH_SIDEREAL_RATE * (julianDate - J2000_JD)
                    + GMST_T2_COEFF * T * T
                    - T * T * T / GMST_T3_DIVISOR;
    
    // Normalize to [0, 360)
    gmstInDegrees = std::fmod(gmstInDegrees, 360.0);
    if (gmstInDegrees < 0) gmstInDegrees += 360.0;
    
    return gmstInDegrees * DEGREES_TO_RADIANS;
}

// Convert Earth Centered Inertial (ECI) to Earth Centered Earth Fixed (ECEF) using Greenwich Sidereal Time
Vec3 eciToECEF(const Vec3 &eci, double gst) {
    double cosGST = std::cos(gst);
    double sinGST = std::sin(gst);
    
    return {
         eci.x * cosGST + eci.y * sinGST,
        -eci.x * sinGST + eci.y * cosGST,
         eci.z
    };
}

// Convert ECEF coordinates to geodetic latitude, longitude, and altitude
Geodetic ecefToGeodetic(const Vec3 &ecef) {
    constexpr double a = 6378.137;            // WGS84 semi-major axis (km)
    constexpr double f = 1.0 / 298.257223563; // flattening
    constexpr double e2 = f * (2 - f);        // eccentricity squared
    
    double x = ecef.x, y = ecef.y, z = ecef.z;
    double lon = std::atan2(y, x);
    double p = std::sqrt(x*x + y*y);
    
    // Iterative latitude calculation (Bowring's method)
    double lat = std::atan2(z, p * (1 - e2));  // initial guess
    for (int i = 0; i < 10; ++i) {
        double sinLat = std::sin(lat);
        double N = a / std::sqrt(1 - e2 * sinLat * sinLat);
        lat = std::atan2(z + e2 * N * sinLat, p);
    }
    
    double sinLat = std::sin(lat);
    double N = a / std::sqrt(1 - e2 * sinLat * sinLat);
    double alt = p / std::cos(lat) - N;
    
    return {lat, lon, alt};
}

double Orbit::getMeanAnomalyAtTime(const double julianDate) const {
    constexpr double secondsPerDay = 86400.0;
    
    // Time since epoch in seconds
    double dt = (julianDate - toJulianDate(epoch)) * secondsPerDay;
    
    // Convert revs/day to rad/s
    double n = (meanMotion * 2.0 * std::numbers::pi) / secondsPerDay;
    
    // Mean anomaly at requested time
    double M = (meanAnomaly * DEGREES_TO_RADIANS) + n * dt;
    
    // Normalize M to [0, 2π)
    M = std::fmod(M, 2.0 * std::numbers::pi);
    if (M < 0) M += 2.0 * std::numbers::pi;
    
    return M;
}

// Solve Kepler's equation: M = E - e*sin(E)
// Returns eccentric anomaly in radians
double Orbit::getEccentricAnomalyFromMeanAnomaly(double meanAnomalyInRadians, double tolerance) const {
    double e = eccentricity;
    double M = meanAnomalyInRadians;

    // Normalize M to [0, 2π)
    M = std::fmod(M, 2.0 * std::numbers::pi);
    if (M < 0) M += 2.0 * std::numbers::pi;

    // Initial guess (good for low eccentricity)
    double E = (e < 0.8) ? M : std::numbers::pi;
    
    // Newton-Raphson iteration
    for (int i = 0; i < 50; ++i) {
        double f = E - e * std::sin(E) - M;
        double fPrime = 1.0 - e * std::cos(E);
        double delta = f / fPrime;
        E -= delta;
        
        if (std::abs(delta) < tolerance) break;
    }
    
    return E;
}

// Get true anomaly in radians
double Orbit::getTrueAnomalyFromEccentricAnomaly(double eccentricAnomalyInRadians) const {
    double e = eccentricity;
    double E = eccentricAnomalyInRadians;

    // ν = 2 × atan2(√(1+e) × sin(E/2), √(1-e) × cos(E/2))
    return 2.0 * std::atan2(std::sqrt(1.0 + e) * std::sin(E / 2.0), 
                            std::sqrt(1.0 - e) * std::cos(E / 2.0));
}

// Get true anomaly in radians for a given time
double Orbit::getTrueAnomalyAtTime(const double julianDate) const {
    double M = getMeanAnomalyAtTime(julianDate);
    
    // Solve for eccentric anomaly
    double E = getEccentricAnomalyFromMeanAnomaly(M);
    
    // Convert to true anomaly
    return getTrueAnomalyFromEccentricAnomaly(E);
}

// Get the position of the satellite in the Earth Centered Inertial (ECI) frame at a given time
Vec3 Orbit::getECI(double trueAnomalyInRadians) const {
    constexpr double mu = 398600.4418;  // Earth gravitational parameter (km³/s²)

    double a;      // semi-major axis (km)
    double e;      // eccentricity
    double i;      // inclination (rad)
    double raan;   // right ascension of ascending node (rad)
    double omega;  // argument of perigee (rad)
    double nu;     // true anomaly (rad)
    
    // Convert revs/day to rad/s
    double n = (meanMotion * 2.0 * std::numbers::pi) / 86400.0;

    // Get the length of the semi-major axis using Kepler's third law
    a = std::cbrt(mu / (n * n)); 

    // Get remaining values from orbital elements
    e = eccentricity;
    i = inclination * DEGREES_TO_RADIANS;
    raan = rightAscensionOfAscendingNode * DEGREES_TO_RADIANS;
    omega = argumentOfPerigee * DEGREES_TO_RADIANS;
    nu = trueAnomalyInRadians;

    // Position in perifocal frame
    double p = a * (1 - e * e);  // semi-latus rectum
    double rMag = p / (1 + e * std::cos(nu));
    
    double xPf = rMag * std::cos(nu);
    double yPf = rMag * std::sin(nu);

    // Rotation angles
    double cosO = std::cos(raan),  sinO = std::sin(raan);
    double cosI = std::cos(i),     sinI = std::sin(i);
    double cosW = std::cos(omega), sinW = std::sin(omega);

    // Combined rotation matrix (perifocal to ECI)
    // This is R_z(-Ω) × R_x(-i) × R_z(-ω)
    double r11 = cosO * cosW - sinO * sinW * cosI;
    double r12 = -cosO * sinW - sinO * cosW * cosI;
    double r21 = sinO * cosW + cosO * sinW * cosI;
    double r22 = -sinO * sinW + cosO * cosW * cosI;
    double r31 = sinW * sinI;
    double r32 = cosW * sinI;

    return {
        r11 * xPf + r12 * yPf,
        r21 * xPf + r22 * yPf,
        r31 * xPf + r32 * yPf
    };
}


// Update orbital elements from TLE data
void Orbit::updateFromTLE(const std::string_view &tle) {
    bool firstLineParsed = false;
    bool secondLineParsed = false;
    for (auto line : tle | std::views::split('\n')) {
        std::string_view lineView(line.begin(), line.end());
        if (lineView.starts_with("1 ")) {
            // NORAD ID is columns 3-7
            noradID = toNumber<int>(trimLeft(lineView.substr(2, 5)));
            // Classification is column 8
            classification = lineView[7];
            // Designator is columns 10-17
            designator = std::string(lineView.substr(9, 8));
            // First Derivative of Mean Motion is columns 34-43
            firstDerivativeMeanMotion = toNumber<double>(trimLeft(lineView.substr(33, 10)));
            // Second Derivative of Mean Motion is columns 45-52 (exponential format)
            secondDerivativeMeanMotion = fromExponentialString(trimLeft(lineView.substr(44, 8)));
            // Bstar Drag Term is columns 54-61 (exponential format)
            bstarDragTerm = fromExponentialString(trimLeft(lineView.substr(53, 8)));
            // Time since epoch is columns 19-32
            epoch = parseEpoch(lineView.substr(18, 14));
            // Element Set Number is columns 65-68
            elementSetNumber = toNumber<int>(trimLeft(lineView.substr(64, 4)));
            firstLineParsed = true;
        } else if (lineView.starts_with("2 ")) {
            // Inclination is columns 9-16
            inclination = toNumber<double>(trimLeft(lineView.substr(8, 8)));
            // RAAN is columns 18-25
            rightAscensionOfAscendingNode = toNumber<double>(trimLeft(lineView.substr(17, 8)));
            // Eccentricity is columns 27-33 (decimal implied)
            eccentricity = toNumber<double>("0." + std::string(trimLeft(lineView.substr(26, 7))));
            // Argument of perigee is columns 35-42
            argumentOfPerigee = toNumber<double>(trimLeft(lineView.substr(34, 8)));
            // Mean Anomaly is columns 44-51
            meanAnomaly = toNumber<double>(trimLeft(lineView.substr(43, 8)));
            // Mean Motion is columns 53-63
            meanMotion = toNumber<double>(trimLeft(lineView.substr(52, 11)));
            // Revolution number at epoch is columns 64-68
            revolutionNumberAtEpoch = toNumber<int>(trimLeft(lineView.substr(63, 5)));
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

double Orbit::getRightAscensionOfAscendingNode() const {
    return rightAscensionOfAscendingNode;
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