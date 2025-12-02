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
#include <filesystem>
#include <fstream>

#include <date/date.h>

namespace sattrack {

// Helper function to trim leading spaces from a string_view
std::string_view trimLeft(const std::string_view &str) {
    auto pos = str.find_first_not_of(' ');
    return pos == std::string_view::npos ? "" : str.substr(pos);
}

// Helper function to trim trailing spaces from a string_view
std::string_view trimRight(const std::string_view &str) {
    auto pos = str.find_last_not_of(' ');
    return pos == std::string_view::npos ? "" : str.substr(0, pos + 1);
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

// WGS84 ellipsoid constants (used by multiple functions)
constexpr double WGS84_A = 6378.137;              // Semi-major axis (km) - equatorial radius
constexpr double WGS84_F = 1.0 / 298.257223563;   // Flattening
constexpr double WGS84_E2 = WGS84_F * (2 - WGS84_F);  // Eccentricity squared ≈ 0.00669437999014

// Convert ECEF coordinates to geodetic latitude, longitude, and altitude
Geodetic ecefToGeodetic(const Vec3 &ecef) {
    double x = ecef.x, y = ecef.y, z = ecef.z;
    double lon = std::atan2(y, x);
    double p = std::sqrt(x*x + y*y);

    // Iterative latitude calculation (Bowring's method)
    double lat = std::atan2(z, p * (1 - WGS84_E2));  // initial guess
    for (int i = 0; i < 10; ++i) {
        double sinLat = std::sin(lat);
        double N = WGS84_A / std::sqrt(1 - WGS84_E2 * sinLat * sinLat);
        lat = std::atan2(z + WGS84_E2 * N * sinLat, p);
    }

    double sinLat = std::sin(lat);
    double N = WGS84_A / std::sqrt(1 - WGS84_E2 * sinLat * sinLat);
    double alt = p / std::cos(lat) - N;

    return {lat, lon, alt};
}

// Convert geodetic coordinates to ECEF (Earth-Centered Earth-Fixed)
// See documentation in orbit.hpp for detailed algorithm explanation
Vec3 Geodetic::toECEF() const {
    double sinLat = std::sin(latInRadians);
    double cosLat = std::cos(latInRadians);
    double sinLon = std::sin(lonInRadians);
    double cosLon = std::cos(lonInRadians);

    // Radius of curvature in the prime vertical
    // This is the distance from the surface to the Z-axis along the ellipsoid normal
    double N = WGS84_A / std::sqrt(1.0 - WGS84_E2 * sinLat * sinLat);

    // ECEF coordinates
    // The (1 - e²) factor in Z accounts for the ellipsoid's polar flattening
    return {
        (N + altInKilometers) * cosLat * cosLon,
        (N + altInKilometers) * cosLat * sinLon,
        (N * (1.0 - WGS84_E2) + altInKilometers) * sinLat
    };
}

// Transform ECEF coordinates to ENU (East-North-Up) local tangent plane
// See documentation in orbit.hpp for detailed algorithm explanation
Vec3 ecefToENU(const Vec3& targetECEF, const Geodetic& observer) {
    // Step 1: Get observer's ECEF position and compute difference vector
    Vec3 observerECEF = observer.toECEF();
    Vec3 diff = targetECEF - observerECEF;

    // Step 2: Precompute trig values for rotation matrix
    double sinLat = std::sin(observer.latInRadians);
    double cosLat = std::cos(observer.latInRadians);
    double sinLon = std::sin(observer.lonInRadians);
    double cosLon = std::cos(observer.lonInRadians);

    // Step 3: Apply rotation matrix (ECEF to ENU)
    // This matrix rotates the ECEF frame to align with the local horizon:
    //   - East points along the local latitude circle (toward increasing longitude)
    //   - North points along the local meridian (toward the pole)
    //   - Up points radially outward (normal to the ellipsoid)
    double east  = -sinLon * diff.x + cosLon * diff.y;
    double north = -sinLat * cosLon * diff.x - sinLat * sinLon * diff.y + cosLat * diff.z;
    double up    =  cosLat * cosLon * diff.x + cosLat * sinLon * diff.y + sinLat * diff.z;

    return {east, north, up};
}

// Compute look angles from observer to a target in ECEF coordinates
LookAngles getLookAngles(const Vec3& satECEF, const Geodetic& observer) {
    // Transform to local ENU coordinates
    Vec3 enu = ecefToENU(satECEF, observer);

    // Compute slant range (straight-line distance)
    double range = enu.magnitude();

    // Compute elevation angle
    // elevation = arcsin(Up / range)
    // When satellite is directly overhead, Up = range, so elevation = π/2
    // When satellite is on horizon, Up = 0, so elevation = 0
    double elevation = std::asin(enu.z / range);

    // Compute azimuth angle
    // azimuth = arctan2(East, North)
    // This gives: 0 = North, π/2 = East, π = South, -π/2 = West
    double azimuth = std::atan2(enu.x, enu.y);

    // Normalize azimuth to [0, 2π)
    if (azimuth < 0) {
        azimuth += 2.0 * std::numbers::pi;
    }

    return {azimuth, elevation, range};
}

// Compute look angles from observer to a satellite at a specific time
LookAngles getLookAngles(const Orbit& orbit, const Geodetic& observer, double julianDate) {
    // Get true anomaly at the requested time
    double trueAnomaly = orbit.getTrueAnomalyAtTime(julianDate);

    // Get ECI coordinates
    Vec3 eci = orbit.getECI(trueAnomaly);

    // Convert ECI to ECEF using Greenwich Sidereal Time
    double gst = gmst(julianDate);
    Vec3 ecef = eciToECEF(eci, gst);

    // Compute look angles
    return getLookAngles(ecef, observer);
}

// Check visibility based on look angles
bool isVisible(const LookAngles& angles, double minElevationInRadians) {
    return angles.elevationInRadians >= minElevationInRadians;
}

// Check visibility of a satellite at a specific time
bool isVisible(const Orbit& orbit, const Geodetic& observer,
               double julianDate, double minElevationInRadians) {
    LookAngles angles = getLookAngles(orbit, observer, julianDate);
    return isVisible(angles, minElevationInRadians);
}

// Helper: Get elevation at a specific time
static double getElevationAtTime(const Orbit& orbit, const Geodetic& observer, time_point t) {
    double jd = toJulianDate(t);
    LookAngles angles = getLookAngles(orbit, observer, jd);
    return angles.elevationInRadians;
}

// Helper: Binary search to find precise threshold crossing time
// Returns the time when elevation crosses the threshold
// 'rising' indicates if we're looking for a rise (true) or set (false)
static time_point binarySearchCrossing(
    const Orbit& orbit,
    const Geodetic& observer,
    double threshold,
    time_point low,
    time_point high,
    bool rising) {

    using namespace std::chrono;

    // Refine to within 1 second
    while (duration_cast<seconds>(high - low).count() > 1) {
        time_point mid = low + (high - low) / 2;
        double elev = getElevationAtTime(orbit, observer, mid);

        bool aboveThreshold = elev >= threshold;

        if (rising) {
            // Looking for rise: below threshold -> above threshold
            // If mid is above, the crossing is before mid
            if (aboveThreshold) {
                high = mid;
            } else {
                low = mid;
            }
        } else {
            // Looking for set: above threshold -> below threshold
            // If mid is above, the crossing is after mid
            if (aboveThreshold) {
                low = mid;
            } else {
                high = mid;
            }
        }
    }

    return rising ? high : low;
}

// Helper: Golden section search to find maximum elevation time
// Returns the time when elevation is at maximum between low and high
static time_point goldenSectionSearchMax(
    const Orbit& orbit,
    const Geodetic& observer,
    time_point low,
    time_point high) {

    using namespace std::chrono;

    constexpr double PHI = 1.618033988749895;  // Golden ratio
    constexpr double RESPHI = 2.0 - PHI;       // 1/phi

    auto duration = high - low;
    time_point x1 = low + duration_cast<system_clock::duration>(duration * RESPHI);
    time_point x2 = high - duration_cast<system_clock::duration>(duration * RESPHI);

    double f1 = getElevationAtTime(orbit, observer, x1);
    double f2 = getElevationAtTime(orbit, observer, x2);

    // Refine to within 1 second
    while (duration_cast<seconds>(high - low).count() > 1) {
        if (f1 > f2) {
            high = x2;
            x2 = x1;
            f2 = f1;
            duration = high - low;
            x1 = low + duration_cast<system_clock::duration>(duration * RESPHI);
            f1 = getElevationAtTime(orbit, observer, x1);
        } else {
            low = x1;
            x1 = x2;
            f1 = f2;
            duration = high - low;
            x2 = high - duration_cast<system_clock::duration>(duration * RESPHI);
            f2 = getElevationAtTime(orbit, observer, x2);
        }
    }

    // Return the midpoint
    return low + (high - low) / 2;
}

// Find the next satellite pass over an observer's location
std::optional<PassInfo> findNextPass(
    const Orbit& orbit,
    const Geodetic& observer,
    double minElevationInRadians,
    time_point startTime) {

    using namespace std::chrono;

    // Search parameters
    constexpr auto COARSE_STEP = seconds(60);        // 60-second steps for coarse search
    constexpr auto MAX_SEARCH = hours(48);           // Search up to 48 hours ahead

    auto searchEnd = startTime + MAX_SEARCH;
    auto currentTime = startTime;

    // Track visibility state
    bool wasVisible = isVisible(orbit, observer, toJulianDate(currentTime), minElevationInRadians);

    // If we start during a pass, skip to the end of it first
    if (wasVisible) {
        while (currentTime < searchEnd) {
            currentTime += COARSE_STEP;
            bool nowVisible = isVisible(orbit, observer, toJulianDate(currentTime), minElevationInRadians);
            if (!nowVisible) {
                wasVisible = false;
                break;
            }
        }
        if (wasVisible) {
            // Still visible after 48 hours - probably a GEO satellite or error
            return std::nullopt;
        }
    }

    // Coarse search: find when satellite rises above threshold
    time_point coarseRise, coarseSet;
    bool foundRise = false;
    bool foundSet = false;

    while (currentTime < searchEnd) {
        currentTime += COARSE_STEP;
        bool nowVisible = isVisible(orbit, observer, toJulianDate(currentTime), minElevationInRadians);

        if (!wasVisible && nowVisible) {
            // Found rise
            coarseRise = currentTime - COARSE_STEP;
            foundRise = true;
            wasVisible = true;
        } else if (wasVisible && !nowVisible) {
            // Found set
            coarseSet = currentTime;
            foundSet = true;
            break;
        }

        wasVisible = nowVisible;
    }

    // Handle edge cases
    if (!foundRise) {
        return std::nullopt;  // No pass found in search window
    }

    if (!foundSet) {
        // Pass started but didn't end - use search end as approximate set
        coarseSet = searchEnd;
    }

    // Binary search refinement for precise rise time
    time_point preciseRise = binarySearchCrossing(
        orbit, observer, minElevationInRadians,
        coarseRise, coarseRise + COARSE_STEP, true);

    // Binary search refinement for precise set time
    time_point preciseSet = binarySearchCrossing(
        orbit, observer, minElevationInRadians,
        coarseSet - COARSE_STEP, coarseSet, false);

    // Golden section search for maximum elevation time
    time_point maxTime = goldenSectionSearchMax(orbit, observer, preciseRise, preciseSet);

    // Compute look angles at each key time
    double riseJD = toJulianDate(preciseRise);
    double maxJD = toJulianDate(maxTime);
    double setJD = toJulianDate(preciseSet);

    LookAngles riseAngles = getLookAngles(orbit, observer, riseJD);
    LookAngles maxAngles = getLookAngles(orbit, observer, maxJD);
    LookAngles setAngles = getLookAngles(orbit, observer, setJD);

    return PassInfo{
        .noradID = orbit.getNoradID(),
        .name = orbit.getName(),
        .riseTime = preciseRise,
        .maxElevationTime = maxTime,
        .setTime = preciseSet,
        .riseAngles = riseAngles,
        .maxAngles = maxAngles,
        .setAngles = setAngles
    };
}

// Get geodetic location (lat, lon, alt) of the satellite at a given time
Geodetic Orbit::getGeodeticLocationAtTime(const time_point tp) const {
    double julianDate = toJulianDate(tp);
    return getGeodeticLocationAtTime(julianDate);
}   

// Get geodetic location (lat, lon, alt) of the satellite at a given time
Geodetic Orbit::getGeodeticLocationAtTime(const double julianDate) const {
    // Get true anomaly at the requested time
    double trueAnomaly = getTrueAnomalyAtTime(julianDate);

    // Get ECI coordinates (the satellite's position in inertial frame)
    Vec3 eci = getECI(trueAnomaly);
    
    // Get GMST at the requested time (the earth's rotation angle)
    double gst = gmst(julianDate);

    // Convert ECI to ECEF (the satellite's position in earth-fixed frame)
    Vec3 ecef = eciToECEF(eci, gst);
    
    // Convert ECEF to geodetic coordinates (latitude/longitude/altitude)
    return ecefToGeodetic(ecef);
}

// Update orbital elements from TLE data with a seperate name
void Orbit::updateFromTLE(const std::string_view &name, const std::string_view &tle) {
    updateFromTLE(tle);
    this->name = std::string(name);
}

// Update orbital elements from TLE data
void Orbit::updateFromTLE(const std::string_view &tle) {
    bool firstLineParsed = false;
    bool secondLineParsed = false;
    for (auto line : tle | std::views::split('\n')) {
        std::string_view lineView(line.begin(), line.end());
        lineView = trimLeft(trimRight(lineView));
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
        } else if (!firstLineParsed && !secondLineParsed && !lineView.empty()) {
            // If neither line has been parsed and this line is not empty, assume it's the name
            name = std::string(lineView);
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

std::string Orbit::getName() const {
    return name;
}

void Orbit::printInfo(std::ostream &os) const {
    os << getName() << std::endl;
    os << "  NORAD ID: " << getNoradID() << std::endl;
    os << "  Classification: " << getClassification() << std::endl;
    os << "  Designator: " << getDesignator() << std::endl;
    auto epoch = std::chrono::floor<std::chrono::seconds>(getEpoch());
    os << "  Epoch: " << date::format("%F %T UTC", epoch) << std::endl;
    os << "  First Derivative of Mean Motion: " << getFirstDerivativeMeanMotion() << std::endl;
    os << "  Second Derivative of Mean Motion: " << getSecondDerivativeMeanMotion() << std::endl;
    os << "  Bstar Drag Term: " << getBstarDragTerm() << std::endl;
    os << "  Element Set Number: " << getElementSetNumber() << std::endl;
    os << "  Inclination: " << getInclination() << " deg" << std::endl;
    os << "  Right Ascension of Ascending Node: " << getRightAscensionOfAscendingNode() << " deg" << std::endl;
    os << "  Eccentricity: " << getEccentricity() << std::endl;
    os << "  Argument of Perigee: " << getArgumentOfPerigee() << " deg" << std::endl;
    os << "  Mean Anomaly: " << getMeanAnomaly() << " deg" << std::endl;
    os << "  Mean Motion: " << getMeanMotion() << " revs per day" << std::endl;
    os << "  Revolution Number at Epoch: " << getRevolutionNumberAtEpoch() << std::endl;
    os << std::endl;
}

// Calculate TLE line checksum (mod 10 sum of digits, with '-' counting as 1)
int calculateChecksum(const std::string& line) {
    int sum = 0;
    for (char c : line) {
        if (c >= '0' && c <= '9') {
            sum += (c - '0');
        } else if (c == '-') {
            sum += 1;
        }
    }
    return sum % 10;
}

// Format a value in TLE exponential notation (e.g., " 00000+0" or " 15237-3" or "-12345-6")
// The format is: [sign]NNNNN[sign]E where NNNNN is 5 digits of mantissa and E is exponent
std::string toTLEExponential(double value) {
    if (value == 0.0) {
        return " 00000+0";
    }

    char sign = (value >= 0) ? ' ' : '-';
    value = std::abs(value);

    // Get exponent (power of 10)
    int exponent = static_cast<int>(std::floor(std::log10(value)));

    // Normalize mantissa to be in range [0.1, 1.0)
    double mantissa = value / std::pow(10.0, exponent + 1);

    // Convert mantissa to 5-digit integer
    int mantissaInt = static_cast<int>(std::round(mantissa * 100000));

    // Handle rounding overflow
    if (mantissaInt >= 100000) {
        mantissaInt = 10000;
        exponent++;
    }

    // Format exponent with sign
    char expSign = (exponent + 1 >= 0) ? '+' : '-';
    int expAbs = std::abs(exponent + 1);

    std::ostringstream ss;
    ss << sign << std::setw(5) << std::setfill('0') << mantissaInt << expSign << expAbs;
    return ss.str();
}

// Format first derivative of mean motion for TLE (e.g., " .00008010" or "-.00012345")
std::string formatFirstDerivative(double value) {
    char sign = (value >= 0) ? ' ' : '-';
    value = std::abs(value);

    std::ostringstream ss;
    ss << sign << '.' << std::setw(8) << std::setfill('0')
       << static_cast<int>(std::round(value * 100000000));
    return ss.str();
}

// Get standard 3-line TLE representation of the orbital elements
std::string Orbit::getTLE() const {
    using namespace std::chrono;

    std::ostringstream tleStream;

    // Line 0: Name
    tleStream << name << '\n';

    // Line 1
    // Compute epoch in TLE format: YYDDD.DDDDDDDD
    auto epochDays = floor<days>(epoch);
    year_month_day ymd{epochDays};
    int year = static_cast<int>(ymd.year());
    int twoDigitYear = year % 100;

    // Day of year (1-366)
    auto yearStart = sys_days{ymd.year()/January/1};
    int dayOfYear = (epochDays - yearStart).count() + 1;

    // Fractional part of day
    auto timeOfDay = epoch - epochDays;
    double fracDay = duration_cast<duration<double, std::ratio<86400>>>(timeOfDay).count();

    // Format epoch string (14 chars: YYDDD.DDDDDDDD)
    std::ostringstream epochStr;
    epochStr << std::setw(2) << std::setfill('0') << twoDigitYear
             << std::setw(3) << std::setfill('0') << dayOfYear
             << '.' << std::setw(8) << std::setfill('0')
             << static_cast<long>(std::round(fracDay * 100000000));

    // Build line 1 without checksum (68 chars), then add checksum
    // TLE Line 1 format (columns 1-indexed):
    // Col 01: Line number (1)
    // Col 02: Space
    // Col 03-07: NORAD Catalog Number (5 digits)
    // Col 08: Classification (U/C/S)
    // Col 09: Space
    // Col 10-17: International Designator (8 chars, left-justified, space-padded)
    // Col 18: Space
    // Col 19-32: Epoch (14 chars: YYDDD.DDDDDDDD)
    // Col 33: Space
    // Col 34-43: First Derivative Mean Motion (10 chars: sX.XXXXXXXX)
    // Col 44: Space
    // Col 45-52: Second Derivative Mean Motion (8 chars: sNNNNNsE)
    // Col 53: Space
    // Col 54-61: BSTAR (8 chars: sNNNNNsE)
    // Col 62: Space
    // Col 63: Ephemeris Type (always 0)
    // Col 64: Space
    // Col 65-68: Element Set Number (4 digits, right-justified, space-padded)
    // Col 69: Checksum
    std::ostringstream line1;
    line1 << "1 "
          << std::setw(5) << std::setfill('0') << noradID
          << classification << ' '
          << std::left << std::setw(8) << std::setfill(' ') << designator << ' '
          << epochStr.str() << ' '
          << formatFirstDerivative(firstDerivativeMeanMotion) << ' '
          << toTLEExponential(secondDerivativeMeanMotion) << ' '
          << toTLEExponential(bstarDragTerm) << ' '
          << "0 "
          << std::right << std::setw(4) << std::setfill(' ') << (elementSetNumber % 10000);

    std::string line1Str = line1.str();
    int checksum1 = calculateChecksum(line1Str);
    tleStream << line1Str << checksum1 << '\n';

    // Build line 2 without checksum (68 chars), then add checksum
    // TLE Line 2 format (columns 1-indexed):
    // Col 01: Line number (2)
    // Col 02: Space
    // Col 03-07: NORAD Catalog Number (5 digits)
    // Col 08: Space
    // Col 09-16: Inclination (8 chars: XXX.XXXX)
    // Col 17: Space
    // Col 18-25: RAAN (8 chars: XXX.XXXX)
    // Col 26: Space
    // Col 27-33: Eccentricity (7 digits, implied decimal)
    // Col 34: Space
    // Col 35-42: Argument of Perigee (8 chars: XXX.XXXX)
    // Col 43: Space
    // Col 44-51: Mean Anomaly (8 chars: XXX.XXXX)
    // Col 52: Space
    // Col 53-63: Mean Motion (11 chars: XX.XXXXXXXX)
    // Col 64-68: Revolution Number at Epoch (5 digits)
    // Col 69: Checksum
    int eccInt = static_cast<int>(std::round(eccentricity * 10000000));

    std::ostringstream line2;
    line2 << "2 "
          << std::setw(5) << std::setfill('0') << noradID << ' '
          << std::right << std::setw(8) << std::setfill(' ')
          << std::fixed << std::setprecision(4) << inclination << ' '
          << std::right << std::setw(8) << std::setfill(' ')
          << std::fixed << std::setprecision(4) << rightAscensionOfAscendingNode << ' '
          << std::setw(7) << std::setfill('0') << eccInt << ' '
          << std::right << std::setw(8) << std::setfill(' ')
          << std::fixed << std::setprecision(4) << argumentOfPerigee << ' '
          << std::right << std::setw(8) << std::setfill(' ')
          << std::fixed << std::setprecision(4) << meanAnomaly << ' '
          << std::right << std::setw(11) << std::setfill(' ')
          << std::fixed << std::setprecision(8) << meanMotion
          << std::setw(5) << std::setfill('0') << (revolutionNumberAtEpoch % 100000);

    std::string line2Str = line2.str();
    int checksum2 = calculateChecksum(line2Str);
    tleStream << line2Str << checksum2 << '\n';

    return tleStream.str();
}

// Load TLE database from file using the standard 3-line TLE format
void loadTLEDatabase(const std::string &filepath, std::map<int, Orbit> &database, std::ostream &logStream) {
    logStream << "Loading TLE database from file: " << filepath << std::endl;
    if (!std::filesystem::exists(filepath)) {
        logStream << "TLE database file does not exist: " + filepath << std::endl;
        return;
    }
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open TLE database file: " + filepath);
    }
    loadTLEDatabase(file, database, logStream);
}

// Load TLE database from stream using the standard 3-line TLE format
void loadTLEDatabase(std::istream &s, std::map<int, Orbit> &database, std::ostream &logStream) {
    bool haveFirstLine = false;
    bool haveSecondLine = false;
    std::string line, line1, line2, nameLine;
    int entriesLoaded = 0;
    logStream << "Loading TLE database... " << std::flush;
    while (std::getline(s, line)) {
        if (line.empty()) continue;

        if (line.starts_with("1 ")) {
            line1 = line;
            haveFirstLine = true;
        } else if (line.starts_with("2 ")) {
            line2 = line;
            haveSecondLine = true;
        } else {
            nameLine = line; // Assume it's the name
        }

        if (haveFirstLine && haveSecondLine) {
            Orbit orbit;

            std::ostringstream tleStream;
            tleStream << nameLine << '\n' << line1 << '\n' << line2 << '\n';

            orbit.updateFromTLE(tleStream.str());

            database[orbit.getNoradID()] = orbit;

            haveFirstLine = false;
            haveSecondLine = false;
            line1.clear();
            line2.clear();
            nameLine.clear();
            entriesLoaded++;
        }
    }
    logStream << " done." << std::endl;
    logStream << "Loaded " << entriesLoaded << " TLE entries." << std::endl;
    logStream << std::flush;
}

// Save TLE database to file using the standard 3-line TLE format
void saveTLEDatabase(const std::string &filepath, const std::map<int, Orbit> &database, std::ostream &logStream) {
    std::ofstream file(filepath);
    logStream << "Saving TLE database to file: " << filepath << std::endl;
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for writing: " + filepath);
    }
    saveTLEDatabase(file, database, logStream);
}

// Save TLE database to stream using the standard 3-line TLE format
void saveTLEDatabase(std::ostream &s, const std::map<int, Orbit> &database, std::ostream &logStream) {
    logStream << "Saving TLE database... " << std::flush;
    int count = 0;
    for (const auto &[id, orbit] : database) {
        s << orbit.getTLE();
        count++;
    }
    logStream << " done." << std::endl;
    logStream << "Saved " << count << " TLE entries." << std::endl;
    logStream << std::flush;
}

} // namespace sattrack