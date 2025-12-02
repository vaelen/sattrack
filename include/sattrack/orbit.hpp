/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_ORBIT_HPP
#define __SATTRACK_ORBIT_HPP

#include <chrono>
#include <cmath>
#include <optional>
#include <string>
#include <string_view>
#include <iostream>
#include <map>

namespace sattrack {

using time_point = std::chrono::system_clock::time_point;

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
constexpr double RADIANS_TO_DEGREES = 180.0 / std::numbers::pi;

/**
 * 3D vector in Cartesian coordinates.  
 */
struct Vec3 {
    double x, y, z;

    Vec3 operator+(const Vec3& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    Vec3 operator-(const Vec3& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    Vec3 operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    double magnitude() const {
        return std::sqrt(x*x + y*y + z*z);
    }

    double dot(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3 cross(const Vec3& other) const {
        return {
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        };
    }

    /**
     * Returns a unit vector (magnitude = 1) in the same direction as this vector.
     *
     * The normalized vector is computed as: v / |v|
     *
     * WARNING: If the vector has zero magnitude, this will result in division
     * by zero, producing NaN or Inf values. Check magnitude() > 0 before calling
     * if the vector might be zero.
     */
    Vec3 normalize() const {
        double mag = magnitude();
        return {x / mag, y / mag, z / mag};
    }
};

/**
 * Geodetic coordinates representing a position on or above Earth's surface.
 *
 * Uses the WGS84 (World Geodetic System 1984) reference ellipsoid, which is
 * the standard for GPS and most modern mapping applications.
 */
struct Geodetic {
    double latInRadians;      ///< Geodetic latitude (-π/2 to +π/2, positive = North)
    double lonInRadians;      ///< Longitude (-π to +π, positive = East)
    double altInKilometers;   ///< Altitude above the WGS84 ellipsoid surface

    /**
     * Converts geodetic coordinates to Earth-Centered Earth-Fixed (ECEF) coordinates.
     *
     * ECEF is a Cartesian coordinate system with:
     * - Origin at Earth's center of mass
     * - X-axis pointing toward the intersection of the equator and prime meridian (0°N, 0°E)
     * - Y-axis pointing toward 0°N, 90°E
     * - Z-axis pointing toward the North Pole
     *
     * The conversion uses the WGS84 ellipsoid parameters:
     * - Semi-major axis (equatorial radius): a = 6378.137 km
     * - Flattening: f = 1/298.257223563
     * - Eccentricity squared: e² = f(2-f) ≈ 0.00669437999014
     *
     * Algorithm:
     * 1. Compute the radius of curvature in the prime vertical (N):
     *    N = a / √(1 - e² × sin²(lat))
     *    This accounts for Earth's ellipsoidal shape - N is larger at the equator
     *    than at the poles.
     *
     * 2. Compute ECEF coordinates:
     *    X = (N + alt) × cos(lat) × cos(lon)
     *    Y = (N + alt) × cos(lat) × sin(lon)
     *    Z = (N × (1 - e²) + alt) × sin(lat)
     *
     * The (1 - e²) factor in Z accounts for the ellipsoid's polar flattening.
     *
     * @return Vec3 containing ECEF coordinates in kilometers
     */
    Vec3 toECEF() const;
};

/**
 * Look angles from an observer to a target (typically a satellite).
 *
 * These angles describe where to point an antenna or telescope to track
 * a satellite from a ground station.
 */
struct LookAngles {
    double azimuthInRadians;      ///< Compass direction (0 = North, π/2 = East, π = South, 3π/2 = West)
    double elevationInRadians;    ///< Angle above horizon (0 = horizon, π/2 = directly overhead)
    double rangeInKilometers;     ///< Slant range (straight-line distance) to the target
};

/**
 * Information about a single satellite pass over a ground station.
 *
 * A "pass" is the period during which a satellite is visible above the
 * observer's horizon (or above a specified minimum elevation threshold).
 */
struct PassInfo {
    int noradID;                  ///< NORAD Catalog ID of the satellite
    std::string name;             ///< Name of the satellite
    time_point riseTime;          ///< When the satellite rises above the elevation threshold
    time_point maxElevationTime;  ///< When the satellite reaches its highest point
    time_point setTime;           ///< When the satellite sets below the elevation threshold
    LookAngles riseAngles;        ///< Azimuth/elevation at rise
    LookAngles maxAngles;         ///< Azimuth/elevation at maximum elevation
    LookAngles setAngles;         ///< Azimuth/elevation at set
};

class Orbit {
public:
    Orbit() = default;
    ~Orbit() = default;

    void updateFromTLE(const std::string_view &name, const std::string_view &tle);
    void updateFromTLE(const std::string_view &tle);

    std::string getName() const;
    int getNoradID() const;
    char getClassification() const;
    std::string getDesignator() const;
    time_point getEpoch() const;
    double getFirstDerivativeMeanMotion() const;
    double getSecondDerivativeMeanMotion() const;
    double getBstarDragTerm() const;
    int getElementSetNumber() const;

    double getInclination() const;
    double getRightAscensionOfAscendingNode() const;
    double getEccentricity() const;
    double getArgumentOfPerigee() const;
    double getMeanAnomaly() const;
    double getMeanMotion() const;
    int getRevolutionNumberAtEpoch() const;
    
    double getMeanAnomalyAtTime(const double julianDate) const;
    double getEccentricAnomalyFromMeanAnomaly(double meanAnomalyInRadians, double tolerance = 1e-12) const;
    double getTrueAnomalyFromEccentricAnomaly(double eccentricAnomalyInRadians) const;
    double getTrueAnomalyAtTime(const double julianDate) const;
    Vec3 getECI(double trueAnomalyInRadians) const;
    Geodetic getGeodeticLocationAtTime(const time_point tp) const;
    Geodetic getGeodeticLocationAtTime(const double julianDate) const;

    void printInfo(std::ostream &os) const;
    std::string getTLE() const;
private:
// Satellite Identification
    std::string name;

// First Line - Satellite Identification
    int noradID;
    char classification;
    std::string designator;
    time_point epoch;
    double firstDerivativeMeanMotion;
    double secondDerivativeMeanMotion;
    double bstarDragTerm;
    int elementSetNumber;

// Second Line - Orbital Elements
    double inclination;
    double rightAscensionOfAscendingNode;
    double eccentricity;
    double argumentOfPerigee;
    double meanAnomaly;
    double meanMotion;
    int revolutionNumberAtEpoch;
};

void loadTLEDatabase(std::istream &s, std::map<int, Orbit> &database, std::ostream &logStream = std::cerr);
void loadTLEDatabase(const std::string &filepath, std::map<int, Orbit> &database, std::ostream &logStream = std::cerr);
void saveTLEDatabase(std::ostream &s, const std::map<int, Orbit> &database, std::ostream &logStream = std::cerr);
void saveTLEDatabase(const std::string &filepath, const std::map<int, Orbit> &database, std::ostream &logStream = std::cerr);

// TLE formatting utilities
int calculateChecksum(const std::string& line);
std::string toTLEExponential(double value);
std::string formatFirstDerivative(double value);

// ============================================================================
// Coordinate System Transformations and Time Functions
// ============================================================================

/**
 * Converts a time_point to Julian Date.
 *
 * Julian Date (JD) is a continuous count of days since the beginning of the
 * Julian Period (January 1, 4713 BC in the proleptic Julian calendar). It is
 * the standard time representation used in astronomy and orbital mechanics
 * because it provides a uniform time scale without the complications of
 * calendars, leap years, or time zones.
 *
 * Key reference points:
 * - J2000.0 epoch (2000-01-01 12:00:00 TT): JD 2451545.0
 * - Unix epoch (1970-01-01 00:00:00 UTC): JD 2440587.5
 *
 * The conversion adds the Unix epoch's Julian Date to the number of days
 * elapsed since the Unix epoch.
 *
 * @param tp The time_point to convert (assumed to be UTC)
 * @return The corresponding Julian Date as a double
 */
double toJulianDate(time_point tp);

/**
 * Computes Greenwich Mean Sidereal Time (GMST) for a given Julian Date.
 *
 * Sidereal time measures Earth's rotation relative to the stars (rather than
 * the Sun). GMST specifically measures the angle between the prime meridian
 * and the vernal equinox (First Point of Aries).
 *
 * GMST is essential for converting between Earth-Centered Inertial (ECI) and
 * Earth-Centered Earth-Fixed (ECEF) coordinate systems, as it describes how
 * much Earth has rotated at any given instant.
 *
 * Algorithm (IAU 1982 model):
 * 1. Compute Julian centuries T since J2000.0:
 *    T = (JD - 2451545.0) / 36525.0
 *
 * 2. Compute GMST in degrees using the polynomial:
 *    GMST = 280.46061837 + 360.98564736629 × (JD - 2451545.0)
 *           + 0.000387933 × T² - T³/38710000
 *
 * The coefficients account for:
 * - 280.46061837°: GMST at J2000.0 epoch
 * - 360.98564736629°/day: Earth's sidereal rotation rate
 * - Higher-order terms: precession and long-term variations
 *
 * @param julianDate The Julian Date for which to compute GMST
 * @return GMST in radians, normalized to [0, 2π)
 */
double gmst(double julianDate);

/**
 * Converts Earth-Centered Inertial (ECI) to Earth-Centered Earth-Fixed (ECEF).
 *
 * ECI and ECEF are both Cartesian coordinate systems centered at Earth's
 * center of mass, but they differ in how they handle Earth's rotation:
 *
 * - ECI (Inertial): Axes are fixed relative to the stars. The X-axis points
 *   toward the vernal equinox, Z-axis toward the celestial north pole.
 *   Used for orbital mechanics because Newton's laws apply directly.
 *
 * - ECEF (Earth-Fixed): Axes rotate with Earth. X-axis points toward 0°N, 0°E
 *   (intersection of equator and prime meridian), Z-axis toward the North Pole.
 *   Used for ground-based positions and GPS.
 *
 * The transformation is a rotation about the Z-axis by the Greenwich Sidereal
 * Time (GST), which represents how much Earth has rotated:
 *
 *     | x_ecef |   |  cos(GST)  sin(GST)  0 |   | x_eci |
 *     | y_ecef | = | -sin(GST)  cos(GST)  0 | × | y_eci |
 *     | z_ecef |   |     0         0      1 |   | z_eci |
 *
 * Note: Z component is unchanged because Earth rotates about its polar axis.
 *
 * @param eci Position vector in ECI coordinates (km)
 * @param gst Greenwich Sidereal Time in radians
 * @return Position vector in ECEF coordinates (km)
 */
Vec3 eciToECEF(const Vec3 &eci, double gst);

/**
 * Converts Earth-Centered Earth-Fixed (ECEF) coordinates to geodetic.
 *
 * Geodetic coordinates (latitude, longitude, altitude) describe a position
 * relative to a reference ellipsoid (WGS84). This is the inverse of
 * Geodetic::toECEF().
 *
 * The conversion is non-trivial because:
 * - Earth is an ellipsoid, not a sphere
 * - Geodetic latitude is the angle between the ellipsoid normal and the
 *   equatorial plane, NOT the angle from Earth's center
 *
 * Algorithm (Bowring's iterative method):
 * 1. Compute longitude directly: λ = atan2(y, x)
 *
 * 2. Compute horizontal distance from Z-axis: p = √(x² + y²)
 *
 * 3. Initial latitude guess: φ₀ = atan2(z, p × (1 - e²))
 *
 * 4. Iterate to refine latitude (typically converges in 2-3 iterations):
 *    - N = a / √(1 - e² × sin²(φ))  (radius of curvature)
 *    - φ = atan2(z + e² × N × sin(φ), p)
 *
 * 5. Compute altitude: h = p / cos(φ) - N
 *
 * Uses WGS84 ellipsoid parameters:
 * - Semi-major axis: a = 6378.137 km
 * - Flattening: f = 1/298.257223563
 *
 * @param ecef Position in ECEF coordinates (km)
 * @return Geodetic coordinates (latitude and longitude in radians, altitude in km)
 */
Geodetic ecefToGeodetic(const Vec3 &ecef);

// ============================================================================
// Visibility and Look Angle Functions
// ============================================================================

/**
 * Transforms a position from ECEF to ENU (East-North-Up) coordinates.
 *
 * ENU is a local tangent plane coordinate system centered at the observer:
 * - East (E): Points toward geographic East (tangent to latitude circle)
 * - North (N): Points toward geographic North (tangent to meridian)
 * - Up (U): Points radially outward from Earth's center (normal to ellipsoid)
 *
 * This transformation is essential for computing look angles because it
 * converts the satellite's position into a reference frame aligned with
 * the observer's local horizon.
 *
 * Algorithm:
 * 1. Compute the difference vector: Δ = target_ECEF - observer_ECEF
 *    This gives the vector from observer to target in ECEF coordinates.
 *
 * 2. Apply the rotation matrix to transform from ECEF to ENU:
 *    The rotation depends on the observer's latitude (φ) and longitude (λ):
 *
 *    | E |   | -sin(λ)        cos(λ)         0      |   | Δx |
 *    | N | = | -sin(φ)cos(λ) -sin(φ)sin(λ)  cos(φ) | × | Δy |
 *    | U |   |  cos(φ)cos(λ)  cos(φ)sin(λ)  sin(φ) |   | Δz |
 *
 * This rotation matrix is derived by:
 * - First rotating about Z-axis by -λ to align X with the local meridian
 * - Then rotating about the new Y-axis by (π/2 - φ) to align Z with local up
 *
 * @param targetECEF The target position in ECEF coordinates (km)
 * @param observer The observer's geodetic location
 * @return Vec3 with (East, North, Up) components in kilometers
 */
Vec3 ecefToENU(const Vec3& targetECEF, const Geodetic& observer);

/**
 * Computes look angles from an observer to a target given in ECEF coordinates.
 *
 * @param satECEF The satellite position in ECEF coordinates (km)
 * @param observer The observer's geodetic location
 * @return LookAngles containing azimuth, elevation, and range
 */
LookAngles getLookAngles(const Vec3& satECEF, const Geodetic& observer);

/**
 * Computes look angles from an observer to a satellite at a specific time.
 *
 * This is a convenience function that combines orbital propagation with
 * look angle calculation. It:
 * 1. Computes the satellite's position at the given Julian date
 * 2. Transforms to ECEF coordinates
 * 3. Computes the look angles from the observer
 *
 * @param orbit The satellite's orbital elements
 * @param observer The observer's geodetic location
 * @param julianDate The time for which to compute look angles
 * @return LookAngles containing azimuth, elevation, and range
 */
LookAngles getLookAngles(const Orbit& orbit, const Geodetic& observer, double julianDate);

/**
 * Checks if a satellite is visible based on its look angles.
 *
 * A satellite is considered "visible" if its elevation angle is at or above
 * the specified minimum elevation threshold.
 *
 * @param angles The pre-computed look angles to the satellite
 * @param minElevationInRadians The minimum elevation threshold (typically 0 to 10°)
 * @return true if the satellite's elevation meets or exceeds the threshold
 */
bool isVisible(const LookAngles& angles, double minElevationInRadians);

/**
 * Checks if a satellite is visible from an observer at a specific time.
 *
 * Combines orbital propagation, coordinate transformation, and visibility
 * checking into a single convenience function.
 *
 * @param orbit The satellite's orbital elements
 * @param observer The observer's geodetic location
 * @param julianDate The time to check visibility
 * @param minElevationInRadians The minimum elevation threshold
 * @return true if the satellite is visible at the specified time
 */
bool isVisible(const Orbit& orbit, const Geodetic& observer,
               double julianDate, double minElevationInRadians);

/**
 * Finds the next satellite pass over an observer's location.
 *
 * A "pass" begins when the satellite rises above the minimum elevation
 * threshold and ends when it sets below that threshold. This function
 * searches forward in time from the start time to find the next complete pass.
 *
 * Algorithm:
 *
 * 1. Coarse search: Steps through time in 60-second increments for up to 48
 *    hours, detecting when elevation transitions above the threshold (rise).
 *    If starting during an existing pass, first searches for the end of that
 *    pass before looking for the next one.
 *
 * 2. Binary search refinement: Once approximate rise and set times are found,
 *    uses binary search to find the precise crossing times within 1 second
 *    accuracy. The search refines the time interval until the gap is < 1 second.
 *
 * 3. Golden section search: Finds the time of maximum elevation between
 *    rise and set times using the golden section algorithm, which efficiently
 *    locates the maximum of a unimodal function without derivatives.
 *
 * Special cases handled:
 * - Starting during an existing pass: Finds the end of current pass, then
 *   continues searching for the next complete pass.
 * - GEO satellites: May not find a pass if the satellite doesn't rise/set.
 * - No pass within window: Returns std::nullopt after searching 48 hours.
 *
 * @param orbit The satellite's orbital elements
 * @param observer The observer's geodetic location
 * @param minElevationInRadians The minimum elevation threshold (typically 0 to 10°)
 * @param startTime The time to begin searching from
 * @return PassInfo for the next pass, or std::nullopt if no pass is found
 *         within a 48-hour search window
 */
std::optional<PassInfo> findNextPass(
    const Orbit& orbit,
    const Geodetic& observer,
    double minElevationInRadians,
    time_point startTime);

}

#endif
