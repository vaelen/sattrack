/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 *
 * SGP4 algorithm based on the Vallado reference implementation.
 * See: https://celestrak.org/software/vallado-sw.php
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
#include <stdexcept>

namespace sattrack {

using time_point = std::chrono::system_clock::time_point;

// ============================================================================
// SGP4 Constants
// ============================================================================

// WGS84 / EGM-96 Constants
constexpr double SGP4_MU = 398600.8;            // Earth gravitational parameter (km³/s²)
constexpr double SGP4_RADIUS_EARTH_KM = 6378.135;  // Earth equatorial radius (km)
constexpr double SGP4_J2 = 0.001082616;         // Second gravitational zonal harmonic
constexpr double SGP4_J3 = -0.00000253881;      // Third gravitational zonal harmonic
constexpr double SGP4_J4 = -0.00000165597;      // Fourth gravitational zonal harmonic
constexpr double SGP4_J3OJ2 = SGP4_J3 / SGP4_J2;
constexpr double SGP4_XKE = 0.0743669161331734132;  // sqrt(GM) in Earth radii^1.5/min
constexpr double SGP4_TUMIN = 13.44683969695931;    // Minutes per time unit
constexpr double SGP4_VKMPERSEC = 7.905366149846074; // km/s per velocity unit
constexpr double SGP4_TWO_PI = 2.0 * M_PI;
constexpr double SGP4_X2O3 = 2.0 / 3.0;

// Astronomical constants
constexpr double J2000_JD = 2451545.0;                      // Julian Date of J2000.0 epoch
constexpr double DAYS_PER_JULIAN_CENTURY = 36525.0;         // Days in a Julian century
constexpr double GMST_AT_J2000 = 280.46061837;              // GMST at J2000.0 epoch (degrees)
constexpr double EARTH_SIDEREAL_RATE = 360.98564736629;     // Earth's rotation rate (deg/day)

// IAU polynomial correction coefficients for long-term variations in Earth's rotation
constexpr double GMST_T2_COEFF = 0.000387933;   // Quadratic correction for precession (T² term)
constexpr double GMST_T3_DIVISOR = 38710000.0;  // Cubic correction divisor (T³ term)

// Degree-radian conversion factors
constexpr double DEGREES_TO_RADIANS = M_PI / 180.0;
constexpr double RADIANS_TO_DEGREES = 180.0 / M_PI;

// ============================================================================
// SGP4 Exception Classes
// ============================================================================

/**
 * Base exception class for SGP4 propagation errors.
 */
class SGP4Exception : public std::runtime_error {
public:
    explicit SGP4Exception(const std::string& msg) : std::runtime_error(msg) {}
};

/**
 * Exception thrown when a satellite has decayed (re-entered atmosphere).
 */
class SatelliteDecayedException : public SGP4Exception {
public:
    SatelliteDecayedException() : SGP4Exception("Satellite has decayed") {}
};

/**
 * Exception thrown when orbital elements are invalid.
 */
class InvalidOrbitException : public SGP4Exception {
public:
    explicit InvalidOrbitException(const std::string& msg) : SGP4Exception(msg) {}
};

// ============================================================================
// Basic Data Types
// ============================================================================

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
     */
    Vec3 normalize() const {
        double mag = magnitude();
        return {x / mag, y / mag, z / mag};
    }
};

/**
 * Geodetic coordinates representing a position on or above Earth's surface.
 */
struct Geodetic {
    double latInRadians;      ///< Geodetic latitude (-π/2 to +π/2, positive = North)
    double lonInRadians;      ///< Longitude (-π to +π, positive = East)
    double altInKilometers;   ///< Altitude above the WGS84 ellipsoid surface

    Vec3 toECEF() const;
};

/**
 * Look angles from an observer to a target (typically a satellite).
 */
struct LookAngles {
    double azimuthInRadians;      ///< Compass direction (0 = North, π/2 = East, etc.)
    double elevationInRadians;    ///< Angle above horizon (0 = horizon, π/2 = overhead)
    double rangeInKilometers;     ///< Slant range (straight-line distance) to the target
};

/**
 * Information about a single satellite pass over a ground station.
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

// ============================================================================
// Orbit Class with SGP4 Propagation
// ============================================================================

/**
 * Represents a satellite orbit and provides SGP4/SDP4 propagation.
 *
 * The SGP4 algorithm is used for near-Earth satellites (orbital period < 225 min)
 * and SDP4 is used for deep-space satellites (orbital period >= 225 min).
 *
 * Usage:
 *   Orbit orbit;
 *   orbit.updateFromTLE(tleString);
 *   Vec3 position = orbit.getECI(julianDate);
 *
 * Note: This class is not thread-safe for the same instance. The mutable SGP4
 * state is initialized lazily on first propagation.
 */
class Orbit {
public:
    Orbit() = default;
    ~Orbit() = default;

    /**
     * Update orbital elements from a TLE string with a separate name.
     */
    void updateFromTLE(const std::string_view &name, const std::string_view &tle);

    /**
     * Update orbital elements from a TLE string (name from TLE line 0).
     */
    void updateFromTLE(const std::string_view &tle);

    // Accessors for TLE data
    std::string getName() const;
    int getNoradID() const;
    char getClassification() const;
    std::string getDesignator() const;
    time_point getEpoch() const;
    double getFirstDerivativeMeanMotion() const;
    double getSecondDerivativeMeanMotion() const;
    double getBstarDragTerm() const;
    int getElementSetNumber() const;

    // Accessors for orbital elements
    double getInclination() const;
    double getRightAscensionOfAscendingNode() const;
    double getEccentricity() const;
    double getArgumentOfPerigee() const;
    double getMeanAnomaly() const;
    double getMeanMotion() const;
    int getRevolutionNumberAtEpoch() const;

    /**
     * Get position in Earth-Centered Inertial (ECI/TEME) coordinates.
     *
     * Uses SGP4/SDP4 propagation to compute the satellite position at the
     * specified Julian Date. The returned coordinates are in the True Equator
     * Mean Equinox (TEME) reference frame, in kilometers.
     *
     * @param julianDate The time for which to compute position
     * @return Position vector in TEME coordinates (km)
     * @throws SatelliteDecayedException if the satellite has decayed
     * @throws InvalidOrbitException if orbital elements are invalid
     */
    Vec3 getECI(double julianDate) const;

    /**
     * Get velocity in Earth-Centered Inertial (ECI/TEME) coordinates.
     *
     * @param julianDate The time for which to compute velocity
     * @return Velocity vector in TEME coordinates (km/s)
     */
    Vec3 getVelocity(double julianDate) const;

    /**
     * Get geodetic location (latitude, longitude, altitude) at a given time.
     */
    Geodetic getGeodeticLocationAtTime(const time_point tp) const;
    Geodetic getGeodeticLocationAtTime(const double julianDate) const;

    /**
     * Print orbital element information to a stream.
     */
    void printInfo(std::ostream &os) const;

    /**
     * Get TLE string representation.
     */
    std::string getTLE() const;

private:
    // ========================================================================
    // TLE Data
    // ========================================================================

    // Satellite Identification
    std::string name;

    // First Line - Satellite Identification
    int noradID = 0;
    char classification = 'U';
    std::string designator;
    time_point epoch;
    double firstDerivativeMeanMotion = 0.0;
    double secondDerivativeMeanMotion = 0.0;
    double bstarDragTerm = 0.0;
    int elementSetNumber = 0;

    // Second Line - Orbital Elements (in degrees, except eccentricity)
    double inclination = 0.0;
    double rightAscensionOfAscendingNode = 0.0;
    double eccentricity = 0.0;
    double argumentOfPerigee = 0.0;
    double meanAnomaly = 0.0;
    double meanMotion = 0.0;  // revolutions per day
    int revolutionNumberAtEpoch = 0;

    // ========================================================================
    // SGP4 State Variables (computed during initialization)
    // ========================================================================

    mutable bool sgp4Initialized_ = false;

    // Epoch in Julian Date (split for precision)
    mutable double jdsatepoch_ = 0.0;      // Integer part
    mutable double jdsatepochF_ = 0.0;     // Fractional part

    // Method flag: 'n' = near-earth (SGP4), 'd' = deep-space (SDP4)
    mutable char method_ = 'n';

    // Initialization flags
    mutable bool isimp_ = false;           // Simple drag flag
    mutable int irez_ = 0;                 // Resonance flag (0=none, 1=1-day, 2=0.5-day)

    // Common orbital parameters
    mutable double a_ = 0.0;               // Semi-major axis (Earth radii)
    mutable double alta_ = 0.0;            // Altitude at apogee
    mutable double altp_ = 0.0;            // Altitude at perigee
    mutable double argpo_ = 0.0;           // Argument of perigee (rad)
    mutable double bstar_ = 0.0;           // Drag term
    mutable double ecco_ = 0.0;            // Eccentricity
    mutable double inclo_ = 0.0;           // Inclination (rad)
    mutable double mo_ = 0.0;              // Mean anomaly (rad)
    mutable double no_kozai_ = 0.0;        // Mean motion (Kozai, rad/min)
    mutable double no_unkozai_ = 0.0;      // Mean motion (un-Kozai'd, rad/min)
    mutable double nodeo_ = 0.0;           // Right ascension (rad)
    mutable double gsto_ = 0.0;            // Greenwich sidereal time at epoch
    mutable double cosio2_ = 0.0;          // cosine of inclination squared
    mutable double eccsq_ = 0.0;           // Eccentricity squared

    // Near-earth coefficients
    mutable double aycof_ = 0.0;
    mutable double con41_ = 0.0;
    mutable double cc1_ = 0.0, cc4_ = 0.0, cc5_ = 0.0;
    mutable double d2_ = 0.0, d3_ = 0.0, d4_ = 0.0;
    mutable double delmo_ = 0.0;
    mutable double eta_ = 0.0;
    mutable double argpdot_ = 0.0;
    mutable double omgcof_ = 0.0;
    mutable double sinmao_ = 0.0;
    mutable double t2cof_ = 0.0, t3cof_ = 0.0, t4cof_ = 0.0, t5cof_ = 0.0;
    mutable double x1mth2_ = 0.0;
    mutable double x7thm1_ = 0.0;
    mutable double mdot_ = 0.0;
    mutable double nodedot_ = 0.0;
    mutable double xlcof_ = 0.0;
    mutable double xmcof_ = 0.0;
    mutable double nodecf_ = 0.0;

    // Deep space coefficients
    mutable double e3_ = 0.0, ee2_ = 0.0;
    mutable double peo_ = 0.0, pgho_ = 0.0, pho_ = 0.0, pinco_ = 0.0, plo_ = 0.0;
    mutable double se2_ = 0.0, se3_ = 0.0, sgh2_ = 0.0, sgh3_ = 0.0, sgh4_ = 0.0;
    mutable double sh2_ = 0.0, sh3_ = 0.0, si2_ = 0.0, si3_ = 0.0, sl2_ = 0.0;
    mutable double sl3_ = 0.0, sl4_ = 0.0;
    mutable double xgh2_ = 0.0, xgh3_ = 0.0, xgh4_ = 0.0;
    mutable double xh2_ = 0.0, xh3_ = 0.0;
    mutable double xi2_ = 0.0, xi3_ = 0.0;
    mutable double xl2_ = 0.0, xl3_ = 0.0, xl4_ = 0.0;
    mutable double xlamo_ = 0.0;
    mutable double zmol_ = 0.0, zmos_ = 0.0;
    mutable double atime_ = 0.0;
    mutable double xli_ = 0.0, xni_ = 0.0;

    // Resonance coefficients
    mutable double d2201_ = 0.0, d2211_ = 0.0;
    mutable double d3210_ = 0.0, d3222_ = 0.0;
    mutable double d4410_ = 0.0, d4422_ = 0.0;
    mutable double d5220_ = 0.0, d5232_ = 0.0, d5421_ = 0.0, d5433_ = 0.0;
    mutable double del1_ = 0.0, del2_ = 0.0, del3_ = 0.0;
    mutable double dedt_ = 0.0, didt_ = 0.0, dmdt_ = 0.0;
    mutable double dnodt_ = 0.0, domdt_ = 0.0;

    // ========================================================================
    // SGP4 Private Methods
    // ========================================================================

    /**
     * Initialize SGP4 state from orbital elements.
     * Called lazily on first propagation.
     */
    void initializeSGP4() const;

    /**
     * Core SGP4/SDP4 propagation.
     * @param tsince Minutes since epoch
     * @param r Output position vector (km)
     * @param v Output velocity vector (km/s)
     */
    void propagateSGP4(double tsince, double r[3], double v[3]) const;

    /**
     * Initialize deep space (SDP4) coefficients.
     */
    void initializeDeepSpace(
        double tc, double& snodm, double& cnodm, double& sinim, double& cosim,
        double& sinomm, double& cosomm, double& day, double& em, double& emsq,
        double& gam, double& rtemsq, double& s1, double& s2, double& s3,
        double& s4, double& s5, double& s6, double& s7, double& ss1,
        double& ss2, double& ss3, double& ss4, double& ss5, double& ss6,
        double& ss7, double& sz1, double& sz2, double& sz3, double& sz11,
        double& sz12, double& sz13, double& sz21, double& sz22, double& sz23,
        double& sz31, double& sz32, double& sz33, double& nm, double& z1,
        double& z2, double& z3, double& z11, double& z12, double& z13,
        double& z21, double& z22, double& z23, double& z31, double& z32,
        double& z33
    ) const;

    /**
     * Apply deep space secular effects.
     */
    void deepSpaceSecular(
        double t, double& em, double& argpm, double& inclm,
        double& nodem, double& mm, double& nm
    ) const;

    /**
     * Apply deep space periodic effects.
     */
    void deepSpacePeriodic(
        double t, double& em, double& inclm, double& nodem,
        double& argpm, double& mm
    ) const;
};

// ============================================================================
// TLE Database Functions
// ============================================================================

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
 */
double toJulianDate(time_point tp);

/**
 * Computes Greenwich Mean Sidereal Time (GMST) for a given Julian Date.
 * @return GMST in radians, normalized to [0, 2π)
 */
double gmst(double julianDate);

/**
 * Converts Earth-Centered Inertial (ECI) to Earth-Centered Earth-Fixed (ECEF).
 */
Vec3 eciToECEF(const Vec3 &eci, double gst);

/**
 * Converts ECEF coordinates to geodetic.
 */
Geodetic ecefToGeodetic(const Vec3 &ecef);

// ============================================================================
// Visibility and Look Angle Functions
// ============================================================================

/**
 * Transforms a position from ECEF to ENU (East-North-Up) coordinates.
 */
Vec3 ecefToENU(const Vec3& targetECEF, const Geodetic& observer);

/**
 * Computes look angles from an observer to a target given in ECEF coordinates.
 */
LookAngles getLookAngles(const Vec3& satECEF, const Geodetic& observer);

/**
 * Computes look angles from an observer to a satellite at a specific time.
 */
LookAngles getLookAngles(const Orbit& orbit, const Geodetic& observer, double julianDate);

/**
 * Checks if a satellite is visible based on its look angles.
 */
bool isVisible(const LookAngles& angles, double minElevationInRadians);

/**
 * Checks if a satellite is visible from an observer at a specific time.
 */
bool isVisible(const Orbit& orbit, const Geodetic& observer,
               double julianDate, double minElevationInRadians);

/**
 * Finds the next satellite pass over an observer's location.
 */
std::optional<PassInfo> findNextPass(
    const Orbit& orbit,
    const Geodetic& observer,
    double minElevationInRadians,
    time_point startTime);

}

#endif
