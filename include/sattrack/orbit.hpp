/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_ORBIT_HPP
#define __SATTRACK_ORBIT_HPP

#include <chrono>
#include <cmath>
#include <string>
#include <string_view>
#include <iostream>
#include <map>

namespace sattrack {

using time_point = std::chrono::system_clock::time_point;

struct Vec3 {
    double x, y, z;
    
    Vec3 operator+(const Vec3& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }
    
    Vec3 operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }
    
    double magnitude() const {
        return std::sqrt(x*x + y*y + z*z);
    }
};

// ECEF to geodetic (lat/lon/altitude)
// Uses WGS84 ellipsoid, iterative method
struct Geodetic {
    double latInRadians;
    double lonInRadians;
    double altInKilometers;
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

Geodetic ecefToGeodetic(const Vec3 &ecef);
double toJulianDate(time_point tp);
double gmst(double julianDate);
Vec3 eciToECEF(const Vec3 &eci, double gst);
Geodetic ecefToGeodetic(const Vec3 &ecef);

void loadTLEDatabase(std::istream &s, std::map<int, Orbit> &database, std::ostream &logStream = std::cerr);
void loadTLEDatabase(const std::string &filepath, std::map<int, Orbit> &database, std::ostream &logStream = std::cerr);
void saveTLEDatabase(std::ostream &s, const std::map<int, Orbit> &database, std::ostream &logStream = std::cerr);
void saveTLEDatabase(const std::string &filepath, const std::map<int, Orbit> &database, std::ostream &logStream = std::cerr);

// TLE formatting utilities
int calculateChecksum(const std::string& line);
std::string toTLEExponential(double value);
std::string formatFirstDerivative(double value);

}

#endif
