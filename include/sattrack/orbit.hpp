/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_ORBIT_HPP
#define __SATTRACK_ORBIT_HPP

#include <chrono>
#include <string>
#include <string_view>

namespace sattrack {

using time_point = std::chrono::system_clock::time_point;

class Orbit {
public:
    Orbit() = default;
    ~Orbit() = default;

    void updateFromTLE(const std::string_view &tle);

    int getNoradID() const;
    char getClassification() const;
    std::string getDesignator() const;
    time_point getEpoch() const;
    double getFirstDerivativeMeanMotion() const;
    double getSecondDerivativeMeanMotion() const;
    double getBstarDragTerm() const;
    int getElementSetNumber() const;

    double getInclination() const;
    double getRAAN() const;
    double getEccentricity() const;
    double getArgumentOfPerigee() const;
    double getMeanAnomaly() const;
    double getMeanMotion() const;
    int getRevolutionNumberAtEpoch() const;

private:
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
    double raan;
    double eccentricity;
    double argumentOfPerigee;
    double meanAnomaly;
    double meanMotion;
    int revolutionNumberAtEpoch;
};

}

#endif
