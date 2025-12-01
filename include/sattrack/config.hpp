/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_CONFIG_HPP
#define __SATTRACK_CONFIG_HPP

#include <string>
#include <set>
#include <optional>
#include <chrono>

namespace sattrack {

using time_point = std::chrono::system_clock::time_point;

class Config {
public:
    // Empty constructor
    Config() = default;
    ~Config() = default;

    bool hasAPIKey();
    void clearAPIKey();
    std::string getAPIKey();
    void setAPIKey(const std::string &apiKey);

    double getLongitude();
    void setLongitude(const double l);

    double getLatitude();
    void setLatitude(const double l);

    double getAltitude();
    void setAltitude(const double a);

    void addSatellite(const int noradID);
    
    template<typename Container>
    std::enable_if<std::is_same_v<typename Container::value_type, int>>
    addAllSatellites(const Container &noradIDs);

    void removeSatellite(const int noradID);

    template<typename Container>
    std::enable_if<std::is_same_v<typename Container::value_type, int>> 
    removeAllSatellites(const Container &noradIDs);

    void clearSatellites();
    std::set<int> getSatellites();

    bool hasSatellites();

    int getDays();
    void setDays(const int days);

    int getMinimumElevation();
    void setMinimumElevation(const int degrees);

    bool getVerbose();
    void setVerbose(bool);

    time_point getTime();
    void setTime(const time_point tp);

private:
    std::optional<std::string> apiKey;
    double longitude;
    double latitude;
    double altitude;
    std::set<int> satellites;
    int days;
    int minimumElevation;
    bool verbose;
    time_point time;
};

}

#endif