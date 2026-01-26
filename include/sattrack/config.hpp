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

    int getDays();
    void setDays(const int days);

    int getMinimumElevation();
    void setMinimumElevation(const int degrees);

    bool getVerbose();
    void setVerbose(bool);

    int getHorizon();
    void setHorizon(int degrees);

    time_point getTime();
    void setTime(const time_point tp);

    std::string getGPSSerialPort();
    void setGPSSerialPort(const std::string &port);

    int getGPSBaudRate();
    void setGPSBaudRate(int baudRate);

    std::string getRotatorSerialPort();
    void setRotatorSerialPort(const std::string &port);

    int getRotatorBaudRate();
    void setRotatorBaudRate(int baudRate);

    std::string getRadioSerialPort();
    void setRadioSerialPort(const std::string &port);

    int getRadioBaudRate();
    void setRadioBaudRate(int baudRate);

private:
    std::optional<std::string> apiKey;
    double longitude;
    double latitude;
    double altitude;
    int days;
    int minimumElevation;
    bool verbose;
    time_point time;
    int horizon;
    std::string gpsSerialPort;
    int gpsBaudRate;
    std::string rotatorSerialPort;
    int rotatorBaudRate;
    std::string radioSerialPort;
    int radioBaudRate;
};

}

#endif