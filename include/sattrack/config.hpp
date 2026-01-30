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

constexpr const char* DEFAULT_GPS_SERIAL_PORT      = "/dev/ttyS1";
constexpr int         DEFAULT_GPS_BAUD_RATE        = 9600;
constexpr const char* DEFAULT_ROTATOR_SERIAL_PORT  = "/dev/ttyS2";
constexpr int         DEFAULT_ROTATOR_BAUD_RATE    = 9600;
constexpr const char* DEFAULT_RADIO_SERIAL_PORT    = "/dev/ttyS4";
constexpr int         DEFAULT_RADIO_BAUD_RATE      = 38400;
constexpr int         DEFAULT_STATUS_INTERVAL_SECONDS = 300; // 5 minutes
constexpr int         DEFAULT_LCD_I2C_BUS            = 1;
constexpr int         DEFAULT_LCD_I2C_ADDRESS        = 0x27;

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

    int getStatusIntervalSeconds();
    void setStatusIntervalSeconds(int seconds);

    int getLCDI2CBus();
    void setLCDI2CBus(int bus);

    int getLCDI2CAddress();
    void setLCDI2CAddress(int address);
    
private:
    std::optional<std::string> apiKey;
    double longitude = 0.0;
    double latitude = 0.0;
    double altitude = 0.0;
    int days = 0;
    int minimumElevation = 0;
    bool verbose = false;
    time_point time;
    int horizon = 0;
    std::string gpsSerialPort = DEFAULT_GPS_SERIAL_PORT;
    int gpsBaudRate = DEFAULT_GPS_BAUD_RATE;
    std::string rotatorSerialPort = DEFAULT_ROTATOR_SERIAL_PORT;
    int rotatorBaudRate = DEFAULT_ROTATOR_BAUD_RATE;
    std::string radioSerialPort = DEFAULT_RADIO_SERIAL_PORT;
    int radioBaudRate = DEFAULT_RADIO_BAUD_RATE;
    int statusIntervalSeconds = DEFAULT_STATUS_INTERVAL_SECONDS;
    int lcdI2CBus = DEFAULT_LCD_I2C_BUS;
    int lcdI2CAddress = DEFAULT_LCD_I2C_ADDRESS;
};

}

#endif