/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack.hpp>
#include <CLI/CLI.hpp>
#include <date/date.h>
#include <filesystem>
#include <fstream>
#include <format>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>

namespace fs = std::filesystem;

/** Replace ~ with HOME directory */
std::string expandTilde(const std::string &path) {
    if (!path.empty() && path[0] == '~') {
        const char *home = std::getenv("HOME");
        if (home) {
            return std::string(home) + path.substr(1);
        }
    }
    return path;
}

std::string formatTimestampUTC(const long timestamp) {
    auto timePoint = std::chrono::system_clock::time_point(
        std::chrono::seconds(timestamp)
    );
    auto truncated = std::chrono::floor<std::chrono::seconds>(timePoint);
    return date::format("%F %T UTC", truncated);
}

std::chrono::system_clock::time_point toTimePoint(const long timestamp) {
    return std::chrono::system_clock::time_point(
        std::chrono::seconds(timestamp)
    );
}

void printPasses(sattrack::Config &config) {
    try {
        constexpr std::string_view rowFormat = "{:^25} {:^25} {:^25} {:^14} {:^18} {:^12} {:^12} {:^12} {:^10}";
        constexpr std::string_view satFormat = " {:<5} {:<17}";
        constexpr const char* timeFormat = "%F %T UTC";
        constexpr std::string_view durationFormat = "{:>3}m {:>2}s";
        constexpr std::string_view etaFormat = "{:>2}h {:>2}m {:>2}s";
        constexpr std::string_view azFormat = "{:>6.2f} {:<3}";
        constexpr std::string_view elFormat = "{:<5.2f}";

        std::string sep25(25, '-');
        std::string sep18(18, '-');
        std::string sep15(15, '-');
        std::string sep14(14, '-');
        std::string sep12(12, '-');
        std::string sep10(10, '-');

        struct RadioPassWrapper {
            n2yo::SatelliteInfo info;
            n2yo::RadioPass pass;
        };

        std::vector<RadioPassWrapper> passes;

        std::cout << "Loading Passes..." << std::flush;
        for (auto noradID : config.getSatellites()) {
            auto response = n2yo::getRadioPasses(config.getAPIKey(),
                noradID, config.getLatitude(), config.getLongitude(), config.getAltitude(),
                config.getDays(), config.getMinimumElevation(), config.getVerbose());
            for (auto pass: response.passes) {
                passes.emplace_back(response.info, pass);
            }
            std::cout << "." << std::flush;
        }
        std::cout << " Done." << std::endl << std::endl << std::flush;

        // Sort by startUTC (ascending)
        auto comparator = [](const RadioPassWrapper& a, const RadioPassWrapper& b) {
            return a.pass.startUTC < b.pass.startUTC;
        };
        std::sort(passes.begin(), passes.end(), comparator);

        std::cout << "Upcoming Passes:" << std::endl;
        std::cout << std::format(rowFormat, "Satellite", "Start", "End", "Duration", "Starts In", "Start Az", "End Az", "Max Az", "Max Elev") << std::endl;
        std::cout << std::format(rowFormat, sep25, sep25, sep25, sep14, sep18, sep12, sep12, sep12, sep10) << std::endl;
        for (auto wrapper: passes) {
            auto info = wrapper.info;
            auto pass = wrapper.pass;
            // Calculate all our times and durations
            auto startTime = toTimePoint(pass.startUTC);
            auto endTime = toTimePoint(pass.endUTC);
            auto duration = std::chrono::seconds(pass.endUTC - pass.startUTC);
            auto durationMins = std::chrono::duration_cast<std::chrono::minutes>(duration).count();
            auto durationSecs = std::chrono::duration_cast<std::chrono::seconds>(duration % std::chrono::minutes(1)).count();
            auto now = std::chrono::system_clock::now();
            auto nowSeconds = std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count();
            auto eta = std::chrono::seconds(pass.startUTC - nowSeconds);
            auto etaHours = std::chrono::duration_cast<std::chrono::hours>(eta).count();
            auto etaMins = std::chrono::duration_cast<std::chrono::minutes>(eta % std::chrono::hours(1)).count();
            auto etaSecs = std::chrono::duration_cast<std::chrono::seconds>(eta % std::chrono::minutes(1)).count();

            // Output the pass data
            auto startTimeSeconds = std::chrono::floor<std::chrono::seconds>(startTime);
            auto endTimeSeconds = std::chrono::floor<std::chrono::seconds>(endTime);
            std::cout << std::format(rowFormat,
                std::format(satFormat, info.id, info.name),
                date::format(timeFormat, startTimeSeconds),
                date::format(timeFormat, endTimeSeconds),
                std::format(durationFormat, durationMins, durationSecs),
                std::format(etaFormat, etaHours, etaMins, etaSecs),
                std::format(azFormat, pass.startAz, pass.startAzCompass),
                std::format(azFormat, pass.endAz, pass.endAzCompass),
                std::format(azFormat, pass.maxAz, pass.maxAzCompass),
                std::format(elFormat, pass.maxEl)) << std::endl;
        }
        std::cout << std::endl;
    } catch (const std::exception &err) {
        std::cerr << err.what() << std::endl;
        std::exit(1);
    }
}


/** Program entry point */
int main(int argc, char* argv[]) {

    sattrack::Config config;
    config.setLatitude(0.0);
    config.setLongitude(0.0);
    config.setAltitude(0.0);
    config.setDays(1);
    config.setMinimumElevation(10);
    config.setVerbose(false);

    auto configFile = expandTilde("~/.sattrack.toml");

    CLI::App app{"SatTrack"};
    argv = app.ensure_utf8(argv);

    app.set_config("--config", configFile, "Read configuration from this file (default: " + configFile + ").");

    std::vector<int> ids;

    app.add_option_function<std::string>("--key",
        [&config](const std::string &key) { config.setAPIKey(key); },
        "N2YO API key");
    app.add_option("--id", ids, "Norad ID of satellite to track (can be listed more than once)");
    app.add_option_function<double>("--lat", 
        [&config](const double l) { config.setLatitude(l); },
        "The latitude of the ground station (in decimal format)");
    app.add_option_function<double>("--long", 
        [&config](const double l) { config.setLongitude(l); },
        "The longitude of the ground station (in decimal format)");
    app.add_option_function<double>("--alt", 
        [&config](const double e) { config.setAltitude(e); },
        "Altitude above sea level in meters");
    app.add_flag_function("-v,--verbose", 
        [&config](const int64_t v) { config.setVerbose(v > 0); },
        "Display debugging information");
    
    app.require_subcommand();
    app.ignore_case();

    auto keyCommand = app.add_subcommand("key", "Set N2YO API key");
    keyCommand->add_option_function<std::string>("key", 
        [&config](const std::string &key) { config.setAPIKey(key); },
        "N2YO API key");

    auto addCommand = app.add_subcommand("add", "Add satellite");
    addCommand->add_option_function<int>("id", 
        [&config](const int id) { config.addSatellite(id); },
        "Norad ID of satellite to start tracking");

    auto removeCommand = app.add_subcommand("remove", "Remove satellite");
    removeCommand->add_option_function<int>("id", 
        [&config](const int id) { config.removeSatellite(id); },
        "Norad ID of satellite to stop tracking");

    auto tleCommand = app.add_subcommand("tle", "Display TLE");
    auto parseTLECommand = tleCommand->add_subcommand("parse", "Parse TLE and display orbital elements");

    auto passesCommand = app.add_subcommand("passes", "Display future passes");
    passesCommand->add_option_function<int>("--days", 
        [&config](const int days) { config.setDays(days); },
        "Number of days worth of passes to display (default 1, max 10)");
    passesCommand->add_option_function<int>("--elev", 
        [&config](const int elev) { config.setMinimumElevation(elev); },
        "Minimum pass elevation above the horizon in degrees (default 10, max 90)");

    auto geoCommand = app.add_subcommand("geo", "Get geodetic location (lat/long/alt) of the satellite at given time");
    geoCommand->add_option_function<int>("id", 
        [&config](const int id) { config.addSatellite(id); },
        "Norad ID of satellite to track");
    geoCommand->add_option_function<std::string>("time",
        [&config](const std::string &timeStr) {
            // Parse time string into time_point using date.h
            std::istringstream in(timeStr);
            std::chrono::system_clock::time_point tp;
            in >> date::parse("%Y-%m-%d %H:%M:%S", tp);
            if (in.fail()) {
                throw std::invalid_argument("Invalid time format (expected YYYY-MM-DD HH:MM:SS UTC): " + timeStr);
            }
            config.setTime(tp);
        }
    );

    app.parse_complete_callback(
        [&config, &ids](void) {
            if (!config.hasAPIKey()) {
                std::cerr << "Please provide your N2YO API key." << std::endl;
                std::exit(1);
            }

            for (auto id: ids) {
                config.addSatellite(id);
            }

            if (!config.hasSatellites()) {
                std::cerr << "Please provide at least one satellite to track." << std::endl;
                std::exit(1);
            }
        }
    );

    tleCommand->final_callback([&config, tleCommand](void) {
        // Skip if a subcommand was invoked
        if (!tleCommand->get_subcommands().empty()) {
            return;
        }
        try {
            for (auto noradID : config.getSatellites()) {
                auto response = n2yo::getTLE(config.getAPIKey(), noradID, config.getVerbose());
                std::cout << response.info.name << std::endl;
                std::cout << response.tle << std::endl << std::endl;
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    parseTLECommand->final_callback([&config](void) {
        try {
            for (auto noradID : config.getSatellites()) {
                auto response = n2yo::getTLE(config.getAPIKey(), noradID, config.getVerbose());
                std::cout << response.info.name << std::endl;
                sattrack::Orbit orbit;
                orbit.updateFromTLE(response.tle);
                std::cout << "  NORAD ID: " << orbit.getNoradID() << std::endl;
                std::cout << "  Classification: " << orbit.getClassification() << std::endl;
                std::cout << "  Designator: " << orbit.getDesignator() << std::endl;
                auto epoch = std::chrono::floor<std::chrono::seconds>(orbit.getEpoch());
                std::cout << "  Epoch: " << date::format("%F %T UTC", epoch) << std::endl;
                std::cout << "  First Derivative of Mean Motion: " << orbit.getFirstDerivativeMeanMotion() << std::endl;
                std::cout << "  Second Derivative of Mean Motion: " << orbit.getSecondDerivativeMeanMotion() << std::endl;
                std::cout << "  Bstar Drag Term: " << orbit.getBstarDragTerm() << std::endl;
                std::cout << "  Element Set Number: " << orbit.getElementSetNumber() << std::endl;
                std::cout << "  Inclination: " << orbit.getInclination() << " deg" << std::endl;
                std::cout << "  Right Ascension of Ascending Node: " << orbit.getRightAscensionOfAscendingNode() << " deg" << std::endl;
                std::cout << "  Eccentricity: " << orbit.getEccentricity() << std::endl;
                std::cout << "  Argument of Perigee: " << orbit.getArgumentOfPerigee() << " deg" << std::endl;
                std::cout << "  Mean Anomaly: " << orbit.getMeanAnomaly() << " deg" << std::endl;
                std::cout << "  Mean Motion: " << orbit.getMeanMotion() << " revs per day" << std::endl;
                std::cout << "  Revolution Number at Epoch: " << orbit.getRevolutionNumberAtEpoch() << std::endl;
                std::cout << std::endl;
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    passesCommand->final_callback([&config](void) {
        printPasses(config);
    });

    CLI11_PARSE(app, argc, argv);

    return 0;
}
