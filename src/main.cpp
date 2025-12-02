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
#include <map>

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

void printPasses(const std::vector<int> &noradIDs, sattrack::Config &config) {
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
        for (auto noradID : noradIDs) {
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
    config.setTime(std::chrono::system_clock::now());

    auto tleFilename = expandTilde("~/.sattrack.tle");

    auto configFile = expandTilde("~/.sattrack.toml");

    CLI::App app{"SatTrack"};
    argv = app.ensure_utf8(argv);

    app.set_config("--config", configFile, "Read configuration from this file (default: " + configFile + ").");

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
    
    app.ignore_case();

    // n2yo command - access N2YO APIs
    auto n2yoCommand = app.add_subcommand("n2yo", "Get satellite info from N2YO");
    n2yoCommand->add_option_function<std::string>("--key",
        [&config](const std::string &key) { config.setAPIKey(key); },
        "N2YO API key");

    std::vector<int> n2yoIDs;

    auto n2yoTleCommand = n2yoCommand->add_subcommand("tle", "Display raw TLE data from N2YO");
    n2yoTleCommand->add_option_function<std::string>("--key",
        [&config](const std::string &key) { config.setAPIKey(key); },
        "N2YO API key");
    n2yoTleCommand->add_option("id", n2yoIDs, "Norad ID(s) of satellite(s) (ie. 25544)");
    
    auto n2yoInfoCommand = n2yoCommand->add_subcommand("info", "Display a satellite's orbital elements");
    n2yoInfoCommand->add_option_function<std::string>("--key",
        [&config](const std::string &key) { config.setAPIKey(key); },
        "N2YO API key");
    n2yoInfoCommand->add_option("id", n2yoIDs, "Norad ID(s) of satellite(s) (ie. 25544)");

    auto n2yoUpdateCommand = n2yoCommand->add_subcommand("update", "Update local TLE data from N2YO");
    n2yoUpdateCommand->add_option_function<std::string>("--key",
        [&config](const std::string &key) { config.setAPIKey(key); },
        "N2YO API key");
    n2yoUpdateCommand->add_option("id", n2yoIDs, "Norad ID(s) of satellite(s) (ie. 25544)");

    auto n2yoPassesCommand = n2yoCommand->add_subcommand("passes", "Display future passes");
    n2yoPassesCommand->add_option_function<std::string>("--key",
        [&config](const std::string &key) { config.setAPIKey(key); },
        "N2YO API key");
    n2yoPassesCommand->add_option_function<int>("--days",
        [&config](const int days) { config.setDays(days); },
        "Number of days worth of passes to display (default 1, max 10)");
    n2yoPassesCommand->add_option_function<int>("--elev",
        [&config](const int elev) { config.setMinimumElevation(elev); },
        "Minimum pass elevation above the horizon in degrees (default 10, max 90)");
    n2yoPassesCommand->add_option("id", n2yoIDs, "Norad ID(s) of satellite(s)");

    auto satInfoCommand = app.add_subcommand("info", "View satellite information from the local TLE database");
    
    std::vector<int> satInfoIDs;
    satInfoCommand->add_option("id", satInfoIDs, "Norad ID(s) of satellite(s) (ie. 25544)");

    auto updateCommand = app.add_subcommand("update", "Update local TLE data from Celestrak");
    
    std::vector<std::string> groups;
    updateCommand->add_option<std::vector<std::string>>("group", groups, "Celestrak TLE group(s) to download (default: active)");

    auto geoCommand = app.add_subcommand("geo", "Get geodetic location (lat/long/alt) of the satellite at given time");
    
    std::vector<int> geoIDs;
    geoCommand->add_option("id", geoIDs, "Norad ID(s) of satellite(s) (ie. 25544)");
    
    geoCommand->add_option_function<std::string>("--time",
        [&config](const std::string &timeStr) {
            // Parse time string into time_point using date.h
            std::istringstream in(timeStr);
            std::chrono::system_clock::time_point tp;
            in >> date::parse("%Y-%m-%d %H:%M:%S", tp);
            if (in.fail()) {
                throw std::invalid_argument("Invalid time format (expected YYYY-MM-DD HH:MM:SS UTC): " + timeStr);
            }
            config.setTime(tp);
        }, "Time at which to get geodetic location (format: YYYY-MM-DD HH:MM:SS UTC)"
    );

    // Command callbacks

    n2yoCommand->final_callback([n2yoCommand](void) {
        if (n2yoCommand->get_subcommands().empty()) {
           std::cerr << n2yoCommand->help() << std::endl;
           std::exit(1);
        }
    });

    n2yoTleCommand->final_callback([n2yoTleCommand, &config, &n2yoIDs](void) {
        if (!config.hasAPIKey()) {
            std::cerr << "Please provide your N2YO API key." << std::endl;
            std::cerr << n2yoTleCommand->help() << std::endl;
            std::exit(1);
        }
        if (n2yoIDs.empty()) {
            std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
            std::cerr << n2yoTleCommand->help() << std::endl;
            std::exit(1);
        }
        try {
            for (auto noradID : n2yoIDs) {
                auto response = n2yo::getTLE(config.getAPIKey(), noradID, config.getVerbose());
                std::cout << response.info.name << std::endl;
                std::cout << response.tle << std::endl << std::endl;
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    n2yoInfoCommand->final_callback([n2yoInfoCommand, &config, &n2yoIDs](void) {
        if (!config.hasAPIKey()) {
            std::cerr << "Please provide your N2YO API key." << std::endl;
            std::cerr << n2yoInfoCommand->help() << std::endl;
            std::exit(1);
        }
        if (n2yoIDs.empty()) {
            std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
            std::cerr << n2yoInfoCommand->help() << std::endl;
            std::exit(1);
        }
        try {
            for (auto noradID : n2yoIDs) {
                auto response = n2yo::getTLE(config.getAPIKey(), noradID, config.getVerbose());
                sattrack::Orbit orbit;
                orbit.updateFromTLE(response.info.name, response.tle);
                orbit.printInfo(std::cout);
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    n2yoPassesCommand->final_callback([n2yoPassesCommand, &config, &n2yoIDs](void) {
        if (!config.hasAPIKey()) {
            std::cerr << "Please provide your N2YO API key." << std::endl;
            std::cerr << n2yoPassesCommand->help() << std::endl;
            std::exit(1);
        }
        if (n2yoIDs.empty()) {
            std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
            std::cerr << n2yoPassesCommand->help() << std::endl;
            std::exit(1);
        }
        printPasses(n2yoIDs, config);
    });

    satInfoCommand->final_callback([satInfoCommand, &satInfoIDs, &tleFilename](void) {
        try {
            if (satInfoIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << satInfoCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Orbit> satellites;
            loadTLEDatabase(tleFilename, satellites);
            for (auto noradID : satInfoIDs) {
                if (!satellites.contains(noradID)) {
                    std::cerr << "Satellite with Norad ID " << noradID << " not found in the local TLE database." << std::endl;
                    continue;
                }
                auto orbit = satellites[noradID];
                orbit.printInfo(std::cout);
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    updateCommand->final_callback([&config, &groups, &tleFilename](void) {
        try {
            std::map<int, sattrack::Orbit> satellites;
            loadTLEDatabase(tleFilename, satellites);

            if (groups.empty()) {
                groups.push_back("active");
            }

            for (auto group : groups) {
                std::cout << "Downloading TLE data for group: " << group << std::endl;
                auto response = celestrak::getTLE(group, config.getVerbose());
                for (auto tle : response.entries) {
                    sattrack::Orbit orbit;
                    orbit.updateFromTLE(tle);
                    satellites[orbit.getNoradID()] = orbit;
                }
            }

            saveTLEDatabase(tleFilename, satellites);
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    geoCommand->final_callback([geoCommand, &config, &tleFilename, &geoIDs](void) {
        try {
            if (geoIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << geoCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Orbit> satellites;
            loadTLEDatabase(tleFilename, satellites);
            for (auto noradID : geoIDs) {
                if (!satellites.contains(noradID)) {
                    std::cerr << "Satellite with Norad ID " << noradID << " not found in the local TLE database." << std::endl;
                    continue;
                }
                auto orbit = satellites[noradID];
                double julianDate = sattrack::toJulianDate(config.getTime());
                auto geo = orbit.getGeodeticLocationAtTime(julianDate);
                std::cout << "Satellite: " << orbit.getName() << std::endl;
                std::cout << "  Latitude: " << geo.latInRadians * (180.0 / M_PI) << " deg" << std::endl;
                std::cout << "  Longitude: " << geo.lonInRadians * (180.0 / M_PI) << " deg" << std::endl;
                std::cout << "  Altitude: " << geo.altInKilometers << " km" << std::endl;
                std::cout << std::endl;
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    CLI11_PARSE(app, argc, argv);

    if (app.get_subcommands().empty()) {
        std::cerr << app.help() << std::endl;
        std::exit(1);
    }

    return 0;
}
