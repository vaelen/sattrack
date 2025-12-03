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
#include <thread>

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

/** Convert azimuth in radians to compass direction string */
std::string azimuthToCompass(double azimuthRad) {
    double deg = azimuthRad * sattrack::RADIANS_TO_DEGREES;
    if (deg < 0) deg += 360.0;
    if (deg >= 360.0) deg -= 360.0;

    if (deg < 11.25) return "N";
    if (deg < 33.75) return "NNE";
    if (deg < 56.25) return "NE";
    if (deg < 78.75) return "ENE";
    if (deg < 101.25) return "E";
    if (deg < 123.75) return "ESE";
    if (deg < 146.25) return "SE";
    if (deg < 168.75) return "SSE";
    if (deg < 191.25) return "S";
    if (deg < 213.75) return "SSW";
    if (deg < 236.25) return "SW";
    if (deg < 258.75) return "WSW";
    if (deg < 281.25) return "W";
    if (deg < 303.75) return "WNW";
    if (deg < 326.25) return "NW";
    if (deg < 348.75) return "NNW";
    return "N";
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

sattrack::PassInfo convertN2YOPassToPassInfo(
    const n2yo::RadioPass &n2yoPass,
    const n2yo::SatelliteInfo &satInfo) {

    using namespace sattrack;
    using namespace std::chrono;

    PassInfo passInfo;

    passInfo.noradID = satInfo.id;
    passInfo.name = satInfo.name;
    passInfo.riseTime = toTimePoint(n2yoPass.startUTC);
    passInfo.setTime = toTimePoint(n2yoPass.endUTC);
    passInfo.maxElevationTime = toTimePoint(n2yoPass.maxUTC);
    passInfo.riseAngles.azimuthInRadians = n2yoPass.startAz * DEGREES_TO_RADIANS;
    passInfo.riseAngles.elevationInRadians = 0.0;  // N2YO does not provide elevation at rise
    passInfo.riseAngles.rangeInKilometers = 0.0;  // N2YO does not provide range at rise
    passInfo.maxAngles.azimuthInRadians = n2yoPass.maxAz * DEGREES_TO_RADIANS;
    passInfo.maxAngles.elevationInRadians = n2yoPass.maxEl * DEGREES_TO_RADIANS;
    passInfo.maxAngles.rangeInKilometers = 0.0;  // N2YO does not provide range at max elevation
    passInfo.setAngles.azimuthInRadians = n2yoPass.endAz * DEGREES_TO_RADIANS;
    passInfo.setAngles.elevationInRadians = 0.0;  // N2YO does not provide elevation at set
    passInfo.setAngles.rangeInKilometers = 0.0;  // N2YO does not provide range at set

    return passInfo;
}

void printPasses(std::vector<sattrack::PassInfo> passes) {
    using namespace sattrack;
    using namespace std::chrono;

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

        // Sort by riseTime (ascending)
        auto comparator = [](const sattrack::PassInfo& a, const sattrack::PassInfo& b) {
            return a.riseTime < b.riseTime;
        };
        std::sort(passes.begin(), passes.end(), comparator);

        std::cout << "Upcoming Passes:" << std::endl;
        std::cout << std::format(rowFormat, "Satellite", "Start", "End", "Duration", "Starts In", "Start Az", "End Az", "Max Az", "Max Elev") << std::endl;
        std::cout << std::format(rowFormat, sep25, sep25, sep25, sep14, sep18, sep12, sep12, sep12, sep10) << std::endl;
        for (const auto& pass : passes) {
            // Calculate all our times and durations
            auto startTime = pass.riseTime;
            auto endTime = pass.setTime;
            auto duration = duration_cast<seconds>(endTime - startTime);
            auto durationMins = duration_cast<minutes>(duration).count();
            auto durationSecs = duration_cast<seconds>(duration % minutes(1)).count();
            auto now = system_clock::now();
            auto eta = duration_cast<seconds>(startTime - now);
            auto etaHours = duration_cast<hours>(eta).count();
            auto etaMins = duration_cast<minutes>(eta % hours(1)).count();
            auto etaSecs = duration_cast<seconds>(eta % minutes(1)).count();

            // Output the pass data
            auto startTimeSeconds = std::chrono::floor<std::chrono::seconds>(startTime);
            auto endTimeSeconds = std::chrono::floor<std::chrono::seconds>(endTime);
            std::cout << std::format(rowFormat,
                std::format(satFormat, pass.noradID, pass.name),
                date::format(timeFormat, startTimeSeconds),
                date::format(timeFormat, endTimeSeconds),
                std::format(durationFormat, durationMins, durationSecs),
                std::format(etaFormat, etaHours, etaMins, etaSecs),
                std::format(azFormat, pass.riseAngles.azimuthInRadians * RADIANS_TO_DEGREES, azimuthToCompass(pass.riseAngles.azimuthInRadians)),
                std::format(azFormat, pass.setAngles.azimuthInRadians * RADIANS_TO_DEGREES, azimuthToCompass(pass.setAngles.azimuthInRadians)),
                std::format(azFormat, pass.maxAngles.azimuthInRadians * RADIANS_TO_DEGREES, azimuthToCompass(pass.maxAngles.azimuthInRadians)),
                std::format(elFormat, pass.maxAngles.elevationInRadians * RADIANS_TO_DEGREES)) << std::endl;
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
    config.setHorizon(0);

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
        "Filters out passes whose maximum elevation is below this number, in degrees (default 10)");
    n2yoPassesCommand->add_option("id", n2yoIDs, "Norad ID(s) of satellite(s)");

    auto satInfoCommand = app.add_subcommand("info", "View satellite information from the local TLE database");
    
    std::vector<int> satInfoIDs;
    satInfoCommand->add_option("id", satInfoIDs, "Norad ID(s) of satellite(s) (ie. 25544)");

    auto tleCommand = app.add_subcommand("tle", "View raw TLE data from the local TLE database");
    
    std::vector<int> satTLEIDs;
    tleCommand->add_option("id", satTLEIDs, "Norad ID(s) of satellite(s) (ie. 25544)");

    auto updateCommand = app.add_subcommand("update", "Update local TLE data from Celestrak");
    
    std::vector<std::string> groups;
    updateCommand->add_option<std::vector<std::string>>("group", groups, "Celestrak TLE group(s) to download (default: active)");

    auto geoCommand = app.add_subcommand("geo", "Get geodetic location (lat/long/alt) of the satellite at given time");

    std::vector<int> geoIDs;
    geoCommand->add_option("id", geoIDs, "Norad ID(s) of satellite(s) (ie. 25544)");

    // Look command - get current look angles for antenna pointing
    auto lookCommand = app.add_subcommand("look", "Get look angles (azimuth/elevation/range) for antenna pointing");

    std::vector<int> lookIDs;
    lookCommand->add_option("id", lookIDs, "Norad ID(s) of satellite(s) (ie. 25544)");
    lookCommand->add_option_function<std::string>("--time",
        [&config](const std::string &timeStr) {
            std::istringstream in(timeStr);
            std::chrono::system_clock::time_point tp;
            in >> date::parse("%Y-%m-%d %H:%M:%S", tp);
            if (in.fail()) {
                throw std::invalid_argument("Invalid time format (expected YYYY-MM-DD HH:MM:SS UTC): " + timeStr);
            }
            config.setTime(tp);
        }, "Time at which to get look angles (format: YYYY-MM-DD HH:MM:SS UTC)"
    );

    // Visible command - check if satellite is currently visible
    auto visibleCommand = app.add_subcommand("visible", "Check if satellite is visible from ground station");

    std::vector<int> visibleIDs;
    visibleCommand->add_option("id", visibleIDs, "Norad ID(s) of satellite(s) (ie. 25544)");
    visibleCommand->add_option_function<int>("--elev",
        [&config](const int elev) { config.setMinimumElevation(elev); },
        "Minimum elevation above the horizon in degrees (default 10)");
    visibleCommand->add_option_function<std::string>("--time",
        [&config](const std::string &timeStr) {
            std::istringstream in(timeStr);
            std::chrono::system_clock::time_point tp;
            in >> date::parse("%Y-%m-%d %H:%M:%S", tp);
            if (in.fail()) {
                throw std::invalid_argument("Invalid time format (expected YYYY-MM-DD HH:MM:SS UTC): " + timeStr);
            }
            config.setTime(tp);
        }, "Time at which to check visibility (format: YYYY-MM-DD HH:MM:SS UTC)"
    );

    // Track command - real-time tracking output
    auto trackCommand = app.add_subcommand("track", "Real-time tracking output (updates every interval)");

    std::vector<int> trackIDs;
    int trackInterval = 1;
    int trackDuration = 60;
    trackCommand->add_option("id", trackIDs, "Norad ID(s) of satellite(s) (ie. 25544)");
    trackCommand->add_option("--interval", trackInterval, "Update interval in seconds (default 1)");
    trackCommand->add_option("--duration", trackDuration, "Duration to track in seconds (default 60)");

    // Passes command - local pass prediction (no N2YO)
    auto passesCommand = app.add_subcommand("passes", "Predict satellite passes (local calculation)");

    std::vector<int> passesIDs;
    passesCommand->add_option("id", passesIDs, "Norad ID(s) of satellite(s) (ie. 25544)");
    passesCommand->add_option_function<int>("--days",
        [&config](const int days) { config.setDays(days); },
        "Number of days to search for passes (default 1)");
    passesCommand->add_option_function<int>("--elev",
        [&config](const int elev) { config.setMinimumElevation(elev); },
        "Filters out passes whose maximum elevation is below this number, in degrees (default 10)");
passesCommand->add_option_function<int>("--horizon",
        [&config](const int elev) { config.setHorizon(elev); },
        "Passes start and end when the satellite crosses this elevation, in degrees (default 0)");
        
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
        try {
            struct RadioPassWrapper {
                n2yo::SatelliteInfo info;
                n2yo::RadioPass pass;
            };

            std::vector<sattrack::PassInfo> passes;

            std::cerr << "Fetching Passes from N2YO..." << std::flush;
            for (auto noradID : n2yoIDs) {
                auto response = n2yo::getRadioPasses(config.getAPIKey(),
                    noradID, config.getLatitude(), config.getLongitude(), config.getAltitude(),
                    config.getDays(), config.getMinimumElevation(), config.getVerbose());
                for (auto pass: response.passes) {
                    passes.emplace_back(convertN2YOPassToPassInfo(pass, response.info));
                }
            }
            std::cerr << " Done." << std::endl << std::endl << std::flush;
            printPasses(passes);
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
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

    tleCommand->final_callback([tleCommand, &satTLEIDs, &tleFilename](void) {
        try {
            if (satTLEIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << tleCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Orbit> satellites;
            loadTLEDatabase(tleFilename, satellites);
            for (auto noradID : satTLEIDs) {
                if (!satellites.contains(noradID)) {
                    std::cerr << "Satellite with Norad ID " << noradID << " not found in the local TLE database." << std::endl;
                    continue;
                }
                auto orbit = satellites[noradID];
                std::cout << orbit.getTLE() << std::endl;
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
                std::cout << "Downloading TLE data for group: " << group << "..." << std::flush;
                auto response = celestrak::getTLE(group, config.getVerbose());
                for (auto tle : response.entries) {
                    sattrack::Orbit orbit;
                    orbit.updateFromTLE(tle);
                    satellites[orbit.getNoradID()] = orbit;
                }
                std::cout << " Done." << std::endl << std::flush;
            }

            saveTLEDatabase(tleFilename, satellites);
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    geoCommand->final_callback([geoCommand, &config, &tleFilename, &geoIDs](void) {
        using namespace sattrack;
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
                std::cout << "  Latitude: " << geo.latInRadians * RADIANS_TO_DEGREES << " deg" << std::endl;
                std::cout << "  Longitude: " << geo.lonInRadians * RADIANS_TO_DEGREES << " deg" << std::endl;
                std::cout << "  Altitude: " << geo.altInKilometers << " km" << std::endl;
                std::cout << std::endl;
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    lookCommand->final_callback([lookCommand, &config, &tleFilename, &lookIDs](void) {
        using namespace sattrack;
        try {
            if (lookIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << lookCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Orbit> satellites;
            loadTLEDatabase(tleFilename, satellites);

            // Create observer from config (convert degrees to radians, meters to km)
            sattrack::Geodetic observer{
                config.getLatitude() * DEGREES_TO_RADIANS,
                config.getLongitude() * DEGREES_TO_RADIANS,
                config.getAltitude() / 1000.0
            };

            double julianDate = sattrack::toJulianDate(config.getTime());

            for (auto noradID : lookIDs) {
                if (!satellites.contains(noradID)) {
                    std::cerr << "Satellite with Norad ID " << noradID << " not found in the local TLE database." << std::endl;
                    continue;
                }
                auto orbit = satellites[noradID];
                auto angles = sattrack::getLookAngles(orbit, observer, julianDate);

                double azDeg = angles.azimuthInRadians * RADIANS_TO_DEGREES ;
                double elDeg = angles.elevationInRadians * RADIANS_TO_DEGREES;

                std::cout << "Satellite: " << orbit.getName() << std::endl;
                std::cout << "  Azimuth:   " << std::format("{:6.2f}", azDeg) << " deg (" << azimuthToCompass(angles.azimuthInRadians) << ")" << std::endl;
                std::cout << "  Elevation: " << std::format("{:6.2f}", elDeg) << " deg" << std::endl;
                std::cout << "  Range:     " << std::format("{:6.1f}", angles.rangeInKilometers) << " km" << std::endl;
                std::cout << std::endl;
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    visibleCommand->final_callback([visibleCommand, &config, &tleFilename, &visibleIDs](void) {
        using namespace sattrack;
        try {
            if (visibleIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << visibleCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Orbit> satellites;
            loadTLEDatabase(tleFilename, satellites);

            sattrack::Geodetic observer{
                config.getLatitude() * DEGREES_TO_RADIANS,
                config.getLongitude() * DEGREES_TO_RADIANS,
                config.getAltitude() / 1000.0
            };

            double julianDate = sattrack::toJulianDate(config.getTime());
            double minElevRad = config.getMinimumElevation() * DEGREES_TO_RADIANS;

            for (auto noradID : visibleIDs) {
                if (!satellites.contains(noradID)) {
                    std::cerr << "Satellite with Norad ID " << noradID << " not found in the local TLE database." << std::endl;
                    continue;
                }
                auto orbit = satellites[noradID];
                auto angles = sattrack::getLookAngles(orbit, observer, julianDate);
                bool visible = sattrack::isVisible(angles, minElevRad);

                double azDeg = angles.azimuthInRadians * RADIANS_TO_DEGREES;
                double elDeg = angles.elevationInRadians * RADIANS_TO_DEGREES;

                std::cout << "Satellite: " << orbit.getName() << std::endl;
                std::cout << "  Visible:   " << (visible ? "YES" : "NO") << " (min elevation: " << config.getMinimumElevation() << " deg)" << std::endl;
                std::cout << "  Azimuth:   " << std::format("{:6.2f}", azDeg) << " deg (" << azimuthToCompass(angles.azimuthInRadians) << ")" << std::endl;
                std::cout << "  Elevation: " << std::format("{:6.2f}", elDeg) << " deg" << std::endl;
                std::cout << "  Range:     " << std::format("{:6.1f}", angles.rangeInKilometers) << " km" << std::endl;
                std::cout << std::endl;
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    trackCommand->final_callback([trackCommand, &config, &tleFilename, &trackIDs, &trackInterval, &trackDuration](void) {
        using namespace sattrack;
        try {
            if (trackIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << trackCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Orbit> satellites;
            loadTLEDatabase(tleFilename, satellites);

            sattrack::Geodetic observer{
                config.getLatitude() * DEGREES_TO_RADIANS,
                config.getLongitude() * DEGREES_TO_RADIANS,
                config.getAltitude() / 1000.0
            };

            // Verify all satellites exist
            for (auto noradID : trackIDs) {
                if (!satellites.contains(noradID)) {
                    std::cerr << "Satellite with Norad ID " << noradID << " not found in the local TLE database." << std::endl;
                    std::exit(1);
                }
            }

            constexpr std::string_view headerFormat = "{:<20} {:>10} {:>5} {:>12} {:>10}";
            constexpr std::string_view rowFormat = "{:<20} {:>7.2f} {:>3} {:>8.2f} {:>10.1f}";

            std::cout << std::format(headerFormat, "Satellite", "Azimuth", "", "Elevation", "Range (km)") << std::endl;
            std::cout << std::string(60, '-') << std::endl;

            auto startTime = std::chrono::system_clock::now();
            auto endTime = startTime + std::chrono::seconds(trackDuration);

            while (std::chrono::system_clock::now() < endTime) {
                auto now = std::chrono::system_clock::now();
                double julianDate = sattrack::toJulianDate(now);

                for (auto noradID : trackIDs) {
                    auto& orbit = satellites[noradID];
                    auto angles = sattrack::getLookAngles(orbit, observer, julianDate);

                    double azDeg = angles.azimuthInRadians * RADIANS_TO_DEGREES;
                    double elDeg = angles.elevationInRadians * RADIANS_TO_DEGREES;

                    std::cout << std::format(rowFormat,
                        orbit.getName().substr(0, 20),
                        azDeg,
                        azimuthToCompass(angles.azimuthInRadians),
                        elDeg,
                        angles.rangeInKilometers) << std::endl;
                }

                if (trackIDs.size() > 1) {
                    std::cout << std::endl;
                }

                std::this_thread::sleep_for(std::chrono::seconds(static_cast<long>(trackInterval)));
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    passesCommand->final_callback([passesCommand, &config, &tleFilename, &passesIDs](void) {
        using namespace sattrack;
        try {
            if (passesIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << passesCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Orbit> satellites;
            loadTLEDatabase(tleFilename, satellites);

            sattrack::Geodetic observer{
                config.getLatitude() * DEGREES_TO_RADIANS,
                config.getLongitude() * DEGREES_TO_RADIANS,
                config.getAltitude() / 1000.0
            };

            double minElevRad = config.getHorizon() * DEGREES_TO_RADIANS;
            auto searchDuration = std::chrono::hours(24 * config.getDays());

            std::vector<PassInfo> passes;

            for (auto noradID : passesIDs) {
                if (!satellites.contains(noradID)) {
                    std::cerr << "Satellite with Norad ID " << noradID << " not found in the local TLE database." << std::endl;
                    continue;
                }
                auto& orbit = satellites[noradID];

                std::cout << "Passes for " << orbit.getName() << ":" << std::endl;
                std::cout << std::endl;

                auto searchTime = std::chrono::system_clock::now();
                auto endTime = searchTime + searchDuration;
                int passCount = 0;

                while (searchTime < endTime) {
                    auto pass = sattrack::findNextPass(orbit, observer, minElevRad, searchTime);
                    if (!pass.has_value()) {
                        break;
                    }

                    if (pass.value().maxAngles.elevationInRadians * RADIANS_TO_DEGREES >= config.getMinimumElevation()) {
                        passCount++;
                        passes.push_back(*pass);
                    }

                    // Move past this pass
                    searchTime = pass->setTime + std::chrono::minutes(1);
                }

                if (passCount == 0) {
                    std::cerr << "No passes found for " << orbit.getName() << "(" << noradID << ")" << std::endl;
                }
            }

            printPasses(passes);

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
