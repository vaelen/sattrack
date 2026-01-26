/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack.hpp>
#include <CLI/CLI.hpp>
#include <date/date.h>
#include <filesystem>
#include <fstream>
#include <spdlog/fmt/fmt.h>
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
        std::cout << fmt::format(rowFormat, "Satellite", "Start", "End", "Duration", "Starts In", "Start Az", "End Az", "Max Az", "Max Elev") << std::endl;
        std::cout << fmt::format(rowFormat, sep25, sep25, sep25, sep14, sep18, sep12, sep12, sep12, sep10) << std::endl;
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
            std::cout << fmt::format(rowFormat,
                fmt::format(satFormat, pass.noradID, pass.name),
                date::format(timeFormat, startTimeSeconds),
                date::format(timeFormat, endTimeSeconds),
                fmt::format(durationFormat, durationMins, durationSecs),
                fmt::format(etaFormat, etaHours, etaMins, etaSecs),
                fmt::format(azFormat, pass.riseAngles.azimuthInRadians * RADIANS_TO_DEGREES, azimuthToCompass(pass.riseAngles.azimuthInRadians)),
                fmt::format(azFormat, pass.setAngles.azimuthInRadians * RADIANS_TO_DEGREES, azimuthToCompass(pass.setAngles.azimuthInRadians)),
                fmt::format(azFormat, pass.maxAngles.azimuthInRadians * RADIANS_TO_DEGREES, azimuthToCompass(pass.maxAngles.azimuthInRadians)),
                fmt::format(elFormat, pass.maxAngles.elevationInRadians * RADIANS_TO_DEGREES)) << std::endl;
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

    ///// Global options /////

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

    ///// Satellite Info Command - view satellite information from local TLE database /////

    auto satInfoCommand = app.add_subcommand("info", "View satellite information from the local TLE database");
    
    std::vector<int> satInfoIDs;
    satInfoCommand->add_option("id", satInfoIDs, "Norad ID(s) of satellite(s) (ie. 25544)");

    satInfoCommand->final_callback([satInfoCommand, &satInfoIDs, &tleFilename](void) {
        try {
            if (satInfoIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << satInfoCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Satellite> satellites;
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

    ///// TLE Command - view raw TLE data from local TLE database /////

    auto tleCommand = app.add_subcommand("tle", "View raw TLE data from the local TLE database");
    
    std::vector<int> satTLEIDs;
    tleCommand->add_option("id", satTLEIDs, "Norad ID(s) of satellite(s) (ie. 25544)");

    tleCommand->final_callback([tleCommand, &satTLEIDs, &tleFilename](void) {
        try {
            if (satTLEIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << tleCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Satellite> satellites;
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

    ///// Update command - fetch latest TLE data from Celestrak /////

    auto updateCommand = app.add_subcommand("update", "Update local TLE data from Celestrak");
    
    std::vector<std::string> groups;
    updateCommand->add_option<std::vector<std::string>>("group", groups, "Celestrak TLE group(s) to download (default: active)");

    updateCommand->final_callback([&config, &groups, &tleFilename](void) {
        try {
            std::map<int, sattrack::Satellite> satellites;
            loadTLEDatabase(tleFilename, satellites);

            if (groups.empty()) {
                groups.push_back("active");
            }

            for (auto group : groups) {
                std::cout << "Downloading TLE data for group: " << group << "..." << std::flush;
                auto response = celestrak::getTLE(group, config.getVerbose());
                for (auto tle : response.entries) {
                    sattrack::Satellite orbit;
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

    ///// Geo Command - get geodetic location of satellite /////

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

    geoCommand->final_callback([geoCommand, &config, &tleFilename, &geoIDs](void) {
        using namespace sattrack;
        try {
            if (geoIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << geoCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Satellite> satellites;
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

    ///// Look command - get current look angles for antenna pointing /////
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

    lookCommand->final_callback([lookCommand, &config, &tleFilename, &lookIDs](void) {
        using namespace sattrack;
        try {
            if (lookIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << lookCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Satellite> satellites;
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
                std::cout << "  Azimuth:   " << fmt::format("{:6.2f}", azDeg) << " deg (" << azimuthToCompass(angles.azimuthInRadians) << ")" << std::endl;
                std::cout << "  Elevation: " << fmt::format("{:6.2f}", elDeg) << " deg" << std::endl;
                std::cout << "  Range:     " << fmt::format("{:6.1f}", angles.rangeInKilometers) << " km" << std::endl;
                std::cout << std::endl;
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    ///// Visible command - check if satellite is currently visible from ground station /////
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

    visibleCommand->final_callback([visibleCommand, &config, &tleFilename, &visibleIDs](void) {
        using namespace sattrack;
        try {
            if (visibleIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << visibleCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Satellite> satellites;
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
                std::cout << "  Azimuth:   " << fmt::format("{:6.2f}", azDeg) << " deg (" << azimuthToCompass(angles.azimuthInRadians) << ")" << std::endl;
                std::cout << "  Elevation: " << fmt::format("{:6.2f}", elDeg) << " deg" << std::endl;
                std::cout << "  Range:     " << fmt::format("{:6.1f}", angles.rangeInKilometers) << " km" << std::endl;
                std::cout << std::endl;
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    ///// Passes command - local pass prediction /////
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

    passesCommand->final_callback([passesCommand, &config, &tleFilename, &passesIDs](void) {
        using namespace sattrack;
        try {
            if (passesIDs.empty()) {
                std::cerr << "Please provide at least one satellite's Norad ID." << std::endl;
                std::cerr << passesCommand->help() << std::endl;
                std::exit(1);
            }
            std::map<int, sattrack::Satellite> satellites;
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

    ///// Daemon command - run sattrack as a background daemon service /////
    auto daemonCommand = app.add_subcommand("daemon", "Run SatTrack as a background daemon service");

    std::vector<int> daemonNoradIDs;
    daemonCommand->add_option("id", daemonNoradIDs, "Norad ID(s) of satellite(s) to track (ie. 25544)");    

    daemonCommand->add_option_function<std::string>("--radio",
        [&config](const std::string& radioInfoString) {
            // Parse radio info string0
            sattrack::Radio radio = sattrack::parseRadioInfo(radioInfoString);
            // Store the radio info somewhere accessible to the daemon
            // For this example, we'll just print it
            std::cout << "Configured radio for Norad ID " << radio.noradID << ":" << std::endl;
            std::cout << "  Type: " << (radio.type == sattrack::RadioType::UPLINK ? "Uplink" : "Downlink") << std::endl;
            std::cout << "  Frequency: " << radio.baseFrequencyInKHz << " kHz" << std::endl;
            std::cout << "  Mode: " << radio.mode << std::endl;
            std::cout << "  Modulation: " << radio.modulation << std::endl;
            std::cout << "  Packeting: " << radio.packeting << std::endl;
         }, 
        "Radio configuration in the format NORAD_ID:TYPE:FREQUENCY_KHZ:MODE:MODULATION:PACKETING "
        "(ie. 25544:DOWNLINK:145800:FM:GMSK9600:APRS)");

    daemonCommand->add_option_function<int>("--horizon",
        [&config](const int elev) { config.setHorizon(elev); },
        "Passes start and end when the satellite crosses this elevation, in degrees (default 0)");
    daemonCommand->add_option_function<int>("--elev",
        [&config](const int elev) { config.setMinimumElevation(elev); },
        "Filters out passes whose maximum elevation is below this number, in degrees (default 10)");

    daemonCommand->final_callback([](void) {
        try {
            sattrack::Daemon daemon;
            daemon.start();
            daemon.wait();
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    ///// End of command definitions /////

    CLI11_PARSE(app, argc, argv);

    if (app.get_subcommands().empty()) {
        std::cerr << app.help() << std::endl;
        std::exit(1);
    }

    return 0;
}
