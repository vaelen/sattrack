#include <sattrack.hpp>
#include <CLI/CLI.hpp>
#include <filesystem>
#include <fstream>

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

/** Program entry point */
int main(int argc, char* argv[]) {

    sattrack::Config config;
    config.setLatitude(0.0);
    config.setLongitude(0.0);
    config.setAltitude(0.0);
    config.setDays(1);
    config.setMinimumElevation(10);

    auto configFile = expandTilde("~/.sattrack.toml");

    CLI::App app{"SatTrack"};
    argv = app.ensure_utf8(argv);

    app.set_config("--config", configFile, "Read configuration from this file (default: " + configFile + ").");

    app.add_option_function<std::string>("--key",
        [&config](const std::string &key) { config.setAPIKey(key); },
        "N2YO API key");
    app.add_option_function<int>("--id", 
        [&config](const int id) { config.addSatellite(id); },
        "Norad ID of satellite to track (can be listed more than once)");
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

    auto passesCommand = app.add_subcommand("passes", "Display future passes");
    passesCommand->add_option_function<int>("--days", 
        [&config](const int days) { config.setDays(days); },
        "Number of days worth of passes to display (default 1, max 10)");
    passesCommand->add_option_function<int>("--elev", 
        [&config](const int elev) { config.setMinimumElevation(elev); },
        "Minimum pass elevation above the horizon in degrees (default 10, max 90)");

    app.parse_complete_callback(
        [&config](void) {
            if (!config.hasAPIKey()) {
                std::cerr << "Please provide your N2YO API key." << std::endl;
                std::exit(1);
            }

            if (!config.hasSatellites()) {
                std::cerr << "Please provide at least one satellite to track." << std::endl;
                std::exit(1);
            }
        }
    );

    tleCommand->final_callback([&config](void) {
        try {
            for (auto noradID : config.getSatellites()) {
                auto response = n2yo::getTLE(config.getAPIKey(), noradID, config.getVerbose());
                std::cout << "TLE for " << response.info.name << " (" << noradID << "):" << std::endl;
                std::cout << response.tle << std::endl << std::endl;
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    passesCommand->final_callback([&config](void) {
        try {
            for (auto noradID : config.getSatellites()) {
                auto response = n2yo::getRadioPasses(config.getAPIKey(), 
                    noradID, config.getLatitude(), config.getLongitude(), config.getAltitude(), 
                    config.getDays(), config.getMinimumElevation(), config.getVerbose());
                std::cout << "Upcoming Passes for " << response.info.name << " (" << noradID << "):" << std::endl;
                std::cout << "\tStart (UTC)\tEnd (UTC)\tStart Azimuth\tEnd Azimuth\tMax Azimuth\tMax Elevation" << std::endl;
                for (auto pass: response.passes) {
                    std::cout << '\t' << pass.startUTC;
                    std::cout << '\t' << pass.endUTC;
                    std::cout << '\t' << pass.startAz << "(" << pass.startAzCompass << ")";
                    std::cout << '\t' << pass.endAz << "(" << pass.endAzCompass << ")";
                    std::cout << '\t' << pass.maxAz << "(" << pass.maxAzCompass << ")";
                    std::cout << '\t' << pass.maxEl;
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }
        } catch (const std::exception &err) {
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
    });

    CLI11_PARSE(app, argc, argv);

    //return sattrack::startGUI(argc, argv);

    return 0;
}
