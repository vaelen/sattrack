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

    std::string configFile = expandTilde("~/.sattrack.toml");

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

    CLI11_PARSE(app, argc, argv);

    if (!config.hasAPIKey()) {
        std::cerr << "Please provide your N2YO API key." << std::endl;
        return 1;
    }

    if (!config.hasSatellites()) {
        std::cerr << "Please provide at least one satellite to track." << std::endl;
        return 1;
    }

    try {
        for (auto noradID : config.getSatellites()) {
            auto response = n2yo::getTLE(config.getAPIKey(), noradID);
            std::cout << "TLE for " << response.info.name << " (" << noradID << "):" << std::endl;
            std::cout << response.tle << std::endl << std::endl;
        }
    } catch (const std::exception &err) {
        std::cerr << err.what() << std::endl;
        return 1;
    }

    //return sattrack::startGUI(argc, argv);

    return 0;
}
