#include <sattrack.hpp>
#include <CLI/CLI.hpp>
#include <toml++/toml.hpp>
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

void writeConfig(const std::string &configPath, const toml::table &config) {
    std::ofstream file(configPath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not create confg file: " + configPath);
    }

    file << config << std::endl;
}

/** Creates a default configuration file. */
toml::table createDefaultConfig(const std::string &configPath) {
    // Create default configuration
    toml::table config;
    config.insert("api_key", "");
    config.insert("latitude", 0.0);
    config.insert("longitude", 0.0);

    // Create parent directories if needed
    fs::path path(configPath);
    if (path.has_parent_path()) {
        fs::create_directories(path.parent_path());
    }

    // Write configuration to file
    writeConfig(configPath, config);

    std::cout << "Created default configuration file: " << configPath << std::endl;

    return config;
}

/** Program entry point */
int main(int argc, char* argv[]) {

    std::string configFile = "~/.sattrack.toml";
    std::optional<std::string> apiKey;

    CLI::App app{"SatTrack"};
    argv = app.ensure_utf8(argv);

    app.add_option("-c,--config", configFile, "Read configuration from this file (default: " + configFile + ").");
    app.add_option("-k,--key", apiKey, "N2YO API key");

    CLI11_PARSE(app, argc, argv);

    configFile = expandTilde(configFile);

    std::error_code ec;
    bool exists = fs::exists(configFile, ec);

    toml::table config;
    if (!ec && exists) {
        // Load settings
        try {
            config = toml::parse_file(configFile);
            std::cout << "Loaded config file: " << configFile << std::endl;
        } catch (const toml::parse_error &err) {
            std::cerr << "Error loading config file:" << std::endl << err << std::endl;
            return 1;
        }
    } else if (ec) {
        // An error occurrect
        std::cerr << "Error: " << ec.message() << std::endl;
        return 1;
    } else {
        // File doesn't exist, create an empty file
        config = createDefaultConfig(configFile);
    }

    if (!apiKey.has_value()) {
        std::optional<std::string> apiKey = config["api_key"].value<std::string>();
    }
    if (!apiKey.has_value()) {
        std::cerr << "Please provide your N2YO API key." << std::endl;
        return 1;
    }

    try {
        std::string tle = sattrack::getTLE(apiKey.value(), 25544);
        std::cout << "TLE:" << std::endl;
        std::cout << tle << std::endl;
    } catch (const std::exception &err) {
        std::cerr << err.what() << std::endl;
        return 1;
    }

    //return sattrack::startGUI(argc, argv);

    return 0;
}
