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

    CLI::App app{"SatTrack"};
    argv = app.ensure_utf8(argv);

    std::string filename = "default";
    app.add_option("-c,--config", configFile, "Read configuration from this file (default: " + configFile + ").");

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


    return sattrack::startGUI(argc, argv);

}
