#include <sattrack/config.hpp>

namespace sattrack {

bool Config::hasAPIKey() {
    return apiKey.has_value();
}

void Config::clearAPIKey() {
    apiKey.reset();
}

std::string Config::getAPIKey() {
    return apiKey.value_or("");
}

void Config::setAPIKey(const std::string &newKey) {
    apiKey = newKey;
}

double Config::getLatitude() {
    return latitude;
}

void Config::setLatitude(const double l) {
    latitude = l;
}

double Config::getLongitude() {
    return longitude;
}

void Config::setLongitude(const double l) {
    longitude = l;
}

void Config::addSatellite(const int noradID) {
    satellites.insert(noradID);
}

template<typename Container>
std::enable_if<std::is_same_v<typename Container::value_type, int>> 
Config::addAllSatellites(const Container &noradIDs) {
    for (int noradID : noradIDs) {
        satellites.insert(noradID);
    }
}

void Config::removeSatellite(const int noradID) {
    satellites.erase(noradID);
}

template<typename Container>
std::enable_if<std::is_same_v<typename Container::value_type, int>> 
Config::removeAllSatellites(const Container &noradIDs) {
    for (int noradID : noradIDs) {
        satellites.erase(noradID);
    }
}

void Config::clearSatellites() {
    satellites.clear();
}

std::set<int> Config::getSatellites() {
    return satellites;
}

bool Config::hasSatellites() {
    return !satellites.empty();
}

}