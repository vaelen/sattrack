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

double Config::getAltitude() {
    return altitude;
}

void Config::setAltitude(const double a) {
    altitude = a;
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

int Config::getDays() {
    if (days > 0 && days <= 10) {
        return days;
    }
    if (days > 10) {
        return 10;
    }
    return 1;
}

void Config::setDays(const int d) {
    if (d > 0 && d <= 10) {
        days = d;
    } else if (d > 10) {
        days = 10;
    } else {
        days = 1;
    }
}

int Config::getMinimumElevation() {
    if (minimumElevation >= 0 && minimumElevation <= 90) {
        return minimumElevation;
    }
    if (minimumElevation > 90) {
        return 90;
    }
    return 0;
}

void Config::setMinimumElevation(const int degrees) {
    if (degrees >=0 && degrees <= 90) {
        minimumElevation = degrees;
    } else if (degrees > 90) {
        minimumElevation = 90;
    } else {
        minimumElevation = 0;
    }
}

bool Config::getVerbose() {
    return verbose;
}

void Config::setVerbose(bool v) {
    verbose = v;
}

}