#include <sattrack/n2yo.hpp>
#include <curlpp/cURLpp.hpp>
#include <curlpp/Easy.hpp>
#include <curlpp/Options.hpp>
#include <rapidjson/document.h>
#include <sstream>
#include <stdexcept>

namespace sattrack {

std::string getTLE(const std::string& apiKey, int noradId) {
    // Initialize curlpp
    curlpp::Cleanup cleaner;
    curlpp::Easy request;

    // Build the API URL
    std::ostringstream urlBuilder;
    urlBuilder << "https://api.n2yo.com/rest/v1/satellite/tle/"
               << noradId << "&apiKey=" << apiKey;
    std::string url = urlBuilder.str();

    // Set up the request
    request.setOpt(new curlpp::options::Url(url));
    request.setOpt(new curlpp::options::FollowLocation(true));

    // Perform the request and capture response
    std::ostringstream responseStream;
    request.setOpt(new curlpp::options::WriteStream(&responseStream));

    try {
        request.perform();
    } catch (curlpp::RuntimeError& e) {
        throw std::runtime_error(std::string("HTTP request failed: ") + e.what());
    } catch (curlpp::LogicError& e) {
        throw std::runtime_error(std::string("HTTP logic error: ") + e.what());
    }

    // Parse the JSON response
    std::string responseStr = responseStream.str();
    rapidjson::Document doc;
    doc.Parse(responseStr.c_str());

    // Check for parse errors
    if (doc.HasParseError()) {
        throw std::runtime_error("Failed to parse JSON response from N2YO API");
    }

    // Check if the response contains an error
    if (doc.HasMember("error")) {
        std::string errorMsg = "N2YO API error: ";
        if (doc["error"].IsString()) {
            errorMsg += doc["error"].GetString();
        } else {
            errorMsg += "Unknown error";
        }
        throw std::runtime_error(errorMsg);
    }

    // Extract TLE data
    if (!doc.HasMember("tle") || !doc["tle"].IsString()) {
        throw std::runtime_error("Response does not contain TLE data");
    }

    return doc["tle"].GetString();
}

} // namespace sattrack