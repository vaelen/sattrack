/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack/n2yo.hpp>
#include <curlpp/cURLpp.hpp>
#include <curlpp/Easy.hpp>
#include <curlpp/Options.hpp>
#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <sstream>
#include <stdexcept>
#include <iomanip>

namespace n2yo {

// API Documentation: https://www.n2yo.com/api/

#define BASE_URI "https://api.n2yo.com/rest/v1/satellite/"

void printDocument(rapidjson::Document &doc) {
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    doc.Accept(writer);

    std::cout << buffer.GetString() << std::endl;
}

rapidjson::Document doGet(std::string url, bool debug = false) {
    // Initialize curlpp
    curlpp::Cleanup cleaner;
    curlpp::Easy request;

    if (debug) {
        std::cout << url << std::endl;
    }

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

    if (debug) {
        printDocument(doc);
    }

    return doc;
}

TLEResponse getTLE(const std::string &apiKey, int noradId, bool debug) {
    // Build the API URL
    std::ostringstream urlBuilder;
    urlBuilder << BASE_URI << "tle/"
               << noradId << "&apiKey=" << apiKey;
    std::string url = urlBuilder.str();

    auto doc = doGet(url, debug);

    // Extract TLE data
    if (!doc.HasMember("tle") || !doc["tle"].IsString()) {
        throw std::runtime_error("Response does not contain TLE data");
    }

    auto info = doc["info"].GetObject();

    TLEResponse response {
        .info = {
            .id = info["satid"].GetInt(),
            .name = info["satname"].GetString()
        },
        .tle = doc["tle"].GetString()
    };

    return response;
}

PositionsResponse getPositions(const std::string& apiKey, int noradId,
                               double observerLat, double observerLng,
                               double observerAlt, int seconds, 
                               bool debug) {
    // Build the API URL
    std::ostringstream urlBuilder;
    urlBuilder << std::fixed << std::setprecision(5);
    urlBuilder << BASE_URI << "positions/"
               << noradId << "/"
               << observerLat << "/"
               << observerLng << "/"
               << observerAlt << "/"
               << seconds << "/"
               << "&apiKey=" << apiKey;
    std::string url = urlBuilder.str();

    auto doc = doGet(url, debug);

    // Extract satellite info
    if (!doc.HasMember("info") || !doc["info"].IsObject()) {
        throw std::runtime_error("Response does not contain satellite info");
    }

    auto info = doc["info"].GetObject();

    // Extract positions array
    if (!doc.HasMember("positions") || !doc["positions"].IsArray()) {
        throw std::runtime_error("Response does not contain positions data");
    }

    auto positionsArray = doc["positions"].GetArray();
    std::vector<Position> positions;
    positions.reserve(positionsArray.Size());

    for (const auto& pos : positionsArray) {
        positions.push_back({
            .latitude = pos["satlatitude"].GetDouble(),
            .longitude = pos["satlongitude"].GetDouble(),
            .altitude = pos["sataltitude"].GetDouble(),
            .azimuth = pos["azimuth"].GetDouble(),
            .elevation = pos["elevation"].GetDouble(),
            .ra = pos["ra"].GetDouble(),
            .dec = pos["dec"].GetDouble(),
            .timestamp = pos["timestamp"].GetInt64()
        });
    }

    PositionsResponse response{
        .info = {
            .id = info["satid"].GetInt(),
            .name = info["satname"].GetString()
        },
        .positions = std::move(positions)
    };

    return response;
}

RadioPassesResponse getRadioPasses(const std::string& apiKey, int noradId,
                                   double observerLat, double observerLng,
                                   double observerAlt, int days, int minElevation, 
                                   bool debug) {
    // Build the API URL
    std::ostringstream urlBuilder;
    urlBuilder << std::fixed << std::setprecision(5);
    urlBuilder << BASE_URI << "radiopasses/"
               << noradId << "/"
               << observerLat << "/"
               << observerLng << "/"
               << observerAlt << "/"
               << days << "/"
               << minElevation << "/"
               << "&apiKey=" << apiKey;
    std::string url = urlBuilder.str();

    auto doc = doGet(url, debug);

    // Extract satellite info
    if (!doc.HasMember("info") || !doc["info"].IsObject()) {
        throw std::runtime_error("Response does not contain satellite info");
    }

    auto info = doc["info"].GetObject();

    // Extract passes array
    if (!doc.HasMember("passes") || !doc["passes"].IsArray()) {
        throw std::runtime_error("No passes found. Check longitude and latitude.");
    }

    auto passesArray = doc["passes"].GetArray();
    std::vector<RadioPass> passes;
    passes.reserve(passesArray.Size());

    for (const auto& pass : passesArray) {
        passes.push_back({
            .startAz = pass["startAz"].GetDouble(),
            .startAzCompass = pass["startAzCompass"].GetString(),
            .startUTC = pass["startUTC"].GetInt64(),
            .maxAz = pass["maxAz"].GetDouble(),
            .maxAzCompass = pass["maxAzCompass"].GetString(),
            .maxEl = pass["maxEl"].GetDouble(),
            .maxUTC = pass["maxUTC"].GetInt64(),
            .endAz = pass["endAz"].GetDouble(),
            .endAzCompass = pass["endAzCompass"].GetString(),
            .endUTC = pass["endUTC"].GetInt64()
        });
    }

    RadioPassesResponse response{
        .info = {
            .id = info["satid"].GetInt(),
            .name = info["satname"].GetString()
        },
        .passes = std::move(passes)
    };

    return response;
}

} // namespace n2yo