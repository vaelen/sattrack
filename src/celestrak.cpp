/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack/celestrak.hpp>
#include <curlpp/cURLpp.hpp>
#include <curlpp/Easy.hpp>
#include <curlpp/Options.hpp>
#include <sstream>
#include <stdexcept>

namespace celestrak {

// Celestrak GP data URL
// Documentation: https://celestrak.org/NORAD/documentation/gp-data-formats.php
#define BASE_URI "https://celestrak.org/NORAD/elements/gp.php"

std::string doGet(const std::string& url, bool debug = false) {
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

    std::string response = responseStream.str();

    if (debug) {
        std::cout << "Response length: " << response.length() << " bytes" << std::endl;
    }

    return response;
}

// Helper function to trim whitespace from both ends of a string
std::string trim(const std::string& str) {
    auto start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return "";
    }
    auto end = str.find_last_not_of(" \t\r\n");
    return str.substr(start, end - start + 1);
}

// Parse TLE format response into entries
// TLE format: Name line, Line 1, Line 2, repeated
std::vector<std::string> parseTLEResponse(const std::string& response) {
    std::vector<std::string> entries;
    std::istringstream stream(response);
    std::string line;
    std::ostringstream entryStream;

    // Read all non-empty lines
    int lineCount = 0;
    while (std::getline(stream, line)) {
        std::string trimmed = trim(line);
        if (!trimmed.empty()) {
            lineCount++;
            entryStream << trimmed << '\n';
            if (lineCount == 3) {
                entries.push_back(entryStream.str());
                entryStream.str("");
                entryStream.clear();
                lineCount = 0;
            }
        }
    }
    return entries;
}

TLEResponse getTLE(const std::string& group, bool debug) {
    // Use "active" as default if group is empty
    std::string groupName = group.empty() ? "active" : group;

    // Build the URL
    std::ostringstream urlBuilder;
    urlBuilder << BASE_URI << "?GROUP=" << groupName << "&FORMAT=tle";
    std::string url = urlBuilder.str();

    // Fetch the data
    std::string response = doGet(url, debug);

    // Check for error responses
    if (response.find("No GP data found") != std::string::npos) {
        throw std::runtime_error("Celestrak error: No GP data found for group '" + groupName + "'");
    }

    // Parse the TLE data
    std::vector<std::string> entries = parseTLEResponse(response);

    if (entries.empty() && !response.empty()) {
        throw std::runtime_error("Failed to parse TLE data from Celestrak response");
    }

    return TLEResponse{
        .group = groupName,
        .entries = std::move(entries)
    };
}

} // namespace celestrak
