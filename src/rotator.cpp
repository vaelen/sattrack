/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack/rotator.hpp>
#include <spdlog/spdlog.h>
#include <optional>
#include <mutex>

namespace sattrack {

using spdlog::debug;
using spdlog::info;

std::optional<double> Rotator::getAzimuth() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return azimuth_;
}

bool Rotator::hasValidAzimuth() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return azimuth_.has_value();
}

std::optional<double> Rotator::getElevation() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return elevation_;
}

bool Rotator::hasValidElevation() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return elevation_.has_value();
}

/** 
 * The default implementation parses the following formats:
 *      AZ=xxx EL=yyy
 *      +0000 +0000
 */
void Rotator::update(std::string &data) {
    double az = 0.0, el = 0.0;
    std::lock_guard<std::mutex> lock(mutex_);
    
    // Check for +0000 +0000 format first
    if (sscanf(data.c_str(), "%lf %lf", &az, &el) == 2) {
        azimuth_ = az;
        elevation_ = el;

        debug("Parsed rotator data in +0000 +0000 format: AZ={:.2f}, EL={:.2f}", az, el);

        return;
    }

    bool parsed = false;

    // Check for AZ=xxx EL=yyy format
    if (sscanf(data.c_str(), "AZ=%lf", &az) == 1) {
        azimuth_ = az;
        parsed = true;
        debug("Parsed rotator data in AZ=xxx format: AZ={:.2f}", az);
    }

    if (sscanf(data.c_str(), "EL=%lf", &el) == 1) {
        elevation_ = el;
        parsed = true;
        debug("Parsed rotator data in EL=yyy format: EL={:.2f}", el);
    }

    if (!parsed) {
        debug("Failed to parse rotator data: {}", data);
    }
}

std::string Rotator::getStatusCommand() const {
 return "C2\r\n";
}

} // namespace sattrack
