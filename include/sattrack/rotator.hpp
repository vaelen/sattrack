/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_ROTATOR_HPP
#define __SATTRACK_ROTATOR_HPP

#include <mutex>
#include <string>
#include <optional>

namespace sattrack {

/**
 * Rotator class for tracking antenna rotator position.
 *
 * Stores the current azimuth and elevation of the rotator in degrees.
 * All access is thread-safe via mutex protection.
 * 
 * The default implementation can control rotators that use the Yaesu GS-232
 * command set, responding to "C2" status requests with "AZ=xxx EL=yyy" responses.
 */
class Rotator {
public:
    Rotator() = default;
    ~Rotator() = default;

    // Non-copyable, non-movable (due to mutex)
    Rotator(const Rotator&) = delete;
    Rotator& operator=(const Rotator&) = delete;
    Rotator(Rotator&&) = delete;
    Rotator& operator=(Rotator&&) = delete;

    /**
     * Get the current azimuth in degrees.
     * @return Azimuth in degrees (0-360)
     */
    std::optional<double> getAzimuth() const;

    /** 
     * Check if we have a valid azimuth reading.
     * @return true if azimuth is valid, false otherwise
     */
    bool hasValidAzimuth() const;

    /**
     * Get the current elevation in degrees.
     * @return Elevation in degrees (typically -90 to 90)
     */
    std::optional<double> getElevation() const;

    /** 
     * Check if we have a valid elevation reading.
     * @return true if elevation is valid, false otherwise
     */
    bool hasValidElevation() const;

    /**
     * Process a status update message from the rotator.
     * @param data Status update message
     */
    void update(std::string &data);

    /**
     * Get the command string to request status from the rotator.
     * @return Status command string
     */
    std::string getStatusCommand() const;

private:
    mutable std::mutex mutex_;
    std::optional<double> azimuth_;
    std::optional<double> elevation_;
    bool moving_ = false;
};

} // namespace sattrack

#endif
