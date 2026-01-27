/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack/serialport.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/chrono.h>
#include <chrono>
#include <asio.hpp>

using spdlog::debug;
using spdlog::info;
using spdlog::warn;
using spdlog::error;

namespace sattrack {

using asio::error_code;


///// SerialPort Implementation /////

SerialPort::SerialPort(asio::io_context& io, const std::string& name, const std::string& device, const SerialPortOptions& options)
    : name_(name)
    , device_(device)
    , options_(options)
    , port_(io) {}

void SerialPort::start() {
    error_code ec;

    if (port_.is_open()) {
        error_code ec;
        port_.close(ec);
        if (ec) {
            error("Failed to close {} serial port: {}", name_, ec.message());
        } else {
            info("{} serial port closed.", name_);
        }
    }

    port_.open(device_, ec);
    if (ec) {
        error("Failed to open {} serial port: {}", name_, ec.message());
        return;
    }

    port_.set_option(asio::serial_port_base::baud_rate(options_.baudRate), ec);
    if (ec) {
        error("Failed to set {} serial port baud rate: {}", name_, ec.message());
    } else {
        debug("{} serial port baud rate set to {}bps.", name_, options_.baudRate);
    }

    port_.set_option(asio::serial_port_base::character_size(options_.characterSize), ec);
    if (ec) {
        error("Failed to set {} serial port character size: {}", name_, ec.message());
    } else {
        debug("{} serial port character size set to {} bits.", name_, options_.characterSize);
    }

    port_.set_option(asio::serial_port_base::parity(options_.parity), ec);
    if (ec) {
        error("Failed to set {} serial port parity: {}", name_, ec.message());
    } else {

        debug("{} serial port parity set to {}.", name_, ParityToString(options_.parity));
    }

    port_.set_option(asio::serial_port_base::stop_bits(options_.stopBits), ec);
    if (ec) {
        error("Failed to set {} serial port stop bits: {}", name_, ec.message());
    } else {
        debug("{} serial port stop bits set to {}.", name_, StopBitsToString(options_.stopBits));
    }

    port_.set_option(asio::serial_port_base::flow_control(options_.flowControl), ec);
    if (ec) {
        error("Failed to set {} serial port flow control: {}", name_, ec.message());
    } else {
        debug("{} serial port flow control set to {}.", name_, FlowControlToString(options_.flowControl));
    }

    info("{} serial port initialized.", name_);

    readNextPacket();
    started();
}

void SerialPort::readNextPacket() {

    if (!port_.is_open()) {
        error("{} serial port is not open for reading.", name_);
        return;
    }

    asio::async_read_until(port_, readBuffer_, options_.readTerminator,
        [this](const error_code& ec, std::size_t bytes_transferred) {
            if (ec) {
                error("{} serial port read error: {}", name_, ec.message());
                return;
            }

            // Process the received data
            auto bufs = readBuffer_.data();
            std::size_t len = bytes_transferred - options_.readTerminator.length();
            std::string line(asio::buffers_begin(bufs), asio::buffers_begin(bufs) + len);
            readBuffer_.consume(bytes_transferred);

            debug("Received on {} serial port: {}", name_, line);

            processOutput(line);

            // Continue reading
            readNextPacket();
        });
}

std::string SerialPort::ParityToString(asio::serial_port_base::parity::type parity) {
    switch (parity) {
        case asio::serial_port_base::parity::none:
            return "None";
        case asio::serial_port_base::parity::odd:
            return "Odd";
        case asio::serial_port_base::parity::even:
            return "Even";
        default:
            return "Unknown";
    }
}

std::string SerialPort::StopBitsToString(asio::serial_port_base::stop_bits::type stopBits) {
    switch (stopBits) {
        case asio::serial_port_base::stop_bits::one:
            return "1";
        case asio::serial_port_base::stop_bits::onepointfive:
            return "1.5";
        case asio::serial_port_base::stop_bits::two:
            return "2";
        default:
            return "?";
    }

}

std::string SerialPort::FlowControlToString(asio::serial_port_base::flow_control::type flowControl) {
    switch (flowControl) {
        case asio::serial_port_base::flow_control::none:
            return "None";
        case asio::serial_port_base::flow_control::software:
            return "Software";
        case asio::serial_port_base::flow_control::hardware:
            return "Hardware";
        default:
            return "Unknown";
    }
}

void SerialPort::processOutput(std::string &data) {
    // Placeholder for processing incoming data from the serial port
    info("Received data on {} serial port: {}", name_, data);
}

void SerialPort::sendCommand(const std::string& command) {
    if (!port_.is_open()) {
        throw std::runtime_error("Serial port not open: " + name_);
    }
    debug("Sending command on {} serial port: {}", name_, command);
    asio::async_write(port_, asio::buffer(command),
        [this](const error_code& ec, std::size_t /*bytes_transferred*/) {
            if (ec) {
                error("Failed to send command on {} serial port: {}", name_, ec.message());
            }
        });
}

///// GPSerialPort Implementation /////

void GPSSerialPort::processOutput(std::string &data) {
    auto oldPosition = gps_.getPosition();
    auto oldTime = gps_.getUTCTime();

    gps_.update(data);

    auto newPosition = gps_.getPosition();
    auto newTime = gps_.getUTCTime();

    if (newPosition.has_value()) {
        bool positionChanged = false;
        auto newPos = newPosition.value();
        auto newLat = newPos.latInRadians * RADIANS_TO_DEGREES;
        auto newLon = newPos.lonInRadians * RADIANS_TO_DEGREES;
        auto newAlt = newPos.altInKilometers;

        if (!oldPosition.has_value()) {
            positionChanged = true;
        } else {
            auto oldPos = oldPosition.value();
            auto oldLat = oldPos.latInRadians * RADIANS_TO_DEGREES;
            auto oldLon = oldPos.lonInRadians * RADIANS_TO_DEGREES;
            auto oldAlt = oldPos.altInKilometers;

            positionChanged = (std::abs(oldLat - newLat) >= GPS_CHANGE_THRESHOLD_DEGREES) ||
                              (std::abs(oldLon - newLon) >= GPS_CHANGE_THRESHOLD_DEGREES) ||
                              (std::abs(oldAlt - newAlt) >= GPS_CHANGE_THRESHOLD_KM);
        }

        if (positionChanged) {
            info("GPS Position {}: Lat {:.6f}°, Lon {:.6f}°, Alt {:.2f} km",
                 oldPosition.has_value() ? "Updated" : "Acquired", newLat, newLon, newAlt);
        }

    } else if (oldPosition.has_value() && !newPosition.has_value()) {
        warn("GPS Position Lost");
    }

    if (!oldTime.has_value() && newTime.has_value()) {
        info("GPS Time Acquired: {:%Y-%m-%d %H:%M:%S} UTC", newTime.value());
    } else if (oldTime.has_value() && !newTime.has_value()) {
        warn("GPS Time Lost");
    }

}

///// RotatorSerialPort Implementation /////

RotatorSerialPort::RotatorSerialPort(asio::io_context& io, const std::string& name, const std::string& device, const SerialPortOptions& options, Rotator& rotator)
    : SerialPort(io, name, device, options), rotator_(rotator), statusCommandTimer_(io) {
}

void RotatorSerialPort::started() {
    scheduleStatusCommandPoll();
}

void RotatorSerialPort::scheduleStatusCommandPoll() {
    statusCommandTimer_.expires_after(statusCommandInterval_);
    statusCommandTimer_.async_wait([this](const error_code& ec) {
        if (!ec && port_.is_open()) {
            sendCommand(rotator_.getStatusCommand());
            scheduleStatusCommandPoll();
        }
    });
}

void RotatorSerialPort::processOutput(std::string &data) {
    debug("Rotator received data: {}", data);

    auto oldAz = rotator_.getAzimuth();
    auto oldEl = rotator_.getElevation();

    rotator_.update(data);

    auto newAz = rotator_.getAzimuth();
    auto newEl = rotator_.getElevation();

    if (oldAz != newAz || oldEl != newEl) {
        info("Rotator Position Updated: Az {:.2f}°, El {:.2f}°", newAz.value_or(0.0), newEl.value_or(0.0));
    }
}


} // namespace sattrack
