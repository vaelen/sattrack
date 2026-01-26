/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack/daemon/serialport.hpp>
#include <spdlog/spdlog.h>

using spdlog::debug;
using spdlog::info;
using spdlog::warn;
using spdlog::error;

namespace sattrack::daemon {

using boost::system::error_code;

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

            processOutput(line.data(), len);

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

void SerialPort::processOutput(const char* data, std::size_t length) {
    // Placeholder for processing incoming data from the serial port
    std::string receivedData(data, length);
    debug("Received data on {} serial port: {}", name_, receivedData);

}

void SerialPort::sendCommand(const std::string& command) {
    if (!port_.is_open()) {
        throw std::runtime_error("Serial port not open: " + name_);
    }
    asio::async_write(port_, asio::buffer(command),
        [this](const error_code& ec, std::size_t /*bytes_transferred*/) {
            if (ec) {
                error("Failed to send command on {} serial port: {}", name_, ec.message());
            }
        });
}

}