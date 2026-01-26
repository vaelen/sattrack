/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_DAEMON_SERIALPORT_HPP
#define __SATTRACK_DAEMON_SERIALPORT_HPP

#include <string>
#include <asio.hpp>

namespace sattrack::daemon {

struct SerialPortOptions {
    int baudRate = 9600;
    int characterSize = 8;
    asio::serial_port_base::parity::type parity = asio::serial_port_base::parity::none;
    asio::serial_port_base::stop_bits::type stopBits = asio::serial_port_base::stop_bits::one;
    asio::serial_port_base::flow_control::type flowControl = asio::serial_port_base::flow_control::none;
    std::string readTerminator = "\n";
};

class SerialPort {
public:
    SerialPort(asio::io_context& io, const std::string& name, const std::string& device, const SerialPortOptions& options);
    virtual ~SerialPort() = default;

    /** Opens the serial port and starts reading packets */
    void start();

    const std::string& name() const { return name_; }
    const std::string& device() const { return device_; }
    bool isOpen() const { return port_.is_open(); }

    static std::string ParityToString(asio::serial_port_base::parity::type parity);
    static std::string StopBitsToString(asio::serial_port_base::stop_bits::type stopBits);
    static std::string FlowControlToString(asio::serial_port_base::flow_control::type flowControl);

    /** Process incoming data from the serial port */
    virtual void processOutput(const char* data, std::size_t length);
    /** Send a command to the serial port */
    virtual void sendCommand(const std::string& command);

private:
    std::string name_;
    std::string device_;
    SerialPortOptions options_;
    asio::serial_port port_;
    asio::streambuf readBuffer_;

    virtual void readNextPacket();
};

}

#endif
