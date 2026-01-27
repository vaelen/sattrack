/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_SERIALPORT_HPP
#define __SATTRACK_SERIALPORT_HPP

#include <string>
#include <asio.hpp>
#include <sattrack/gps.hpp>
#include <sattrack/rotator.hpp>

namespace sattrack {

constexpr double GPS_CHANGE_THRESHOLD_DEGREES = 0.0001;
constexpr double GPS_CHANGE_THRESHOLD_KM = 0.001;

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

    /** Execite any post-initalization code */
    virtual void started() {}

    /** Process incoming data from the serial port */
    virtual void processOutput(std::string &data);

    /** Send a command to the serial port */
    virtual void sendCommand(const std::string& command);

protected:
    std::string name_;
    std::string device_;
    SerialPortOptions options_;
    asio::serial_port port_;
    asio::streambuf readBuffer_;

    virtual void readNextPacket();
};

class GPSSerialPort : public SerialPort {
public:
    GPSSerialPort(asio::io_context& io, const std::string& name, const std::string& device, const SerialPortOptions& options, GPS& gps)
        : SerialPort(io, name, device, options), gps_(gps) {};
    virtual ~GPSSerialPort() = default;
    void processOutput(std::string &data) override;

protected:
    GPS& gps_;
};

class RotatorSerialPort : public SerialPort {
public:
    RotatorSerialPort(asio::io_context& io, const std::string& name, const std::string& device, const SerialPortOptions& options, Rotator& rotator);
    virtual ~RotatorSerialPort() = default;
    void started() override;
    void processOutput(std::string &data) override;
    void scheduleStatusCommandPoll();

protected:
    Rotator& rotator_;
    asio::steady_timer statusCommandTimer_;
    std::chrono::seconds statusCommandInterval_{5};
};

class RadioSerialPort : public SerialPort {
public:
    RadioSerialPort(asio::io_context& io, const std::string& name, const std::string& device, const SerialPortOptions& options)
        : SerialPort(io, name, device, options) {};

    virtual ~RadioSerialPort() = default;

};

} // namespace sattrack

#endif
