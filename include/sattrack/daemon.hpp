/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_DAEMON_HPP
#define __SATTRACK_DAEMON_HPP

#include <atomic>
#include <memory>
#include <mutex>
#include <thread>
#include <queue>
#include <condition_variable>
#include <boost/asio.hpp>
#include <sattrack/daemon/serialport.hpp>

namespace sattrack {

namespace asio = boost::asio;

using daemon::SerialPortOptions;

// The current status of the daemon
enum class DaemonStatus {
    STOPPED,
    STARTING,
    RUNNING,
    STOPPING
};

enum class DaemonEvent {
    STOP,
    RELOAD_CONFIG
};

class Daemon {
public:
    Daemon() = default;
    ~Daemon();

    DaemonStatus status();
    void start();
    void initSignals();
    void stop();
    void send(DaemonEvent event);
    void wait();

private:
    std::atomic<DaemonStatus> _status = DaemonStatus::STOPPED;
    std::thread eventLoopThread;
    std::mutex eventMutex;
    std::queue<DaemonEvent> eventQueue;
    std::condition_variable eventCV;
    std::thread ioThread;
    asio::io_context io;
    asio::signal_set signals{io, SIGINT, SIGTERM};
    std::unique_ptr<daemon::SerialPort> gpsSerialPort;
    std::unique_ptr<daemon::SerialPort> rotatorSerialPort;
    std::unique_ptr<daemon::SerialPort> radioSerialPort;

    void eventLoop();
};

}

#endif
