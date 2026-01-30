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
#include <asio.hpp>
#include <sattrack/serialport.hpp>
#include <sattrack/gps.hpp>
#include <sattrack/rotator.hpp>
#include <sattrack/config.hpp>
#include <sattrack/lcd.hpp>

namespace sattrack {

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
    Daemon(Config config) : config(std::move(config)) {}
    ~Daemon();

    DaemonStatus status();
    void start();
    void initSignals();
    void initLCD();
    void scheduleLCDUpdate();
    void scheduleStatusTimer();
    void stop();
    void send(DaemonEvent event);
    void wait();

private:
    Config config;
    std::atomic<DaemonStatus> _status = DaemonStatus::STOPPED;
    std::thread eventLoopThread;
    std::mutex eventMutex;
    std::queue<DaemonEvent> eventQueue;
    std::condition_variable eventCV;
    std::thread ioThread;
    asio::io_context io;
    asio::signal_set signals{io, SIGINT, SIGTERM};
    GPS gps;
    Rotator rotator;
    std::unique_ptr<GPSSerialPort> gpsSerialPort;
    std::unique_ptr<RotatorSerialPort> rotatorSerialPort;
    std::unique_ptr<RadioSerialPort> radioSerialPort;
    std::unique_ptr<asio::steady_timer> statusTimer;
    LCD lcd;
    std::unique_ptr<asio::steady_timer> lcdTimer;
    void eventLoop();
};

}

#endif
