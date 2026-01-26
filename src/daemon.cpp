/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack/daemon.hpp>
#include <spdlog/spdlog.h>

using spdlog::debug;
using spdlog::info;
using spdlog::warn;
using spdlog::error;
using spdlog::critical;

namespace sattrack {

Daemon::~Daemon() {
    stop();
    if (eventLoopThread.joinable()) {
        info("Waiting for event loop thread to finish...");
        eventLoopThread.join();
    }
    if (ioThread.joinable()) {
        info("Waiting for IO thread to finish...");
        ioThread.join();
    }
    info("Daemon destroyed.");
}

DaemonStatus Daemon::status() {
    return _status.load();
}

void Daemon::start() {
    DaemonStatus expected = DaemonStatus::STOPPED;
    if (!_status.compare_exchange_strong(expected, DaemonStatus::STARTING)) {
        return;
    }
    info("Starting daemon...");
    eventLoopThread = std::thread([this]{ eventLoop(); });

    initSignals();

    // TODO: Add serial port configuration options

    SerialPortOptions gpsPortOptions;
    gpsPortOptions.readTerminator = "\r\n";
    gpsSerialPort = std::make_unique<daemon::SerialPort>(io, "GPS", "/dev/ttyS1", gpsPortOptions);
    gpsSerialPort->start();

    SerialPortOptions rotatorPortOptions;
    rotatorPortOptions.readTerminator = "\r\n";
    rotatorSerialPort = std::make_unique<daemon::SerialPort>(io, "Rotator", "/dev/ttyS2", rotatorPortOptions);
    rotatorSerialPort->start();

    SerialPortOptions radioPortOptions;
    radioPortOptions.readTerminator = "\r\n";
    radioPortOptions.baudRate = 38400;
    radioSerialPort = std::make_unique<daemon::SerialPort>(io, "Radio", "/dev/ttyS4", radioPortOptions);
    radioSerialPort->start();

    ioThread = std::thread([this]{ io.run(); });
}

void Daemon::initSignals() {
    signals.async_wait([this](auto ec, int sig) {
        if (ec) {
            error("Error receiving signal: {}", ec.message());
        } else {
            info("Received signal {}.", sig);
        }
        stop(); 
    });
}

void Daemon::stop() {
    DaemonStatus expected = DaemonStatus::RUNNING;
    if (!_status.compare_exchange_strong(expected, DaemonStatus::STOPPING)) {
        return;
    }
    info("Stopping daemon...");
    eventCV.notify_all();  // Wake up event loop so it sees the status change
    io.stop();
}

void Daemon::send(DaemonEvent event) {
    {
        std::scoped_lock lock(eventMutex);
        eventQueue.push(event);
    }
    eventCV.notify_one();
}

void Daemon::eventLoop() {
    // Update status to RUNNING
    _status.store(DaemonStatus::RUNNING);
    info("Daemon started.");

    while (status() == DaemonStatus::RUNNING) {
        std::unique_lock eventLock(eventMutex);
        eventCV.wait(eventLock, [this]{
            return !eventQueue.empty() || status() != DaemonStatus::RUNNING;
        });
        if (status() != DaemonStatus::RUNNING) {
            break;
        }
        while (!eventQueue.empty()) {
            DaemonEvent event = eventQueue.front();
            eventQueue.pop();
            eventLock.unlock();

            switch (event) {
                case DaemonEvent::STOP:
                    info("Received STOP event.");
                    stop();
                    break;
                case DaemonEvent::RELOAD_CONFIG:
                    info("Received RELOAD_CONFIG event.");
                    break;
            }

            eventLock.lock();
        }
    }

    // Update status to STOPPED when exiting loop
    _status.store(DaemonStatus::STOPPED);
    info("Daemon stopped.");
    
}

void Daemon::wait() {
    if (eventLoopThread.joinable()) {
        eventLoopThread.join();
    }
    if (ioThread.joinable()) {
        ioThread.join();
    }
}

}
