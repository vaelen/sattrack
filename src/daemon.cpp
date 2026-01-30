/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack/daemon.hpp>
#include <spdlog/spdlog.h>
#include <chrono>
#include <spdlog/fmt/chrono.h>

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

    initLCD();

    initSignals();

    // TODO: Add serial port configuration options

    SerialPortOptions gpsPortOptions;
    gpsPortOptions.readTerminator = "\r\n";
    gpsPortOptions.baudRate = config.getGPSBaudRate();
    gpsSerialPort = std::make_unique<GPSSerialPort>(io, "GPS", config.getGPSSerialPort(), gpsPortOptions, gps);
    gpsSerialPort->start();

    SerialPortOptions rotatorPortOptions;
    rotatorPortOptions.readTerminator = "\r\n";
    rotatorPortOptions.baudRate = config.getRotatorBaudRate();
    rotatorSerialPort = std::make_unique<RotatorSerialPort>(io, "Rotator", config.getRotatorSerialPort(), rotatorPortOptions, rotator);
    rotatorSerialPort->start();

    SerialPortOptions radioPortOptions;
    radioPortOptions.readTerminator = "\r\n";
    radioPortOptions.baudRate = config.getRadioBaudRate();
    radioSerialPort = std::make_unique<RadioSerialPort>(io, "Radio", config.getRadioSerialPort(), radioPortOptions);
    radioSerialPort->start();

    statusTimer = std::make_unique<asio::steady_timer>(io);
    scheduleStatusTimer();

    lcdTimer = std::make_unique<asio::steady_timer>(io);
    scheduleLCDUpdate();

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

void Daemon::initLCD() {
    try {
        lcd.setBus(config.getLCDI2CBus());
        lcd.setAddress(config.getLCDI2CAddress());
        lcd.init();
        lcd.clear();
        lcd.setLine(0, "   Groundstation");
        lcd.setLine(1, "    Starting...");   
    } catch (const std::exception& e) {
        error("Failed to initialize LCD: {}", e.what());
    }
}

void Daemon::scheduleLCDUpdate() {
    lcdTimer->expires_after(std::chrono::seconds(5));
    lcdTimer->async_wait([this](const asio::error_code& ec) {
        if (!ec) {
            try {
                lcd.clear();
                auto gpsTime = gps.getUTCTime();
                if (gpsTime.has_value()) {
                    lcd.setCursor(0, 0);
                    lcd.print("{}", *gpsTime);
                } else {
                    lcd.setCursor(0, 0);
                    lcd.print("");
                }
                auto gpsPos = gps.getPosition();
                if (gpsPos.has_value()) {
                    lcd.setCursor(0, 1);
                    lcd.print("{:.4f} {:.4f}", 
                        gpsPos->latInRadians * RADIANS_TO_DEGREES, 
                        gpsPos->lonInRadians * RADIANS_TO_DEGREES);
                } else {
                    lcd.setCursor(0, 1);
                    lcd.print("No GPS Fix");
                }
                auto rotAz = rotator.getAzimuth();
                auto rotEl = rotator.getElevation();
                if (rotAz.has_value() && rotEl.has_value()) {
                    lcd.setCursor(0, 2);
                    lcd.print("Az: {:.1f} El: {:.1f}", rotAz.value(), rotEl.value());
                } else {
                    lcd.setCursor(0, 2);
                    lcd.print("Looking for Rotator");
                }
                lcd.setCursor(0, 3);
                lcd.print("     Waiting...");
            } catch (const std::exception& e) {
                error("Failed to update LCD: {}", e.what());
            }
            scheduleLCDUpdate(); // Reschedule the timer
        } else {
            error("LCD update timer error: {}", ec.message());
        }
    });
}

void Daemon::scheduleStatusTimer() {
    statusTimer->expires_after(std::chrono::seconds(config.getStatusIntervalSeconds()));
    statusTimer->async_wait([this](const asio::error_code& ec) {
        if (!ec) {
            info("--- Groundstation Status ---");
            auto gpsPos = gps.getPosition();
            if (gpsPos.has_value()) {
                info("- GPS Position: Lat {:.6f}°, Lon {:.6f}°, Alt {:.2f} km",
                     gpsPos->latInRadians * RADIANS_TO_DEGREES,
                     gpsPos->lonInRadians * RADIANS_TO_DEGREES,
                     gpsPos->altInKilometers);
            } else {
                info("- GPS Position: No fix");
            }
            auto gpsTime = gps.getUTCTime();
            if (gpsTime.has_value()) {
                info("- GPS Time: {}", gpsTime.value());
            } else {
                info("- GPS Time: Unknown");
            }
            auto rotAz = rotator.getAzimuth();
            auto rotEl = rotator.getElevation();
            if (rotAz.has_value() && rotEl.has_value()) {
                info("- Rotator Position: Az {:.2f}°, El {:.2f}°", rotAz.value(), rotEl.value());
            } else {
                info("- Rotator Position: Unknown");
            }
            info("----------------------------");
            scheduleStatusTimer(); // Reschedule the timer
        } else {
            error("Status logging timer error: {}", ec.message());
        }
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
    lcd.clear();
    lcd.setLine(0, "   Groundstation");
    lcd.setLine(1, "      Offline");   
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
