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
    if (_eventLoopThread.joinable()) {
        info("Waiting for event loop thread to finish...");
        _eventLoopThread.join();
    }
    info("Daemon destroyed.");
}

DaemonStatus Daemon::status() {
    std::scoped_lock lock(_statusMutex);
    return _status;
}

void Daemon::start() {
    std::scoped_lock lock(_statusMutex);
    if (_status == DaemonStatus::STOPPED) {
        info("Starting daemon...");
        _status = DaemonStatus::STARTING;
        _eventLoopThread = std::thread([this]{ eventLoop(); });
    }
}

void Daemon::stop() {
    std::scoped_lock lock(_statusMutex);
    if (_status == DaemonStatus::RUNNING) {
        info("Stopping daemon...");
        _status = DaemonStatus::STOPPING;
    }
}

void Daemon::send(DaemonEvent event) {
    {
        std::scoped_lock lock(_eventMutex);
        _eventQueue.push(event);
    }
    _eventCV.notify_one();
}

void Daemon::eventLoop() {
    // Update status to RUNNING
    {
        std::scoped_lock lock(_statusMutex);
        _status = DaemonStatus::RUNNING;
        info("Daemon started.");
    }

    while (status() == DaemonStatus::RUNNING) {
        std::unique_lock eventLock(_eventMutex);
        _eventCV.wait(eventLock, [this]{
            return !_eventQueue.empty() || status() != DaemonStatus::RUNNING;
        });
        while (!_eventQueue.empty()) {
            DaemonEvent event = _eventQueue.front();
            _eventQueue.pop();
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
    {
        std::scoped_lock lock(_statusMutex);
        _status = DaemonStatus::STOPPED;
        info("Daemon stopped.");
    }
    
}

}
