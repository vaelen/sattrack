/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_DAEMON_HPP
#define __SATTRACK_DAEMON_HPP

#include <mutex>
#include <thread>
#include <queue>
#include <condition_variable>

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
    Daemon() = default;
    ~Daemon();

    DaemonStatus status();
    void start();
    void stop();
    void send(DaemonEvent event);

private:
    DaemonStatus _status = DaemonStatus::STOPPED;
    std::mutex _statusMutex;
    std::thread _eventLoopThread;
    std::mutex _eventMutex;
    std::queue<DaemonEvent> _eventQueue;
    std::condition_variable _eventCV;

    void eventLoop();
};

}

#endif
