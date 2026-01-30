/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack/lcd.hpp>
#include <spdlog/spdlog.h>
#include <stdexcept>

using spdlog::debug;
using spdlog::info;
using spdlog::warn;
using spdlog::error;

namespace sattrack {

LCD::~LCD() {
    deinit();
}

void LCD::setBus(int bus) {
    std::lock_guard<std::mutex> lock(mutex_);
    devicePath_ = fmt::format("/dev/i2c-{}", bus);
    config_.i2c_device = devicePath_.c_str();
}

void LCD::setAddress(int address) {
    std::lock_guard<std::mutex> lock(mutex_);
    config_.i2c_addr = address;
}

void LCD::setSize(uint8_t cols, uint8_t rows) {
    std::lock_guard<std::mutex> lock(mutex_);
    config_.cols = cols;
    config_.rows = rows;
}

void LCD::init() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (initialized_) {
        return;
    }

    // Ensure device path is set in config
    config_.i2c_device = devicePath_.c_str();

    debug("Initializing LCD on {} at address 0x{:02X} ({}x{})",
          devicePath_, config_.i2c_addr, config_.cols, config_.rows);

    i2clcd_err_t err = i2clcd_init(&config_, &handle_);
    if (err != I2CLCD_OK) {
        throw std::runtime_error(
            fmt::format("Failed to initialize LCD: {}", i2clcd_strerror(err)));
    }

    initialized_ = true;
    info("LCD initialized successfully");
}

void LCD::deinit() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }

    debug("Deinitializing LCD");
    i2clcd_deinit(handle_);
    handle_ = nullptr;
    initialized_ = false;
}

void LCD::clear() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }
    checkError(i2clcd_clear(handle_), "clear");
}

void LCD::clearLine(uint8_t line) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }
    checkError(i2clcd_clear_line(handle_, line), "clear line");
}

void LCD::home() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }
    checkError(i2clcd_home(handle_), "home");
}

void LCD::display(bool on) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }
    checkError(i2clcd_display(handle_, on), "display");
}

void LCD::backlight(bool on) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }
    checkError(i2clcd_backlight(handle_, on), "backlight");
}

void LCD::setCursor(uint8_t col, uint8_t row) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }
    checkError(i2clcd_set_cursor(handle_, col, row), "set cursor");
}

void LCD::cursor(bool visible) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }
    checkError(i2clcd_cursor(handle_, visible), "cursor");
}

void LCD::blink(bool on) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }
    checkError(i2clcd_blink(handle_, on), "blink");
}

void LCD::putChar(char c) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }
    checkError(i2clcd_putc(handle_, c), "putc");
}

void LCD::print(const std::string& str) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }
    checkError(i2clcd_puts(handle_, str.c_str()), "print");
}

void LCD::setLine(uint8_t line, const std::string& text) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }
    checkError(i2clcd_set_line(handle_, line, text.c_str()), "set line");
}

void LCD::createChar(uint8_t location, const uint8_t charmap[8]) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!initialized_ || handle_ == nullptr) {
        return;
    }
    checkError(i2clcd_create_char(handle_, location, charmap), "create char");
}

bool LCD::isInitialized() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return initialized_;
}

void LCD::checkError(i2clcd_err_t err, const std::string& operation) {
    if (err != I2CLCD_OK) {
        warn("LCD {} failed: {}", operation, i2clcd_strerror(err));
    }
}

} // namespace sattrack
