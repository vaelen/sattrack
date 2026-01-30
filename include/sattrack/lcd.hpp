/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#ifndef __SATTRACK_LCD_HPP
#define __SATTRACK_LCD_HPP

#include <mutex>
#include <string>
#include <cstdint>

#include <i2clcd.h>
#include <spdlog/fmt/fmt.h>

namespace sattrack {

/**
 * LCD class for controlling an I2C LCD display.
 *
 * Wraps the i2clcd C library with a thread-safe C++ interface.
 * Supports HD44780-compatible displays via PCF8574 I2C backpack.
 */
class LCD {
public:
    LCD() = default;
    ~LCD();

    // Non-copyable, non-movable (due to mutex and C handle)
    LCD(const LCD&) = delete;
    LCD& operator=(const LCD&) = delete;
    LCD(LCD&&) = delete;
    LCD& operator=(LCD&&) = delete;

    /**
     * Set the I2C bus number.
     * @param bus I2C bus number (e.g., 1 for /dev/i2c-1)
     */
    void setBus(int bus);

    /**
     * Set the I2C device address.
     * @param address I2C address (typically 0x27 or 0x3F)
     */
    void setAddress(int address);

    /**
     * Set the display size.
     * @param cols Number of columns (typically 16 or 20)
     * @param rows Number of rows (typically 2 or 4)
     */
    void setSize(uint8_t cols, uint8_t rows);

    /**
     * Initialize the LCD display.
     * Must be called after setting bus, address, and size.
     * @throws std::runtime_error if initialization fails
     */
    void init();

    /**
     * Deinitialize the LCD and release resources.
     */
    void deinit();

    /**
     * Clear the entire display.
     */
    void clear();

    /**
     * Clear a specific line.
     * @param line Line number (0-indexed)
     */
    void clearLine(uint8_t line);

    /**
     * Move cursor to home position (0,0).
     */
    void home();

    /**
     * Turn display on or off.
     * @param on true to turn on, false to turn off
     */
    void display(bool on);

    /**
     * Turn backlight on or off.
     * @param on true to turn on, false to turn off
     */
    void backlight(bool on);

    /**
     * Set cursor position.
     * @param col Column (0-indexed)
     * @param row Row (0-indexed)
     */
    void setCursor(uint8_t col, uint8_t row);

    /**
     * Show or hide cursor.
     * @param visible true to show cursor, false to hide
     */
    void cursor(bool visible);

    /**
     * Enable or disable cursor blinking.
     * @param on true to enable blinking, false to disable
     */
    void blink(bool on);

    /**
     * Write a single character to the display.
     * @param c Character to write
     */
    void putChar(char c);

    /**
     * Write a string to the display.
     * @param str String to write
     */
    void print(const std::string& str);

    /**
     * Write formatted text to the display.
     * @param fmt Format string
     * @param args Format arguments
     */
    template<typename... Args>
    void print(fmt::format_string<Args...> format, Args&&... args) {
        print(fmt::format(format, std::forward<Args>(args)...));
    }

    /**
     * Set the content of an entire line.
     * @param line Line number (0-indexed)
     * @param text Text to display on the line
     */
    void setLine(uint8_t line, const std::string& text);

    /**
     * Create a custom character.
     * @param location Character slot (0-7)
     * @param charmap 8-byte character bitmap
     */
    void createChar(uint8_t location, const uint8_t charmap[8]);

    /**
     * Check if the LCD is initialized.
     * @return true if initialized, false otherwise
     */
    bool isInitialized() const;

private:
    mutable std::mutex mutex_;
    std::string devicePath_ = "/dev/i2c-1";
    i2clcd_config_t config_ = I2CLCD_CONFIG_DEFAULT;
    i2clcd_t* handle_ = nullptr;
    bool initialized_ = false;

    /**
     * Check an i2clcd error code and log warnings if necessary.
     * @param err Error code from i2clcd function
     * @param operation Description of the operation for logging
     */
    void checkError(i2clcd_err_t err, const std::string& operation);
};

} // namespace sattrack

#endif
