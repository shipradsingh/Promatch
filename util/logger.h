#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <mutex>
#include <sstream>
#include <chrono>
#include <iomanip>

namespace qrc {

// Log levels in order of increasing severity
enum class LogLevel {
    DEBUG = 0,
    INFO = 1,
    WARN = 2,
    ERROR = 3,
    NONE = 4  // Use to disable all logging
};

class Logger {
public:
    // Get the singleton instance
    static Logger& getInstance() {
        static Logger instance;
        return instance;
    }

    // Set the minimum log level (all logs below this level will be ignored)
    void setLevel(LogLevel level) {
        m_level = level;
    }

    // Get the current log level
    LogLevel getLevel() const {
        return m_level;
    }

    // Set output to file
    void setLogFile(const std::string& filename) {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_fileStream.is_open()) {
            m_fileStream.close();
        }
        m_fileStream.open(filename, std::ios::out | std::ios::app);
        m_useFileOutput = m_fileStream.is_open();
    }

    // Close the log file and revert to console output
    void closeLogFile() {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_fileStream.is_open()) {
            m_fileStream.close();
        }
        m_useFileOutput = false;
    }

    // Log methods for each level
    template<typename... Args>
    void debug(Args&&... args) {
        log(LogLevel::DEBUG, std::forward<Args>(args)...);
    }

    template<typename... Args>
    void info(Args&&... args) {
        log(LogLevel::INFO, std::forward<Args>(args)...);
    }

    template<typename... Args>
    void warn(Args&&... args) {
        log(LogLevel::WARN, std::forward<Args>(args)...);
    }

    template<typename... Args>
    void error(Args&&... args) {
        log(LogLevel::ERROR, std::forward<Args>(args)...);
    }

    // Parse a string to get the log level
    static LogLevel parseLogLevel(const std::string& levelStr) {
        if (levelStr == "DEBUG") return LogLevel::DEBUG;
        if (levelStr == "INFO") return LogLevel::INFO;
        if (levelStr == "WARN") return LogLevel::WARN;
        if (levelStr == "ERROR") return LogLevel::ERROR;
        if (levelStr == "NONE") return LogLevel::NONE;
        return LogLevel::INFO; // Default to INFO if invalid
    }

private:
    // Private constructor to enforce singleton pattern
    Logger() : m_level(LogLevel::INFO), m_useFileOutput(false) {}
    
    // No copying allowed
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    LogLevel m_level;
    std::ofstream m_fileStream;
    std::mutex m_mutex;
    bool m_useFileOutput;

    std::string getLevelString(LogLevel level) {
        switch (level) {
            case LogLevel::DEBUG: return "DEBUG";
            case LogLevel::INFO:  return "INFO ";
            case LogLevel::WARN:  return "WARN ";
            case LogLevel::ERROR: return "ERROR";
            default:              return "UNKNOWN";
        }
    }

    std::string getCurrentTimestamp() {
        auto now = std::chrono::system_clock::now();
        auto now_time_t = std::chrono::system_clock::to_time_t(now);
        auto now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()) % 1000;

        std::stringstream ss;
        ss << std::put_time(std::localtime(&now_time_t), "%Y-%m-%d %H:%M:%S");
        ss << '.' << std::setfill('0') << std::setw(3) << now_ms.count();
        return ss.str();
    }

    // Helper to print individual arguments
    template<typename T>
    void printArg(std::stringstream& ss, T&& arg) {
        ss << arg;
    }

    template<typename T, typename... Args>
    void printArg(std::stringstream& ss, T&& arg, Args&&... args) {
        ss << arg;
        printArg(ss, std::forward<Args>(args)...);
    }

    // Main log function
    template<typename... Args>
    void log(LogLevel level, Args&&... args) {
        if (level < m_level) {
            return; // Skip logging if below the current level
        }

        std::stringstream ss;
        ss << "[" << getCurrentTimestamp() << "] [" << getLevelString(level) << "] ";
        printArg(ss, std::forward<Args>(args)...);

        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_useFileOutput && m_fileStream.is_open()) {
            m_fileStream << ss.str() << std::endl;
        } else {
            std::cout << ss.str() << std::endl;
        }
    }
};

// Create a global logger instance for easier access
inline Logger& log = Logger::getInstance();

// Function to parse command line arguments to set log level
inline void initLogging(int argc, char* argv[]) {
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--log-level" && i + 1 < argc) {
            log.setLevel(Logger::parseLogLevel(argv[i + 1]));
            return;
        }
    }
    // Default to INFO if no flag is provided
    log.setLevel(LogLevel::INFO);
}

} // namespace qrc