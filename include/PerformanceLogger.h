#pragma once

#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <memory>
#include <mutex>

class PerformanceLogger {
public:
    enum class LogLevel {
        DEBUG = 0,
        INFO = 1,
        WARNING = 2,
        ERROR = 3,
        CRITICAL = 4
    };

    static PerformanceLogger& getInstance();
    
    void initialize(const std::string& logFilePath = "blackhole_performance.log");
    void shutdown();
    
    // Core logging functions
    void log(LogLevel level, const std::string& message);
    void logFrame(double frameTime, double fps, const std::string& details = "");
    void logPerformanceDrop(double frameTime, const std::string& cause = "");
    void logGPUError(const std::string& error);
    void logRenderingStage(const std::string& stage, double duration);
    void logOpenGLInfo(); // Call after OpenGL context is created
    
    // Timing utilities
    class Timer {
    public:
        Timer(const std::string& name, PerformanceLogger& logger);
        ~Timer();
        
    private:
        std::string name_;
        PerformanceLogger& logger_;
        std::chrono::high_resolution_clock::time_point start_time_;
    };
    


private:
    PerformanceLogger() = default;
    ~PerformanceLogger() = default;
    PerformanceLogger(const PerformanceLogger&) = delete;
    PerformanceLogger& operator=(const PerformanceLogger&) = delete;
    
    void writeToFile(const std::string& message);
    std::string getCurrentTimestamp();
    std::string logLevelToString(LogLevel level);
    
    std::unique_ptr<std::ofstream> log_file_;
    std::mutex log_mutex_;
    bool initialized_ = false;
    
    // Performance tracking
    struct PerformanceData {
        double total_frame_time = 0.0;
        int frame_count = 0;
        double max_frame_time = 0.0;
        double min_frame_time = 999.0;
        int critical_drops = 0;
        int major_drops = 0;
        int emergency_drops = 0;
    };
    
    PerformanceData perf_data_;
    std::chrono::high_resolution_clock::time_point session_start_;
};

// Convenience macros
#define PERF_LOG_DEBUG(msg) PerformanceLogger::getInstance().log(PerformanceLogger::LogLevel::DEBUG, msg)
#define PERF_LOG_INFO(msg) PerformanceLogger::getInstance().log(PerformanceLogger::LogLevel::INFO, msg)
#define PERF_LOG_WARNING(msg) PerformanceLogger::getInstance().log(PerformanceLogger::LogLevel::WARNING, msg)
#define PERF_LOG_ERROR(msg) PerformanceLogger::getInstance().log(PerformanceLogger::LogLevel::ERROR, msg)
#define PERF_LOG_CRITICAL(msg) PerformanceLogger::getInstance().log(PerformanceLogger::LogLevel::CRITICAL, msg)
#define PERF_TIMER(name) PerformanceLogger::Timer _timer(name, PerformanceLogger::getInstance())
