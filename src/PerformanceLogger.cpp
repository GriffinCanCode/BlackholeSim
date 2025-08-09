#include "PerformanceLogger.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>

#ifdef __APPLE__
#include <OpenGL/gl3.h>
#else
#include <GL/gl.h>
#endif

PerformanceLogger& PerformanceLogger::getInstance() {
    static PerformanceLogger instance;
    return instance;
}

void PerformanceLogger::initialize(const std::string& logFilePath) {
    std::lock_guard<std::mutex> lock(log_mutex_);
    
    if (initialized_) {
        return;
    }
    
    log_file_ = std::make_unique<std::ofstream>(logFilePath, std::ios::out | std::ios::app);
    if (!log_file_->is_open()) {
        std::cerr << "Failed to open performance log file: " << logFilePath << std::endl;
        return;
    }
    
    initialized_ = true;
    session_start_ = std::chrono::high_resolution_clock::now();
    
    // Write session header
    writeToFile("=== NEW PERFORMANCE LOGGING SESSION STARTED ===");
    log(LogLevel::INFO, "Performance logger initialized - logging to " + logFilePath);
    
    // Note: OpenGL system info will be logged later after context is created
}

void PerformanceLogger::shutdown() {
    std::lock_guard<std::mutex> lock(log_mutex_);
    
    if (!initialized_) {
        return;
    }
    
    // Write session summary
    auto session_duration = std::chrono::high_resolution_clock::now() - session_start_;
    double session_seconds = std::chrono::duration<double>(session_duration).count();
    
    writeToFile("\n=== PERFORMANCE SESSION SUMMARY ===");
    writeToFile("Session Duration: " + std::to_string(session_seconds) + " seconds");
    writeToFile("Total Frames: " + std::to_string(perf_data_.frame_count));
    
    if (perf_data_.frame_count > 0) {
        double avg_frame_time = perf_data_.total_frame_time / perf_data_.frame_count;
        writeToFile("Average Frame Time: " + std::to_string(avg_frame_time * 1000.0) + " ms");
        writeToFile("Max Frame Time: " + std::to_string(perf_data_.max_frame_time * 1000.0) + " ms");
        writeToFile("Min Frame Time: " + std::to_string(perf_data_.min_frame_time * 1000.0) + " ms");
        writeToFile("Average FPS: " + std::to_string(1.0 / avg_frame_time));
    }
    
    writeToFile("Critical Drops (>1000ms): " + std::to_string(perf_data_.critical_drops));
    writeToFile("Major Drops (>200ms): " + std::to_string(perf_data_.major_drops));
    writeToFile("Emergency Drops (>100ms): " + std::to_string(perf_data_.emergency_drops));
    writeToFile("=== SESSION ENDED ===\n");
    
    if (log_file_) {
        log_file_->close();
        log_file_.reset();
    }
    
    initialized_ = false;
}

void PerformanceLogger::log(LogLevel level, const std::string& message) {
    if (!initialized_) return;
    
    std::string log_entry = "[" + getCurrentTimestamp() + "] [" + logLevelToString(level) + "] " + message;
    writeToFile(log_entry);
    
    // Also print critical/error messages to console
    if (level >= LogLevel::ERROR) {
        std::cout << log_entry << std::endl;
    }
}

void PerformanceLogger::logFrame(double frameTime, double fps, const std::string& details) {
    if (!initialized_) return;
    
    // Update statistics
    perf_data_.total_frame_time += frameTime;
    perf_data_.frame_count++;
    perf_data_.max_frame_time = std::max(perf_data_.max_frame_time, frameTime);
    perf_data_.min_frame_time = std::min(perf_data_.min_frame_time, frameTime);
    
    // Only log detailed frame info for problematic frames or periodically
    static int frame_counter = 0;
    bool should_log = (frameTime > 0.033) || (++frame_counter % 60 == 0); // Log if >33ms or every 60 frames
    
    if (should_log) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2);
        oss << "FRAME: " << frameTime * 1000.0 << "ms | " << fps << " FPS";
        if (!details.empty()) {
            oss << " | " << details;
        }
        
        LogLevel level = LogLevel::DEBUG;
        if (frameTime > 0.1) level = LogLevel::ERROR;
        else if (frameTime > 0.05) level = LogLevel::WARNING;
        else if (frameTime > 0.033) level = LogLevel::INFO;
        
        log(level, oss.str());
    }
}

void PerformanceLogger::logPerformanceDrop(double frameTime, const std::string& cause) {
    if (!initialized_) return;
    
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);
    oss << "PERFORMANCE DROP: " << frameTime * 1000.0 << "ms";
    if (!cause.empty()) {
        oss << " | Cause: " << cause;
    }
    
    LogLevel level = LogLevel::WARNING;
    if (frameTime > 1.0) {
        level = LogLevel::CRITICAL;
        perf_data_.critical_drops++;
    } else if (frameTime > 0.2) {
        level = LogLevel::ERROR;
        perf_data_.major_drops++;
    } else if (frameTime > 0.1) {
        level = LogLevel::WARNING;
        perf_data_.emergency_drops++;
    }
    
    log(level, oss.str());
    
    // Check for OpenGL errors when we have performance drops
    GLenum error = glGetError();
    if (error != GL_NO_ERROR) {
        logGPUError("OpenGL Error during performance drop: " + std::to_string(error));
    }
}

void PerformanceLogger::logGPUError(const std::string& error) {
    if (!initialized_) return;
    log(LogLevel::ERROR, "GPU ERROR: " + error);
}

void PerformanceLogger::logRenderingStage(const std::string& stage, double duration) {
    if (!initialized_) return;
    
    // Only log stages that take significant time
    if (duration > 0.001) { // > 1ms
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3);
        oss << "STAGE [" << stage << "]: " << duration * 1000.0 << "ms";
        
        LogLevel level = LogLevel::DEBUG;
        if (duration > 0.05) level = LogLevel::ERROR;
        else if (duration > 0.02) level = LogLevel::WARNING;
        else if (duration > 0.01) level = LogLevel::INFO;
        
        log(level, oss.str());
    }
}

void PerformanceLogger::logOpenGLInfo() {
    if (!initialized_) return;
    
    // Log system information - call this after OpenGL context is created
    const GLubyte* renderer = glGetString(GL_RENDERER);
    const GLubyte* vendor = glGetString(GL_VENDOR);
    const GLubyte* version = glGetString(GL_VERSION);
    
    if (renderer) log(LogLevel::INFO, "GPU Renderer: " + std::string(reinterpret_cast<const char*>(renderer)));
    if (vendor) log(LogLevel::INFO, "GPU Vendor: " + std::string(reinterpret_cast<const char*>(vendor)));
    if (version) log(LogLevel::INFO, "OpenGL Version: " + std::string(reinterpret_cast<const char*>(version)));
}

// Timer implementation
PerformanceLogger::Timer::Timer(const std::string& name, PerformanceLogger& logger)
    : name_(name), logger_(logger), start_time_(std::chrono::high_resolution_clock::now()) {
}

PerformanceLogger::Timer::~Timer() {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(end_time - start_time_).count();
    logger_.logRenderingStage(name_, duration);
}

// Private methods
void PerformanceLogger::writeToFile(const std::string& message) {
    if (log_file_ && log_file_->is_open()) {
        *log_file_ << message << std::endl;
        log_file_->flush(); // Ensure immediate write for crash scenarios
    }
}

std::string PerformanceLogger::getCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;
    
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&time_t), "%H:%M:%S");
    oss << '.' << std::setfill('0') << std::setw(3) << ms.count();
    
    return oss.str();
}

std::string PerformanceLogger::logLevelToString(LogLevel level) {
    switch (level) {
        case LogLevel::DEBUG: return "DEBUG";
        case LogLevel::INFO: return "INFO";
        case LogLevel::WARNING: return "WARN";
        case LogLevel::ERROR: return "ERROR";
        case LogLevel::CRITICAL: return "CRIT";
        default: return "UNKNOWN";
    }
}
