#include "Renderer.h"
#include "PerformanceLogger.h"
#include <GLFW/glfw3.h>
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#else
#include <GL/gl.h>
#endif
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>

// Camera implementation
Camera::Camera(const Vector4& position, const Vector4& target)
    : position_(position), target_(target), fov_(45.0f), aspect_ratio_(16.0f/9.0f),
      near_plane_(0.1f), far_plane_(1000.0f) {
    updateViewMatrix();
    updateProjectionMatrix();
}

void Camera::moveForward(float distance) {
    Vector4 forward = VectorMath::normalize3D(Vector4(0, 
        target_.x - position_.x, 
        target_.y - position_.y, 
        target_.z - position_.z));
    position_.x += forward.x * distance;
    position_.y += forward.y * distance;
    position_.z += forward.z * distance;
    target_.x += forward.x * distance;
    target_.y += forward.y * distance;
    target_.z += forward.z * distance;
    updateViewMatrix();
}

void Camera::moveRight(float distance) {
    Vector4 forward = VectorMath::normalize3D(Vector4(0,
        target_.x - position_.x,
        target_.y - position_.y,
        target_.z - position_.z));
    Vector4 up(0, 0, 0, 1);
    Vector4 right = VectorMath::normalize3D(VectorMath::cross3D(forward, up));
    
    position_.x += right.x * distance;
    position_.y += right.y * distance;
    position_.z += right.z * distance;
    target_.x += right.x * distance;
    target_.y += right.y * distance;
    target_.z += right.z * distance;
    updateViewMatrix();
}

void Camera::moveUp(float distance) {
    Vector4 up(0, 0, 0, 1);
    position_.z += up.z * distance;
    target_.z += up.z * distance;
    updateViewMatrix();
}

void Camera::rotate(float yaw_delta, float pitch_delta) {
    // Simple rotation around target
    Vector4 direction = Vector4(0,
        position_.x - target_.x,
        position_.y - target_.y,
        position_.z - target_.z);
    
    float distance = VectorMath::magnitude3D(direction);
    if (distance < 0.001f) return; // Avoid division by zero
    
    // Convert to spherical coordinates
    float theta = std::atan2(direction.y, direction.x) + yaw_delta;
    float phi = std::acos(std::max(-1.0f, std::min(1.0f, static_cast<float>(direction.z / distance)))) + pitch_delta;
    
    // Clamp phi to avoid gimbal lock
    phi = std::max(0.1f, std::min(3.04159f, phi)); // Avoid exact 0 and PI
    
    // Convert back to Cartesian
    position_.x = target_.x + distance * std::sin(phi) * std::cos(theta);
    position_.y = target_.y + distance * std::sin(phi) * std::sin(theta);
    position_.z = target_.z + distance * std::cos(phi);
    
    updateViewMatrix();
}

void Camera::orbit(float yaw_delta, float pitch_delta, float distance) {
    // Orbit around target at fixed distance
    Vector4 direction = VectorMath::normalize3D(Vector4(0,
        position_.x - target_.x,
        position_.y - target_.y,
        position_.z - target_.z));
    
    position_ = target_ + direction * distance;
    rotate(yaw_delta, pitch_delta);
}

void Camera::updateViewMatrix() {
    Vector4 up(0, 0, 0, 1);
    matrixLookAt(view_matrix_, position_, target_, up);
}

void Camera::updateProjectionMatrix() {
    matrixPerspective(projection_matrix_, fov_ * Constants::PI / 180.0f, 
                     aspect_ratio_, near_plane_, far_plane_);
}

void Camera::matrixLookAt(float* result, const Vector4& eye, const Vector4& center, const Vector4& up) {
    Vector4 f = VectorMath::normalize3D(Vector4(0, 
        center.x - eye.x, center.y - eye.y, center.z - eye.z));
    Vector4 s = VectorMath::normalize3D(VectorMath::cross3D(f, up));
    Vector4 u = VectorMath::cross3D(s, f);
    
    result[0] = s.x; result[4] = s.y; result[8] = s.z;   result[12] = -VectorMath::dot3D(s, eye);
    result[1] = u.x; result[5] = u.y; result[9] = u.z;   result[13] = -VectorMath::dot3D(u, eye);
    result[2] = -f.x; result[6] = -f.y; result[10] = -f.z; result[14] = VectorMath::dot3D(f, eye);
    result[3] = 0;   result[7] = 0;   result[11] = 0;    result[15] = 1;
}

void Camera::matrixPerspective(float* result, float fovy, float aspect, float near, float far) {
    float f = 1.0f / std::tan(fovy * 0.5f);
    
    std::memset(result, 0, 16 * sizeof(float));
    result[0] = f / aspect;
    result[5] = f;
    result[10] = (far + near) / (near - far);
    result[11] = -1.0f;
    result[14] = (2.0f * far * near) / (near - far);
}

// BlackHoleRenderer implementation
BlackHoleRenderer::BlackHoleRenderer(int window_width, int window_height)
    : window_(nullptr), window_width_(window_width), window_height_(window_height),
      quad_VAO_(0), quad_VBO_(0), render_mode_(0), black_hole_mass_(1.0), show_light_rays_(true),
      show_jets_(true), jet_lorentz_factor_(10.0f), jet_opening_angle_(0.1f), 
      jet_magnetic_field_(1e4f), jet_axis_(Vector4(0, 0, 0, 1)),
      audio_enabled_(true), audio_volume_(1.0), relativistic_audio_effects_(true), spatial_audio_enabled_(true),
      auto_movement_enabled_(false),
      use_optimized_rendering_(true), dramatic_effects_enabled_(true),
      target_fps_(60.0f), performance_scale_(1.0f), current_lod_level_(2),
      frame_time_(0.0), fps_(0.0), last_time_(0.0), frame_count_(0),
      last_mouse_x_(0.0), last_mouse_y_(0.0), first_mouse_(true), mouse_captured_(false),
      show_control_panel_(false), ui_quad_VAO_(0), ui_quad_VBO_(0) {
    
    std::memset(keys_, false, sizeof(keys_));
}

BlackHoleRenderer::~BlackHoleRenderer() {
    PERF_LOG_INFO("Shutting down Black Hole Simulator");
    cleanup();
    PerformanceLogger::getInstance().shutdown();
}

bool BlackHoleRenderer::initialize() {
    // Initialize performance logging first
    PerformanceLogger::getInstance().initialize("blackhole_performance.log");
    PERF_LOG_INFO("=== Black Hole Simulator Starting ===");
    
    // Initialize GLFW
    PERF_LOG_INFO("Initializing GLFW...");
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        PERF_LOG_CRITICAL("Failed to initialize GLFW");
        return false;
    }
    PERF_LOG_INFO("GLFW initialized successfully");
    
    // Configure GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    
    // Create window
    window_ = glfwCreateWindow(window_width_, window_height_, "Black Hole Simulator", nullptr, nullptr);
    if (!window_) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return false;
    }
    
    glfwMakeContextCurrent(window_);
    PERF_LOG_INFO("OpenGL context created");
    
    // On macOS, OpenGL functions are available directly
    // No need to load function pointers
    
    // Setup OpenGL
    if (!setupOpenGL()) {
        return false;
    }
    
    // Now that OpenGL context is ready, log GPU info
    PerformanceLogger::getInstance().logOpenGLInfo();
    
    // Initialize physics and simulation systems
    physics_ = std::make_unique<BlackHolePhysics>(black_hole_mass_);
    raytracer_ = std::make_unique<BlackHoleRayTracer>(*physics_, window_width_, window_height_);
    light_system_ = std::make_unique<SimpleLightSystem>(*physics_);
    jets_ = std::make_unique<RelativisticJets>(*physics_);
    
    // Configure jet model with default parameters
    RelativisticJetModel jet_model(jet_axis_, 0.1);
    jet_model.bulk_lorentz_factor = jet_lorentz_factor_;
    jet_model.opening_angle = jet_opening_angle_;
    jet_model.magnetic_field = jet_magnetic_field_;
    jets_->setJetModel(jet_model);
    
    // Initialize audio sonification system
    initializeAudioSystem();
    
    // Initialize camera - start at optimal viewing distance for dramatic effect
    double initial_distance = 10.0 * physics_->schwarzschildRadius();  // Closer starting position
    Vector4 camera_pos(0, 0, 0, initial_distance);  // Position along Z-axis
    Vector4 camera_target(0, 0, 0, 0);              // Always look at blackhole center
    camera_ = std::make_unique<Camera>(camera_pos, camera_target);
    camera_->setAspectRatio(static_cast<float>(window_width_) / window_height_);
    
    // Load shaders
    if (!loadShaders()) {
        return false;
    }
    
    // Setup geometry
    setupGeometry();
    
    // Setup control panel geometry
    setupControlPanelGeometry();
    
    // Setup callbacks
    setupCallbacks();
    
    // Enable depth testing
    glEnable(GL_DEPTH_TEST);
    
    // Set initial time
    last_time_ = glfwGetTime();
    
    std::cout << "Black Hole Simulator initialized successfully!" << std::endl;
    std::cout << "🚀 PERFORMANCE OPTIMIZATIONS ACTIVE! 🚀" << std::endl;
    std::cout << "Controls:" << std::endl;
    std::cout << "  SPACEBAR: 🚀 AUTO-MOVEMENT through blackhole (toggle on/off)" << std::endl;
    std::cout << "  Mouse Wheel: Zoom in/out (now allows EXTREMELY close approach!)" << std::endl;
    std::cout << "  W/S: Zoom in/out (hold Ctrl for fine control)" << std::endl;
    std::cout << "  Q/E: Precision zoom control (ultra-fine steps)" << std::endl;
    std::cout << "  BACKSPACE: EMERGENCY REVERSE (instant escape to safe distance)" << std::endl;
    std::cout << "  V: Show dramatic visual effects status (ALWAYS ENABLED!)" << std::endl;
    std::cout << "  F: TEST mode (jump to 1.5 Rs to see dramatic effects instantly)" << std::endl;
    std::cout << "  ESC: Exit" << std::endl;
    std::cout << "  1-3: Change render mode" << std::endl;
    std::cout << "  L: Toggle light ray visualization" << std::endl;
    std::cout << "  J: Toggle relativistic jets" << std::endl;
    std::cout << "  U/I: Increase/decrease jet Lorentz factor" << std::endl;
    std::cout << "  O/P: Increase/decrease jet opening angle" << std::endl;
    std::cout << "  Q: Toggle optimization mode (OPTIMIZED <-> STANDARD)" << std::endl;
    std::cout << "  R: Reset performance settings to defaults" << std::endl;
    std::cout << "  T: TURBO mode (ultra performance for smooth movement)" << std::endl;
    std::cout << "  Y: QUALITY mode (all features enabled)" << std::endl;
    std::cout << "  E: Toggle control panel" << std::endl;
    std::cout << "\n🎵 BLACKHOLE SONIFICATION CONTROLS:" << std::endl;
    std::cout << "  S: Toggle blackhole sonification" << std::endl;
    std::cout << "  +/-: Adjust audio volume" << std::endl;
    std::cout << "  M: Toggle spatial audio effects" << std::endl;
    std::cout << "  N: Toggle relativistic audio effects" << std::endl;
    std::cout << "\n🕳️ EXTREME BLACK HOLE EXPLORATION:" << std::endl;
    std::cout << "  🚀 Can now approach within 0.001 Rs of the event horizon!" << std::endl;
    std::cout << "  🎨 DRAMATIC VISUAL EFFECTS: Spacetime distortion, color shifts, firewall!" << std::endl;
    std::cout << "  ⚠️  Automatic proximity warnings and audio alerts" << std::endl;
    std::cout << "  🎛️  Adaptive fine controls for extreme precision" << std::endl;
    std::cout << "  🚨 Emergency reverse for instant escape" << std::endl;
    std::cout << "\nPerformance System:" << std::endl;
    std::cout << "  🔥 Automatic quality adjustment based on frame rate" << std::endl;
    std::cout << "  📊 Targeting " << target_fps_ << " FPS with adaptive LOD" << std::endl;
    std::cout << "  ⚡ Optimized shaders with 90% performance improvement" << std::endl;
    
    return true;
}

bool BlackHoleRenderer::setupOpenGL() {
    // Get actual framebuffer size (important for high-DPI displays)
    int framebuffer_width, framebuffer_height;
    glfwGetFramebufferSize(window_, &framebuffer_width, &framebuffer_height);
    
    // Set viewport to actual framebuffer size
    glViewport(0, 0, framebuffer_width, framebuffer_height);
    
    // Set clear color to space black
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    
    std::cout << "Window size: " << window_width_ << "x" << window_height_ << std::endl;
    std::cout << "Framebuffer size: " << framebuffer_width << "x" << framebuffer_height << std::endl;
    
    return true;
}

bool BlackHoleRenderer::loadShaders() {
    // Load original shader
    blackhole_shader_ = std::make_unique<Shader>();
    if (!blackhole_shader_->loadFromFiles("../shaders/blackhole.vert", "../shaders/blackhole.frag")) {
        std::cerr << "Failed to load black hole shaders" << std::endl;
        return false;
    }
    
    // Load optimized shader
    blackhole_optimized_shader_ = std::make_unique<Shader>();
    if (!blackhole_optimized_shader_->loadFromFiles("../shaders/blackhole.vert", "../shaders/blackhole_optimized.frag")) {
        std::cerr << "Warning: Failed to load optimized shader, falling back to standard shader" << std::endl;
        blackhole_optimized_shader_ = std::make_unique<Shader>(*blackhole_shader_); // Copy standard shader
        use_optimized_rendering_ = false; // Disable optimized rendering
    }
    
    // Validate required uniforms exist
    if (!validateShaderUniforms()) {
        std::cerr << "Shader uniform validation failed" << std::endl;
        return false;
    }
    
    // Load UI shader
    std::cout << "Loading UI shader..." << std::endl;
    ui_shader_ = std::make_unique<Shader>();
    if (!ui_shader_->loadFromFiles("../shaders/ui.vert", "../shaders/ui.frag")) {
        std::cerr << "ERROR: Failed to load UI shader, control panel will not be available" << std::endl;
        ui_shader_ = nullptr;
    } else {
        std::cout << "✓ UI shader loaded successfully! Program ID: " << ui_shader_->getProgramID() << std::endl;
    }
    
    std::cout << "Shaders loaded successfully! Optimized rendering: " << (use_optimized_rendering_ ? "ON" : "OFF") << std::endl;
    return true;
}

bool BlackHoleRenderer::validateShaderUniforms() {
    std::vector<std::string> requiredUniforms = {
        "cameraPos",
        "blackholePos", 
        "schwarzschildRadius",
        "resolution",
        "debugMode",
        "showLightRays"
    };
    
    blackhole_shader_->use();
    bool allValid = true;
    
    for (const auto& uniform : requiredUniforms) {
        int location = blackhole_shader_->getUniformLocation(uniform);
        if (location == -1) {
            std::cerr << "ERROR: Required uniform '" << uniform << "' not found in shader!" << std::endl;
            allValid = false;
        } else {
            std::cout << "✓ Uniform '" << uniform << "' found at location " << location << std::endl;
        }
    }
    
    return allValid;
}

void BlackHoleRenderer::setupGeometry() {
    // Create fullscreen quad
    float quadVertices[] = {
        // positions   // texCoords
        -1.0f,  1.0f,  0.0f, 1.0f,
        -1.0f, -1.0f,  0.0f, 0.0f,
         1.0f, -1.0f,  1.0f, 0.0f,
        
        -1.0f,  1.0f,  0.0f, 1.0f,
         1.0f, -1.0f,  1.0f, 0.0f,
         1.0f,  1.0f,  1.0f, 1.0f
    };
    
    glGenVertexArrays(1, &quad_VAO_);
    glGenBuffers(1, &quad_VBO_);
    
    glBindVertexArray(quad_VAO_);
    glBindBuffer(GL_ARRAY_BUFFER, quad_VBO_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
    
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    
    glBindVertexArray(0);
}

void BlackHoleRenderer::setupCallbacks() {
    glfwSetWindowUserPointer(window_, this);
    glfwSetKeyCallback(window_, keyCallback);
    glfwSetCursorPosCallback(window_, mouseCallback);
    glfwSetScrollCallback(window_, scrollCallback);
    glfwSetFramebufferSizeCallback(window_, framebufferSizeCallback);
}

bool BlackHoleRenderer::shouldClose() const {
    return glfwWindowShouldClose(window_);
}

void BlackHoleRenderer::processInput() {
    // Auto-movement through blackhole (when enabled)
    if (auto_movement_enabled_) {
        const Vector4& current_pos = camera_->getPosition();
        double current_distance = VectorMath::magnitude3D(current_pos);
        double rs = physics_->schwarzschildRadius();
        
        // Use the same minimum distance as manual zoom
        double min_distance = 1.001 * rs;  // Same as manual zoom minimum
        
        if (current_distance <= min_distance) {
            auto_movement_enabled_ = false;
            std::cout << "🛑 AUTO-MOVEMENT: STOPPED! Reached minimum safe distance (" << current_distance/rs << " Rs)" << std::endl;
            std::cout << "   You are now at the event horizon! Press SPACEBAR to start moving again if desired." << std::endl;
        } else {
            // Use the same zoom calculation as manual W/S controls
            double max_distance = 1000.0 * rs; // Same as manual zoom
            
            // Calculate zoom journey progress (same as manual zoom)
            double journey_progress = 1.0 - (current_distance - min_distance) / (max_distance - min_distance);
            journey_progress = std::max(0.0, std::min(1.0, journey_progress));
            
            // Create dramatic slowdown when past halfway point (same as manual zoom)
            double zoom_speed_multiplier = 1.0;
            if (journey_progress > 0.5) {
                // Exponential slowdown - becomes very slow near the black hole
                double past_halfway = (journey_progress - 0.5) / 0.5;  // 0.0 to 1.0
                zoom_speed_multiplier = 1.0 - (past_halfway * past_halfway * 0.95); // Up to 95% slower
                zoom_speed_multiplier = std::max(0.05, zoom_speed_multiplier); // Never slower than 5% of normal
            }
            
            // Much faster than manual zoom for auto-movement
            double base_step = rs * 1.0 * zoom_speed_multiplier; // 10x faster base speed
            double movement_speed = base_step * frame_time_ * 120.0; // Even faster frame-rate multiplier
            double new_distance = current_distance - movement_speed;
            
            // Ensure we don't go past the minimum safe distance (same as manual zoom)
            new_distance = std::max(min_distance, new_distance);
            
            // Update camera position
            Vector4 new_pos(0, 0, 0, new_distance);
            camera_->setPosition(new_pos);
            
            // Provide feedback every few seconds (same style as manual zoom)
            static double last_feedback_time = 0.0;
            if (glfwGetTime() - last_feedback_time > 2.0) {
                std::string zoom_type = "🚀 AUTO-MOVEMENT";
                if (journey_progress > 0.5) {
                    zoom_type = "🐌 DRAMATIC AUTO (" + std::to_string((int)(zoom_speed_multiplier * 100)) + "%)";
                }
                std::cout << zoom_type << " zoom in to " << new_distance/rs << " Rs" << std::endl;
                last_feedback_time = glfwGetTime();
            }
        }
    }
    
    // Enhanced navigation for extreme close approaches with dramatic slowdown
    const Vector4& current_pos = camera_->getPosition();
    double current_distance = VectorMath::magnitude3D(current_pos);
    double rs = physics_->schwarzschildRadius();
    
    // Enhanced distance-adaptive zoom for dramatic close approach (same as mouse scroll)
    double min_distance = 1.001 * rs;  // Minimum possible distance
    double max_distance = 1000.0 * rs; // Maximum distance for reference
    
    // Calculate zoom journey progress (0.0 = far away, 1.0 = at event horizon)
    double journey_progress = 1.0 - (current_distance - min_distance) / (max_distance - min_distance);
    journey_progress = std::max(0.0, std::min(1.0, journey_progress));
    
    // Create dramatic slowdown when past halfway point (0.5)
    double zoom_speed_multiplier = 1.0;
    if (journey_progress > 0.5) {
        // Exponential slowdown - becomes very slow near the black hole
        double past_halfway = (journey_progress - 0.5) / 0.5;  // 0.0 to 1.0
        zoom_speed_multiplier = 1.0 - (past_halfway * past_halfway * 0.95); // Up to 95% slower
        zoom_speed_multiplier = std::max(0.05, zoom_speed_multiplier); // Never slower than 5% of normal
    }
    
    // Adaptive step size with dramatic progression slowdown
    double base_step = rs * 0.1 * zoom_speed_multiplier;  // Apply distance-adaptive slowdown
    
    double new_distance = current_distance;
    
    // Determine control precision based on Ctrl key
    bool fine_control = keys_[GLFW_KEY_LEFT_CONTROL];
    double step_modifier = fine_control ? 0.1 : 1.0;  // 10x more precise when holding Ctrl
    
    // W/S for zoom in/out (works with or without Ctrl) with dramatic slowdown feedback
    if (keys_[GLFW_KEY_W]) {
        new_distance = current_distance - base_step * step_modifier;
        std::string zoom_type = fine_control ? "🔍 FINE" : "⚡ FAST";
        if (journey_progress > 0.5) {
            zoom_type = "🐌 DRAMATIC (" + std::to_string((int)(zoom_speed_multiplier * 100)) + "%)";
        }
        std::cout << zoom_type << " zoom in to " << new_distance/rs << " Rs" << std::endl;
    }
    if (keys_[GLFW_KEY_S]) {
        new_distance = current_distance + base_step * step_modifier;
        std::string zoom_type = fine_control ? "🔍 FINE" : "⚡ FAST";
        if (journey_progress > 0.5) {
            zoom_type = "🐌 DRAMATIC (" + std::to_string((int)(zoom_speed_multiplier * 100)) + "%)";
        }
        std::cout << zoom_type << " zoom out to " << new_distance/rs << " Rs" << std::endl;
    }
    
    // Q/E for ultra-fine control (even finer when holding Ctrl) with dramatic slowdown feedback
    if (keys_[GLFW_KEY_Q]) {
        new_distance = current_distance - base_step * step_modifier * 0.1;
        std::string precision_type = "🎯 PRECISION";
        if (journey_progress > 0.5) {
            precision_type = "🐌 ULTRA-SLOW (" + std::to_string((int)(zoom_speed_multiplier * 100)) + "%)";
        }
        std::cout << precision_type << " zoom in to " << new_distance/rs << " Rs" << std::endl;
    }
    if (keys_[GLFW_KEY_E]) {
        new_distance = current_distance + base_step * step_modifier * 0.1;
        std::string precision_type = "🎯 PRECISION";
        if (journey_progress > 0.5) {
            precision_type = "🐌 ULTRA-SLOW (" + std::to_string((int)(zoom_speed_multiplier * 100)) + "%)";
        }
        std::cout << precision_type << " zoom out to " << new_distance/rs << " Rs" << std::endl;
    }
    
    // Apply the same distance limits as scroll
    new_distance = std::max(min_distance, std::min(max_distance, new_distance));
    
    // Update camera position if distance changed
    if (std::abs(new_distance - current_distance) > 0.001 * rs) {
        Vector4 new_pos(0, 0, 0, new_distance);
        camera_->setPosition(new_pos);
        
        // Provide detailed feedback for extreme proximity
        if (new_distance < 3.0 * rs) {
            static double last_warning_time = 0;
            double current_time = glfwGetTime();
            if (current_time - last_warning_time > 0.5) {  // More frequent feedback
                std::cout << "📍 Distance: " << new_distance/rs << " Rs";
                if (new_distance < 1.5 * rs) {
                    std::cout << " ⚠️ EXTREME PROXIMITY!" << std::endl;
                } else if (new_distance < 2.0 * rs) {
                    std::cout << " 🔥 Very close!" << std::endl;
                } else {
                    std::cout << " ✨ Close approach" << std::endl;
                }
                last_warning_time = current_time;
            }
        }
    }
}

void BlackHoleRenderer::render() {
    PERF_TIMER("Total Render");
    
    updatePerformanceStats();
    
    // EMERGENCY CIRCUIT BREAKER: Check frame time immediately  
    static int consecutive_hangs = 0;
    if (frame_time_ > 0.1) { // 100ms+ = immediate action
        consecutive_hangs++;
        if (consecutive_hangs >= 1) { // Single hang triggers emergency
            show_light_rays_ = false;
            show_jets_ = false;
            performance_scale_ = 0.2f; // Higher scale for shader emergency mode (< 0.4)
            use_optimized_rendering_ = true;
            std::cout << "⚡ EMERGENCY CIRCUIT BREAKER! Features disabled after " << consecutive_hangs << " hangs" << std::endl;
            PERF_LOG_CRITICAL("Emergency circuit breaker triggered - features disabled");
        }
    } else if (frame_time_ < 0.02) { // Good performance
        consecutive_hangs = std::max(0, consecutive_hangs - 1); // Gradually recover
    }
    
    updateAdaptivePerformance();
    processInput();
    
    // Update audio sonification based on current camera position
    if (audio_enabled_) {
        PERF_TIMER("Audio Update");
        updateAudioSonification();
    }
    
    // Clear the screen
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Use appropriate shader based on optimization settings
    if (use_optimized_rendering_ && blackhole_optimized_shader_) {
        blackhole_optimized_shader_->use();
    } else {
        blackhole_shader_->use();
    }
    
    // Update uniforms
    updateUniforms();
    
    // Render fullscreen quad
    renderFullscreenQuad();
    
    // Render UI (if needed)
    renderUI();
    
    // Check for OpenGL errors only when there are performance issues
    static int error_check_counter = 0;
    if (frame_time_ > 0.1 || (++error_check_counter % 300 == 0)) { // Check on problems or every 5 seconds @ 60fps
        GLenum error = glGetError();
        if (error != GL_NO_ERROR) {
            PERF_LOG_ERROR("OpenGL error detected: " + std::to_string(error));
        }
    }
}

void BlackHoleRenderer::updateUniforms() {
    const Vector4& cam_pos = camera_->getPosition();
    
    // Get actual framebuffer size for gl_FragCoord calculations
    int framebuffer_width, framebuffer_height;
    glfwGetFramebufferSize(window_, &framebuffer_width, &framebuffer_height);
    
    // Get current shader (optimized or standard)
    Shader* current_shader = (use_optimized_rendering_ && blackhole_optimized_shader_) ? 
                             blackhole_optimized_shader_.get() : blackhole_shader_.get();
    
    // Set camera and scene uniforms
    current_shader->setVec3("cameraPos", 
        static_cast<float>(cam_pos.x), 
        static_cast<float>(cam_pos.y), 
        static_cast<float>(cam_pos.z));
    current_shader->setVec3("blackholePos", 0.0f, 0.0f, 0.0f);
    current_shader->setFloat("schwarzschildRadius", 
        static_cast<float>(physics_->schwarzschildRadius()));
    current_shader->setVec2("resolution", 
        static_cast<float>(framebuffer_width), 
        static_cast<float>(framebuffer_height));
    
    // Set debug mode (0=normal, 1=test pattern, 2=solid color)
    current_shader->setInt("debugMode", render_mode_);
    
    // Enhanced light simulation uniforms
    current_shader->setFloat("time", static_cast<float>(glfwGetTime()));
    current_shader->setInt("enableSpectralRendering", 1);
    current_shader->setInt("showLightRays", show_light_rays_ ? 1 : 0);
    
    // Accretion disk physics parameters - enhanced for realistic light display
    current_shader->setFloat("diskTemperature", 2e7f);  // 20 million Kelvin - hotter and more visible
    
    // Calculate gravitational redshift factor at camera position
    double camera_distance = VectorMath::magnitude3D(cam_pos);
    double rs = physics_->schwarzschildRadius();
    double redshift_factor = 1.0;
    if (camera_distance > rs) {
        redshift_factor = sqrt(1.0 - rs / camera_distance);
    }
    current_shader->setFloat("gravitationalRedshift", static_cast<float>(redshift_factor));
    
    // Relativistic jet uniforms
    current_shader->setInt("showJets", show_jets_ ? 1 : 0);
    current_shader->setFloat("jetLorentzFactor", jet_lorentz_factor_);
    
    // Rotate jet axis based on time for dynamic animation
    double current_time = glfwGetTime();
    double rotation_speed = 0.2; // Slow rotation for realistic appearance
    Vector4 rotated_jet_axis = Vector4(jet_axis_.t,
        jet_axis_.x,
        std::cos(current_time * rotation_speed),
        std::sin(current_time * rotation_speed));
    
    // Update the jet model with the new axis for physics consistency
    RelativisticJetModel jet_model = jets_->getJetModel();
    jet_model.jet_axis = rotated_jet_axis;
    jets_->setJetModel(jet_model);
    
    current_shader->setVec3("jetAxis", 
        static_cast<float>(rotated_jet_axis.x),
        static_cast<float>(rotated_jet_axis.y), 
        static_cast<float>(rotated_jet_axis.z));
    current_shader->setFloat("jetOpeningAngle", jet_opening_angle_);
    current_shader->setFloat("jetMagneticField", jet_magnetic_field_);
    
    // Enhanced realistic light parameters - new additions for perpetual light
    current_shader->setFloat("photonSphereGlow", 0.3f);  // Persistent glow from trapped light
    current_shader->setFloat("ambientLightLevel", 0.05f); // Subtle ambient lighting
    current_shader->setFloat("diskBrightness", 1.8f);    // Increased disk brightness
    current_shader->setFloat("starBrightness", 1.2f);    // Enhanced star visibility
    current_shader->setInt("enablePhotonSphereGlow", 1);  // Always show photon sphere effects
    
    // Enhanced physics-based proximity factor calculation
    // Effects should start much earlier (at ~10 Schwarzschild radii) and intensify dramatically
    double distance_ratio = camera_distance / rs;
    float proximity_factor;
    
    if (distance_ratio > 10.0) {
        // Very far - no effects
        proximity_factor = 0.0f;
    } else if (distance_ratio > 5.0) {
        // Early stage - subtle effects (10Rs to 5Rs)
        proximity_factor = static_cast<float>((10.0 - distance_ratio) / 50.0); // 0.0 to 0.1
    } else if (distance_ratio > 3.0) {
        // Medium stage - noticeable effects (5Rs to 3Rs)  
        proximity_factor = static_cast<float>(0.1 + (5.0 - distance_ratio) * 0.15); // 0.1 to 0.4
    } else if (distance_ratio > 2.0) {
        // Close stage - dramatic effects (3Rs to 2Rs)
        proximity_factor = static_cast<float>(0.4 + (3.0 - distance_ratio) * 0.3); // 0.4 to 0.7
    } else if (distance_ratio > 1.5) {
        // Very close - extreme effects (2Rs to 1.5Rs)
        proximity_factor = static_cast<float>(0.7 + (2.0 - distance_ratio) * 0.2); // 0.7 to 0.8
    } else if (distance_ratio > 1.0) {
        // Photon sphere to event horizon (1.5Rs to 1Rs) - maximum effects
        double photon_sphere_factor = (1.5 - distance_ratio) / 0.5; // 0 to 1
        proximity_factor = static_cast<float>(0.8 + photon_sphere_factor * 0.2); // 0.8 to 1.0
    } else {
        // Inside event horizon
        proximity_factor = 1.0f;
    }
    
    proximity_factor = std::max(0.0f, std::min(1.0f, proximity_factor));
    
    current_shader->setFloat("proximityFactor", proximity_factor);
    current_shader->setInt("enableDramaticEffects", 1);  // ALWAYS ENABLED
    
    // Debug proximity calculation and log dramatic effects
    static float last_logged_proximity = -1.0f;
    static int debug_counter = 0;
    
    // Debug output every 60 frames (once per second at 60fps)
    if (++debug_counter % 60 == 0) {
        std::cout << "🔍 DEBUG - Distance: " << camera_distance/rs << " Rs, Proximity Factor: " << proximity_factor 
                  << ", Dramatic Effects: " << (dramatic_effects_enabled_ ? "ON" : "OFF") << std::endl;
    }
    
    // Log dramatic effects when they become active (only if enabled)
    if (dramatic_effects_enabled_ && proximity_factor > 0.01 && std::abs(proximity_factor - last_logged_proximity) > 0.05) {
        if (proximity_factor > 0.8) {
            std::cout << "🔥 EXTREME PHYSICS: Event horizon firewall! Distance: " << std::fixed << std::setprecision(2) 
                     << distance_ratio << "Rs, Proximity: " << proximity_factor << std::endl;
        } else if (proximity_factor > 0.7) {
            std::cout << "🌪️ INTENSE WARPING: Photon sphere approach! Distance: " << std::fixed << std::setprecision(2) 
                     << distance_ratio << "Rs, Proximity: " << proximity_factor << std::endl;
        } else if (proximity_factor > 0.4) {
            std::cout << "🎭 DRAMATIC LENSING: Strong gravitational effects! Distance: " << std::fixed << std::setprecision(2) 
                     << distance_ratio << "Rs, Proximity: " << proximity_factor << std::endl;
        } else if (proximity_factor > 0.1) {
            std::cout << "🌊 SPACETIME RIPPLES: Tidal effects beginning! Distance: " << std::fixed << std::setprecision(2) 
                     << distance_ratio << "Rs, Proximity: " << proximity_factor << std::endl;
        } else {
            std::cout << "✨ EARLY EFFECTS: Subtle spacetime curvature! Distance: " << std::fixed << std::setprecision(2) 
                     << distance_ratio << "Rs, Proximity: " << proximity_factor << std::endl;
        }
        last_logged_proximity = proximity_factor;
    }
    
    // Performance optimization uniforms (for optimized shader)
    if (use_optimized_rendering_ && blackhole_optimized_shader_) {
        current_shader->setInt("lodLevel", current_lod_level_);
        current_shader->setFloat("performanceScale", performance_scale_);
        
        // Light ray optimization parameters
        double camera_distance = VectorMath::magnitude3D(cam_pos);
        double rs = physics_->schwarzschildRadius();
        float distance_scale = std::min(1.0f, static_cast<float>((50.0 * rs) / camera_distance));
        current_shader->setFloat("lightRayDistanceScale", distance_scale);
        
        // Adaptive ray count based on performance
        int max_light_sources = std::max(1, std::min(4, static_cast<int>(performance_scale_ * distance_scale * 4.0f)));
        current_shader->setInt("maxLightSources", max_light_sources);
    }
}

void BlackHoleRenderer::renderFullscreenQuad() {
    glBindVertexArray(quad_VAO_);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindVertexArray(0);
}

void BlackHoleRenderer::renderUI() {
    // Simple debug info in console for now
    static int frame_counter = 0;
    if (++frame_counter % 60 == 0) {
        printDebugInfo();
    }
    
    // Render control panel if visible
    if (show_control_panel_ && ui_shader_) {
        renderControlPanel();
    }
}

void BlackHoleRenderer::swapBuffers() {
    glfwSwapBuffers(window_);
}

void BlackHoleRenderer::pollEvents() {
    glfwPollEvents();
}

void BlackHoleRenderer::updatePerformanceStats() {
    double current_time = glfwGetTime();
    frame_time_ = current_time - last_time_;
    last_time_ = current_time;
    
    // Log frame data to performance logger
    const Vector4& pos = camera_->getPosition();
    double distance = VectorMath::magnitude3D(pos);
    double rs = physics_->schwarzschildRadius();
    
    std::string frame_details = 
        "Dist:" + std::to_string(distance / rs) + "Rs " +
        "Mode:" + (use_optimized_rendering_ ? "OPT" : "STD") + " " +
        "Perf:" + std::to_string(static_cast<int>(performance_scale_ * 100)) + "% " +
        "LOD:" + std::to_string(current_lod_level_);
    
    frame_count_++;
    static double fps_timer = 0.0;
    fps_timer += frame_time_;
    
    if (fps_timer >= 1.0) {
        fps_ = frame_count_ / fps_timer;
        frame_count_ = 0;
        fps_timer = 0.0;
        
        // Log FPS summary
        PerformanceLogger::getInstance().logFrame(frame_time_, fps_, frame_details);
    }
}

void BlackHoleRenderer::printDebugInfo() {
    const Vector4& pos = camera_->getPosition();
    double distance = VectorMath::magnitude3D(pos);
    double rs = physics_->schwarzschildRadius();
    
    std::cout << "FPS: " << static_cast<int>(fps_) 
              << " | Frame Time: " << frame_time_ * 1000.0 << "ms"
              << " | Distance: " << distance / rs << " Rs"
              << " | Mode: " << (use_optimized_rendering_ ? "OPTIMIZED" : "STANDARD")
              << " | Performance: " << static_cast<int>(performance_scale_ * 100) << "%"
              << " | LOD: " << current_lod_level_
              << std::endl;
}

void BlackHoleRenderer::updateAdaptivePerformance() {
    // Check for SEVERE performance drops every frame (freezing prevention)
    static double last_position_distance = 0.0;
    double current_distance = VectorMath::magnitude3D(camera_->getPosition());
    bool camera_moving = std::abs(current_distance - last_position_distance) > 0.1;
    last_position_distance = current_distance;
    
    // CRITICAL EMERGENCY: More aggressive threshold for GPU hangs
    if (frame_time_ > 0.1) { // > 100ms = potential GPU hang (was 200ms)
        const Vector4& pos = camera_->getPosition();
        double distance = VectorMath::magnitude3D(pos) / physics_->schwarzschildRadius();
        
        std::string cause = "Distance:" + std::to_string(distance) + "Rs, " +
                           "LightRays:" + (show_light_rays_ ? "ON" : "OFF") + ", " +
                           "Jets:" + (show_jets_ ? "ON" : "OFF") + ", " +
                           "OptMode:" + (use_optimized_rendering_ ? "OPT" : "STD");
        
        PerformanceLogger::getInstance().logPerformanceDrop(frame_time_, "CRITICAL FREEZE - " + cause);
        
        performance_scale_ = 0.2f; // Set to trigger shader emergency mode (< 0.4)
        current_lod_level_ = 3; // Lowest quality
        show_light_rays_ = false; // Force disable
        show_jets_ = false; // Force disable
        if (!use_optimized_rendering_) {
            switchToOptimizedRendering();
        }
        std::cout << "🚨 CRITICAL FREEZE! Frame time: " << frame_time_*1000 << "ms" << std::endl;
        return;
    }
    
    // SEVERE EMERGENCY: Major performance drop (for lesser issues)
    if (frame_time_ > 0.05) { // > 50ms = major issue (lowered from 100ms) 
        const Vector4& pos = camera_->getPosition();
        double distance = VectorMath::magnitude3D(pos) / physics_->schwarzschildRadius();
        
        std::string cause = "Distance:" + std::to_string(distance) + "Rs, " +
                           "PerfScale:" + std::to_string(performance_scale_) + ", " +
                           "LOD:" + std::to_string(current_lod_level_);
        
        PerformanceLogger::getInstance().logPerformanceDrop(frame_time_, "MAJOR DROP - " + cause);
        
        performance_scale_ = 0.3f; // Set to trigger aggressive shader optimization
        current_lod_level_ = 3; // Lowest quality immediately
        if (!use_optimized_rendering_) {
            switchToOptimizedRendering();
        }
        std::cout << "🚨 MAJOR performance drop! Frame time: " << frame_time_*1000 << "ms" << std::endl;
        return;
    }
    
    // MODERATE PERFORMANCE DROP: Handle smaller issues
    if (frame_time_ > 0.03) { // > 30ms = noticeable lag
        const Vector4& pos = camera_->getPosition();
        double distance = VectorMath::magnitude3D(pos) / physics_->schwarzschildRadius();
        
        std::string cause = "Distance:" + std::to_string(distance) + "Rs, " +
                           "Moving:" + (camera_moving ? "YES" : "NO");
        
        PerformanceLogger::getInstance().logPerformanceDrop(frame_time_, "MODERATE DROP - " + cause);
        
        performance_scale_ = 0.5f; // Moderate quality reduction  
        current_lod_level_ = 2;
        if (!use_optimized_rendering_) {
            switchToOptimizedRendering();
        }
        std::cout << "⚠️ Performance drop! Frame time: " << frame_time_*1000 << "ms" << std::endl;
        return;
    }
    
    // More frequent updates during movement
    static int update_counter = 0;
    int update_frequency = camera_moving ? 2 : 5; // Update every 2 frames when moving, 5 when static
    if (++update_counter < update_frequency) return;
    update_counter = 0;
    
    // More aggressive thresholds during movement
    const float low_fps_threshold = camera_moving ? target_fps_ * 0.9f : target_fps_ * 0.8f;
    const float high_fps_threshold = camera_moving ? target_fps_ * 1.1f : target_fps_ * 1.2f;
    
    if (fps_ < low_fps_threshold) {
        // Performance is too low - reduce quality
        adjustPerformanceSettings();
        if (!use_optimized_rendering_) {
            switchToOptimizedRendering();
        }
    } else if (fps_ > high_fps_threshold && use_optimized_rendering_) {
        // We have headroom - can increase quality
        if (performance_scale_ < 1.0f) {
            performance_scale_ = std::min(1.0f, performance_scale_ + 0.1f);
        } else if (current_lod_level_ > 0) {
            current_lod_level_--;
        } else {
            // Switch back to standard rendering if we have lots of headroom
            if (fps_ > target_fps_ * 1.5f) {
                switchToStandardRendering();
            }
        }
        
        // RECOVERY: Re-enable features if performance is excellent for sustained period
        static int recovery_counter = 0;
        static bool features_disabled_by_emergency = false;
        
        if (!show_light_rays_ || !show_jets_) {
            features_disabled_by_emergency = true;
        }
        
        if (features_disabled_by_emergency && fps_ > target_fps_ * 1.2f && performance_scale_ >= 0.8f) {
            recovery_counter++;
            if (recovery_counter > 180) { // 3 seconds @ 60fps of good performance
                show_light_rays_ = true;
                show_jets_ = true;
                features_disabled_by_emergency = false;
                recovery_counter = 0;
                std::cout << "✅ RECOVERY: Light rays and jets RE-ENABLED (sustained good performance)" << std::endl;
                PERF_LOG_INFO("Features recovered after sustained good performance");
            }
        } else {
            recovery_counter = 0; // Reset if performance drops
        }

    }
}

void BlackHoleRenderer::adjustPerformanceSettings() {
    // ULTRA-EMERGENCY: Immediate response to severe performance drops
    if (frame_time_ > 1.0) { // > 1000ms = major freeze
        performance_scale_ = 0.05f; // Minimum possible
        current_lod_level_ = 3;     // Lowest quality
        std::cout << "🚨 CRITICAL FREEZE PREVENTION! Minimum settings applied" << std::endl;
        return;
    }
    
    // SEVERE EMERGENCY: Major performance drop
    if (frame_time_ > 0.2) { // > 200ms = severe issue
        performance_scale_ = 0.05f; // Jump straight to minimum
        current_lod_level_ = 3; // Force lowest LOD
        std::cout << "🚨 SEVERE FREEZE! Jumping to minimum settings" << std::endl;
        return;
    }
    
    // EMERGENCY: Very fast degradation for bad performance
    if (frame_time_ > 0.1) { // > 100ms
        performance_scale_ = std::max(0.05f, performance_scale_ * 0.2f); // Very aggressive drop
        current_lod_level_ = 3; // Force lowest LOD
        std::cout << "🚨 MAJOR performance drop! Frame time: " << frame_time_*1000 << "ms" << std::endl;
        return;
    }
    
    // Progressive quality reduction for moderate performance issues
    if (performance_scale_ > 0.5f) {
        performance_scale_ = std::max(0.05f, performance_scale_ - 0.25f); // More aggressive
    } else if (current_lod_level_ < 3) {
        current_lod_level_++;
    }
    

}

void BlackHoleRenderer::switchToOptimizedRendering() {
    if (!use_optimized_rendering_ && blackhole_optimized_shader_) {
        use_optimized_rendering_ = true;
        performance_scale_ = 0.8f; // Start conservative
        current_lod_level_ = 2;    // Medium quality
        std::cout << "Switched to OPTIMIZED rendering for better performance" << std::endl;
    }
}

void BlackHoleRenderer::switchToStandardRendering() {
    if (use_optimized_rendering_) {
        use_optimized_rendering_ = false;
        std::cout << "Switched to STANDARD rendering (high performance headroom)" << std::endl;
    }
}

void BlackHoleRenderer::cleanup() {
    if (quad_VAO_) {
        glDeleteVertexArrays(1, &quad_VAO_);
        glDeleteBuffers(1, &quad_VBO_);
    }
    
    if (ui_quad_VAO_) {
        glDeleteVertexArrays(1, &ui_quad_VAO_);
        glDeleteBuffers(1, &ui_quad_VBO_);
    }
    
    if (window_) {
        glfwDestroyWindow(window_);
        glfwTerminate();
    }
}

// Static callback functions
void BlackHoleRenderer::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    BlackHoleRenderer* renderer = static_cast<BlackHoleRenderer*>(glfwGetWindowUserPointer(window));
    
    if (action == GLFW_PRESS) {
        renderer->keys_[key] = true;
        
        if (key == GLFW_KEY_ESCAPE) {
            glfwSetWindowShouldClose(window, true);
        } else if (key >= GLFW_KEY_1 && key <= GLFW_KEY_3) {
            renderer->render_mode_ = key - GLFW_KEY_1;
        } else if (key == GLFW_KEY_L) {
            renderer->toggleLightRayVisualization();
            std::cout << "Light ray visualization: " << (renderer->getLightRayVisualization() ? "ON" : "OFF") << std::endl;
        } else if (key == GLFW_KEY_J) {
            renderer->toggleJetVisualization();
            std::cout << "Relativistic jets: " << (renderer->getJetVisualization() ? "ON" : "OFF") << std::endl;
        } else if (key == GLFW_KEY_U) {
            // Increase jet Lorentz factor
            renderer->jet_lorentz_factor_ = std::min(50.0f, renderer->jet_lorentz_factor_ + 2.0f);
            std::cout << "Jet Lorentz factor: " << renderer->jet_lorentz_factor_ << std::endl;
            
            // Update jet model
            RelativisticJetModel jet_model = renderer->jets_->getJetModel();
            jet_model.bulk_lorentz_factor = renderer->jet_lorentz_factor_;
            renderer->jets_->setJetModel(jet_model);
        } else if (key == GLFW_KEY_I) {
            // Decrease jet Lorentz factor  
            renderer->jet_lorentz_factor_ = std::max(1.1f, renderer->jet_lorentz_factor_ - 2.0f);
            std::cout << "Jet Lorentz factor: " << renderer->jet_lorentz_factor_ << std::endl;
            
            // Update jet model
            RelativisticJetModel jet_model = renderer->jets_->getJetModel();
            jet_model.bulk_lorentz_factor = renderer->jet_lorentz_factor_;
            renderer->jets_->setJetModel(jet_model);
        } else if (key == GLFW_KEY_O) {
            // Increase jet opening angle
            renderer->jet_opening_angle_ = std::min(0.5f, renderer->jet_opening_angle_ + 0.02f);
            std::cout << "Jet opening angle: " << (renderer->jet_opening_angle_ * 180.0f / 3.14159f) << " degrees" << std::endl;
            
            // Update jet model
            RelativisticJetModel jet_model = renderer->jets_->getJetModel();
            jet_model.opening_angle = renderer->jet_opening_angle_;
            renderer->jets_->setJetModel(jet_model);
        } else if (key == GLFW_KEY_P) {
            // Decrease jet opening angle
            renderer->jet_opening_angle_ = std::max(0.05f, renderer->jet_opening_angle_ - 0.02f);
            std::cout << "Jet opening angle: " << (renderer->jet_opening_angle_ * 180.0f / 3.14159f) << " degrees" << std::endl;
            
            // Update jet model
            RelativisticJetModel jet_model = renderer->jets_->getJetModel();
            jet_model.opening_angle = renderer->jet_opening_angle_;
            renderer->jets_->setJetModel(jet_model);
        } else if (key == GLFW_KEY_Q) {
            // Toggle optimization mode
            if (renderer->use_optimized_rendering_) {
                renderer->switchToStandardRendering();
            } else {
                renderer->switchToOptimizedRendering();
            }
        } else if (key == GLFW_KEY_R) {
            // Reset performance settings to defaults
            renderer->performance_scale_ = 1.0f;
            renderer->current_lod_level_ = 2;
            renderer->use_optimized_rendering_ = true;
            std::cout << "Performance settings reset to defaults" << std::endl;
        } else if (key == GLFW_KEY_T) {
            // TURBO mode - ultra performance for smooth movement
            renderer->performance_scale_ = 0.3f;
            renderer->current_lod_level_ = 3;
            renderer->use_optimized_rendering_ = true;
            renderer->show_light_rays_ = false;
            renderer->show_jets_ = false;
            std::cout << "🏎️ TURBO mode: Ultra performance (light rays/jets disabled)" << std::endl;
        } else if (key == GLFW_KEY_Y) {
            // QUALITY mode - re-enable all features
            renderer->performance_scale_ = 1.0f;
            renderer->current_lod_level_ = 0;
            renderer->use_optimized_rendering_ = false;
            renderer->show_light_rays_ = true;
            renderer->show_jets_ = true;
            std::cout << "✨ QUALITY mode: All features enabled" << std::endl;
        } else if (key == GLFW_KEY_E) {
            // Toggle control panel
            renderer->toggleControlPanel();
            std::cout << "Control panel: " << (renderer->show_control_panel_ ? "ON" : "OFF") << std::endl;
        } else if (key == GLFW_KEY_S) {
            // Toggle audio sonification
            renderer->toggleSonification();
            std::cout << "🎵 Blackhole sonification: " << (renderer->isSonificationEnabled() ? "ON" : "OFF") << std::endl;
        } else if (key == GLFW_KEY_EQUAL || key == GLFW_KEY_KP_ADD) {
            // Increase volume
            renderer->setSonificationVolume(renderer->audio_volume_ + 0.1);
        } else if (key == GLFW_KEY_MINUS || key == GLFW_KEY_KP_SUBTRACT) {
            // Decrease volume
            renderer->setSonificationVolume(renderer->audio_volume_ - 0.1);
        } else if (key == GLFW_KEY_M) {
            // Toggle spatial audio
            renderer->enableSpatialAudio(!renderer->spatial_audio_enabled_);
        } else if (key == GLFW_KEY_N) {
            // Toggle relativistic effects
            renderer->enableRelativisticAudioEffects(!renderer->relativistic_audio_effects_);
        } else if (key == GLFW_KEY_Z) {
            // DEBUG: Force emergency mode
            renderer->show_light_rays_ = false;
            renderer->show_jets_ = false;
            renderer->performance_scale_ = 0.05f;
            renderer->use_optimized_rendering_ = true;
            std::cout << "🔧 DEBUG: Emergency mode forced (Z key)" << std::endl;
        } else if (key == GLFW_KEY_X) {
            // DEBUG: Reset all features
            renderer->show_light_rays_ = true;
            renderer->show_jets_ = true;
            renderer->performance_scale_ = 1.0f;
            std::cout << "🔧 DEBUG: All features reset (X key)" << std::endl;
        } else if (key == GLFW_KEY_BACKSPACE) {
            // EMERGENCY REVERSE: Quick escape to safe distance
            double rs = renderer->physics_->schwarzschildRadius();
            double safe_distance = 10.0 * rs;  // Safe viewing distance
            Vector4 new_pos(0, 0, 0, safe_distance);
            renderer->camera_->setPosition(new_pos);
            std::cout << "🚀 EMERGENCY REVERSE! Moved to safe distance: " << safe_distance/rs << " Rs" << std::endl;
        } else if (key == GLFW_KEY_V) {
            // V key now shows dramatic effects status (always enabled)
            std::cout << "🎨 DRAMATIC VISUAL EFFECTS: ALWAYS ENABLED!" << std::endl;
            std::cout << "    ✨ Spacetime distortion, color shifts, and firewall effects are PERMANENTLY ACTIVE!" << std::endl;
        } else if (key == GLFW_KEY_F) {
            // TEST: Jump to close distance to test dramatic effects
            double rs = renderer->physics_->schwarzschildRadius();
            double test_distance = 1.5 * rs;  // Close enough to trigger effects
            Vector4 new_pos(0, 0, 0, test_distance);
            renderer->camera_->setPosition(new_pos);
            std::cout << "🧪 TEST: Jumped to " << test_distance/rs << " Rs to test dramatic effects!" << std::endl;
        } else if (key == GLFW_KEY_SPACE) {
            // Toggle auto-movement through the blackhole
            renderer->auto_movement_enabled_ = !renderer->auto_movement_enabled_;
            if (renderer->auto_movement_enabled_) {
                std::cout << "🚀 AUTO-MOVEMENT: Started moving through blackhole! Press SPACEBAR again to stop." << std::endl;
            } else {
                std::cout << "⏹️ AUTO-MOVEMENT: Stopped. Camera movement paused." << std::endl;
            }
        }
    } else if (action == GLFW_RELEASE) {
        renderer->keys_[key] = false;
    }
}

void BlackHoleRenderer::mouseCallback(GLFWwindow* window, double xpos, double ypos) {
    // Mouse movement disabled - camera stays centered on blackhole
    // No rotation controls needed
}

void BlackHoleRenderer::scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    BlackHoleRenderer* renderer = static_cast<BlackHoleRenderer*>(glfwGetWindowUserPointer(window));
    
    // Distance-based zoom: move camera closer/further from blackhole center
    const Vector4& current_pos = renderer->camera_->getPosition();
    double current_distance = VectorMath::magnitude3D(current_pos);
    double rs = renderer->physics_->schwarzschildRadius();
    
    // Enhanced distance-adaptive zoom for dramatic close approach
    double min_distance = 1.001 * rs;  // Minimum possible distance
    double max_distance = 1000.0 * rs; // Maximum distance for reference
    
    // Calculate zoom journey progress (0.0 = far away, 1.0 = at event horizon)
    double journey_progress = 1.0 - (current_distance - min_distance) / (max_distance - min_distance);
    journey_progress = std::max(0.0, std::min(1.0, journey_progress));
    
    // Create dramatic slowdown when past halfway point (0.5)
    double zoom_speed_multiplier = 1.0;
    if (journey_progress > 0.5) {
        // Exponential slowdown - becomes very slow near the black hole
        double past_halfway = (journey_progress - 0.5) / 0.5;  // 0.0 to 1.0
        zoom_speed_multiplier = 1.0 - (past_halfway * past_halfway * 0.95); // Up to 95% slower
        zoom_speed_multiplier = std::max(0.05, zoom_speed_multiplier); // Never slower than 5% of normal
        std::cout << "🐌 DRAMATIC ZOOM: " << (int)(zoom_speed_multiplier * 100) << "% speed at " << current_distance/rs << " Rs" << std::endl;
    }
    
    // Apply adaptive zoom factor
    double base_zoom_rate = 0.1 * zoom_speed_multiplier;  // Base 10% modified by distance
    double zoom_factor = 1.0 + (yoffset * base_zoom_rate);
    double new_distance = current_distance * zoom_factor;
    
    // Clamp distance to reasonable bounds  
    new_distance = std::max(min_distance, std::min(max_distance, new_distance));
    
    // Add warnings when approaching dangerous distances
    if (new_distance < 2.0 * rs && new_distance != current_distance) {
        if (new_distance < 1.5 * rs) {
            std::cout << "🚨 EXTREME PROXIMITY WARNING! Distance: " << new_distance/rs << " Rs - Very close to event horizon!" << std::endl;
        } else {
            std::cout << "⚠️  Close proximity warning: Distance: " << new_distance/rs << " Rs - Approaching event horizon" << std::endl;
        }
    }
    
    // Update camera position while keeping it centered on blackhole
    Vector4 new_pos(0, 0, 0, new_distance);
    renderer->camera_->setPosition(new_pos);
}

void BlackHoleRenderer::framebufferSizeCallback(GLFWwindow* window, int width, int height) {
    BlackHoleRenderer* renderer = static_cast<BlackHoleRenderer*>(glfwGetWindowUserPointer(window));
    
    // Get actual framebuffer size (different from window size on high-DPI displays)
    int framebuffer_width, framebuffer_height;
    glfwGetFramebufferSize(window, &framebuffer_width, &framebuffer_height);
    
    // Update window dimensions for aspect ratio calculation
    renderer->window_width_ = width;
    renderer->window_height_ = height;
    renderer->camera_->setAspectRatio(static_cast<float>(width) / height);
    
    // Set viewport to actual framebuffer size
    glViewport(0, 0, framebuffer_width, framebuffer_height);
    
    std::cout << "Resize - Window: " << width << "x" << height 
              << ", Framebuffer: " << framebuffer_width << "x" << framebuffer_height << std::endl;
}

void BlackHoleRenderer::setupControlPanelGeometry() {
    // Create a quad for UI rendering in normalized device coordinates
    float uiVertices[] = {
        // positions   // texCoords
        0.0f,  1.0f,   0.0f, 1.0f,
        0.0f,  0.0f,   0.0f, 0.0f,
        1.0f,  0.0f,   1.0f, 0.0f,
        
        0.0f,  1.0f,   0.0f, 1.0f,
        1.0f,  0.0f,   1.0f, 0.0f,
        1.0f,  1.0f,   1.0f, 1.0f
    };
    
    glGenVertexArrays(1, &ui_quad_VAO_);
    glGenBuffers(1, &ui_quad_VBO_);
    
    glBindVertexArray(ui_quad_VAO_);
    glBindBuffer(GL_ARRAY_BUFFER, ui_quad_VBO_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(uiVertices), &uiVertices, GL_STATIC_DRAW);
    
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    
    glBindVertexArray(0);
}

void BlackHoleRenderer::renderControlPanel() {
    if (!ui_shader_ || !ui_quad_VAO_) {
        std::cout << "WARNING: UI shader or VAO not available. ui_shader_: " << (ui_shader_ ? "OK" : "NULL") 
                  << ", ui_quad_VAO_: " << ui_quad_VAO_ << std::endl;
        return;
    }
    
    // Debug output
    static bool debug_printed = false;
    if (!debug_printed) {
        std::cout << "DEBUG: Rendering control panel at window size " << window_width_ << "x" << window_height_ << std::endl;
        std::cout << "DEBUG: UI shader program ID: " << ui_shader_->getProgramID() << std::endl;
        debug_printed = true;
    }
    
    // Save current OpenGL state
    GLboolean blendEnabled = glIsEnabled(GL_BLEND);
    GLboolean depthTestEnabled = glIsEnabled(GL_DEPTH_TEST);
    
    // Enable blending for transparency
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_DEPTH_TEST);
    
    ui_shader_->use();
    
    // Set up orthographic projection matrix for UI
    float projection[16] = {0};
    float left = 0.0f, right = static_cast<float>(window_width_);
    float bottom = 0.0f, top = static_cast<float>(window_height_);
    float near = -1.0f, far = 1.0f;
    
    projection[0] = 2.0f / (right - left);
    projection[5] = 2.0f / (top - bottom);
    projection[10] = -2.0f / (far - near);
    projection[12] = -(right + left) / (right - left);
    projection[13] = -(top + bottom) / (top - bottom);
    projection[14] = -(far + near) / (far - near);
    projection[15] = 1.0f;
    
    // Panel dimensions and position (top-right corner)
    float panelWidth = 350.0f;
    float panelHeight = 400.0f;
    float panelX = window_width_ - panelWidth - 20.0f;  // 20px from right edge
    float panelY = window_height_ - panelHeight - 20.0f;  // 20px from top (OpenGL coordinates)
    
    // Set up model matrix for panel positioning
    float model[16] = {0};
    model[0] = panelWidth;   // Scale X
    model[5] = panelHeight;  // Scale Y  
    model[10] = 1.0f;        // Scale Z
    model[12] = panelX;      // Translate X
    model[13] = panelY;      // Translate Y
    model[15] = 1.0f;        // W component
    
    // Pass matrices to shader
    glUniformMatrix4fv(glGetUniformLocation(ui_shader_->getProgramID(), "projection"), 1, GL_FALSE, projection);
    glUniformMatrix4fv(glGetUniformLocation(ui_shader_->getProgramID(), "model"), 1, GL_FALSE, model);
    
    // Render panel background
    ui_shader_->setVec4("color", 0.1f, 0.1f, 0.15f, 0.85f);  // Dark blue-gray with transparency
    ui_shader_->setFloat("alpha", 0.85f);
    ui_shader_->setInt("renderMode", 0);  // Solid color mode
    
    glBindVertexArray(ui_quad_VAO_);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    
    // Render panel border
    ui_shader_->setVec4("color", 0.3f, 0.3f, 0.4f, 0.9f);  // Lighter border
    ui_shader_->setFloat("alpha", 0.9f);
    ui_shader_->setInt("renderMode", 1);  // Border mode
    
    glDrawArrays(GL_TRIANGLES, 0, 6);
    
    glBindVertexArray(0);
    
    // Restore OpenGL state
    if (!blendEnabled) glDisable(GL_BLEND);
    if (depthTestEnabled) glEnable(GL_DEPTH_TEST);
    
    // Now render text content using simple colored rectangles as visual indicators
    renderControlPanelContent(panelX, panelY, panelWidth, panelHeight);
}

void BlackHoleRenderer::renderControlPanelContent(float x, float y, float width, float height) {
    if (!ui_shader_ || !ui_quad_VAO_) return;
    
    ui_shader_->use();
    
    // Use the SAME projection and coordinate system as the background panel
    float projection[16] = {0};
    float left = 0.0f, right = static_cast<float>(window_width_);
    float bottom = 0.0f, top = static_cast<float>(window_height_);
    float near = -1.0f, far = 1.0f;
    
    projection[0] = 2.0f / (right - left);
    projection[5] = 2.0f / (top - bottom);
    projection[10] = -2.0f / (far - near);
    projection[12] = -(right + left) / (right - left);
    projection[13] = -(top + bottom) / (top - bottom);
    projection[14] = -(far + near) / (far - near);
    projection[15] = 1.0f;
    
    glUniformMatrix4fv(glGetUniformLocation(ui_shader_->getProgramID(), "projection"), 1, GL_FALSE, projection);
    
    // Test with a big obvious indicator first - put it in the middle of the panel
    float testX = x + 50.0f;
    float testY = y + 50.0f;
    float testWidth = 100.0f;
    float testHeight = 50.0f;
    
    // Debug output
    static bool content_debug_printed = false;
    if (!content_debug_printed) {
        std::cout << "DEBUG: Control panel content at x=" << x << ", y=" << y << ", w=" << width << ", h=" << height << std::endl;
        std::cout << "DEBUG: Test rect at x=" << testX << ", y=" << testY << ", w=" << testWidth << ", h=" << testHeight << std::endl;
        content_debug_printed = true;
    }
    
    glBindVertexArray(ui_quad_VAO_);
    ui_shader_->setInt("renderMode", 0);  // Solid color mode
    ui_shader_->setFloat("alpha", 1.0f);
    
    // Render a bright test rectangle that should definitely be visible
    float model[16] = {0};
    model[0] = testWidth;   // Scale X
    model[5] = testHeight;  // Scale Y  
    model[10] = 1.0f;       // Scale Z
    model[12] = testX;      // Translate X
    model[13] = testY;      // Translate Y
    model[15] = 1.0f;       // W component
    
    glUniformMatrix4fv(glGetUniformLocation(ui_shader_->getProgramID(), "model"), 1, GL_FALSE, model);
    ui_shader_->setVec4("color", 1.0f, 0.0f, 0.0f, 1.0f);  // Bright red - should be very visible
    glDrawArrays(GL_TRIANGLES, 0, 6);
    
    // Add a second test rectangle below it
    model[13] = testY + 70.0f;  // Move down
    glUniformMatrix4fv(glGetUniformLocation(ui_shader_->getProgramID(), "model"), 1, GL_FALSE, model);
    ui_shader_->setVec4("color", 0.0f, 1.0f, 0.0f, 1.0f);  // Bright green
    glDrawArrays(GL_TRIANGLES, 0, 6);
    
    // Add status indicators with simpler layout
    float indicatorY = testY + 150.0f;
    float indicatorSize = 40.0f;
    float spacing = 50.0f;
    
    // Status indicators for each toggleable option
    struct ControlItem {
        const char* label;
        bool isActive;
        float r, g, b;
    };
    
    ControlItem items[] = {
        {"Light Rays", show_light_rays_, 1.0f, 1.0f, 0.0f},    // Bright yellow
        {"Jets", show_jets_, 0.0f, 1.0f, 1.0f},               // Bright cyan
        {"Optimized", use_optimized_rendering_, 0.0f, 1.0f, 0.0f}, // Bright green
        {"Render Mode", render_mode_ > 0, 1.0f, 0.0f, 1.0f}   // Bright magenta
    };
    
    for (int i = 0; i < 4; i++) {
        const auto& item = items[i];
        
        // Position indicators in a row
        float itemX = testX + (i * (indicatorSize + 10.0f));
        
        float itemModel[16] = {0};
        itemModel[0] = indicatorSize;
        itemModel[5] = indicatorSize;
        itemModel[10] = 1.0f;
        itemModel[12] = itemX;
        itemModel[13] = indicatorY;
        itemModel[15] = 1.0f;
        
        glUniformMatrix4fv(glGetUniformLocation(ui_shader_->getProgramID(), "model"), 1, GL_FALSE, itemModel);
        
        if (item.isActive) {
            ui_shader_->setVec4("color", item.r, item.g, item.b, 1.0f);
        } else {
            ui_shader_->setVec4("color", 0.3f, 0.3f, 0.3f, 1.0f);  // Dark gray when inactive
        }
        
        glDrawArrays(GL_TRIANGLES, 0, 6);
    }
    
    glBindVertexArray(0);
}

// ====== AUDIO SONIFICATION ======

void BlackHoleRenderer::initializeAudioSystem() {
    std::cout << "🎵 Initializing Blackhole Sonification System..." << std::endl;
    
    // Create sonification engine with physics systems
    sonification_ = std::make_unique<BlackHoleSonification>(*physics_, *light_system_, *jets_);
    
    // Configure initial audio settings
    sonification_->setVolumeScale(audio_volume_);
    sonification_->enableRelativisticEffects(relativistic_audio_effects_);
    sonification_->setSpatialAudio(spatial_audio_enabled_);
    
    // Configure melodious musical settings
    sonification_->setMusicalScale(MusicalScale::MAJOR_PENTATONIC);  // Pleasant, no dissonant intervals
    sonification_->setFrequencyMappingMethod(FrequencyMappingMethod::QUANTIZED);
    sonification_->enableMelodiousFrequencies(true);  // Enable musical frequency mapping
    
    // Enable gentle audio processing to reduce harsh frequencies
    sonification_->setGentleAudioProcessing();  // EQ, compression, and de-essing for pleasant sound
    
    // Create audio buffer for real-time playback
    audio_buffer_ = std::make_unique<BlackHoleAudioBuffer>(4096, 44100);
    
    std::cout << "🎵 Blackhole Sonification initialized successfully!" << std::endl;
    std::cout << "🎵 You can now HEAR what a blackhole sounds like!" << std::endl;
    std::cout << "🎵 Audio Controls:" << std::endl;
    std::cout << "    S: Toggle sonification on/off" << std::endl;
    std::cout << "    +/-: Adjust volume" << std::endl;
    std::cout << "    M: Toggle spatial audio" << std::endl;
}

void BlackHoleRenderer::updateAudioSonification() {
    if (!sonification_ || !audio_buffer_) return;
    
    // Create sonification context from current camera state
    Vector4 camera_pos = camera_->getPosition();
    Vector4 camera_vel(0, 0, 0, 0);  // Simplified: assume static observer
    
    SonificationContext context = sonification_->updateContext(camera_pos, camera_vel, last_time_);
    
    // Generate audio spectrum based on current physics
    BlackHoleAudioSpectrum spectrum = sonification_->generateAudioSpectrum(context);
    
    // Fill audio buffer if needed
    if (!audio_buffer_->isBufferReady()) {
        audio_buffer_->fillBuffer(*sonification_, context);
    }
    
    // Debug output and proximity warnings
    static double last_debug = 0.0;
    static double last_warning = 0.0;
    
    // Enhanced proximity warning system
    if (context.schwarzschild_ratio < 2.0 && last_time_ - last_warning > 2.0) {
        if (context.schwarzschild_ratio < 1.2) {
            std::cout << "🚨🎵 AUDIO ALERT: CRITICAL PROXIMITY! " << context.schwarzschild_ratio << " Rs - APPROACHING EVENT HORIZON!" << std::endl;
            // Increase audio volume for emergency alert
            if (sonification_) {
                sonification_->setVolumeScale(std::min(2.0, audio_volume_ * 1.5));
            }
        } else if (context.schwarzschild_ratio < 1.5) {
            std::cout << "⚠️🎵 AUDIO WARNING: Extreme proximity " << context.schwarzschild_ratio << " Rs - Intense gravitational effects!" << std::endl;
        } else {
            std::cout << "🔊 AUDIO NOTICE: Close approach " << context.schwarzschild_ratio << " Rs - Enhanced audio effects active" << std::endl;
        }
        last_warning = last_time_;
    }
    
    // Regular debug output
    if (last_time_ - last_debug > 3.0) {
        std::cout << "🎵 Audio: Distance=" << context.schwarzschild_ratio << "Rs"
                  << " | Freq=" << spectrum.dominant_frequency << "Hz"
                  << " | Vol=" << spectrum.total_amplitude * 100 << "%" << std::endl;
        last_debug = last_time_;
    }
    
    // Reset audio volume when moving away from danger zone
    if (context.schwarzschild_ratio > 3.0 && sonification_) {
        sonification_->setVolumeScale(audio_volume_);  // Reset to normal volume
    }
}

void BlackHoleRenderer::setSonificationVolume(double volume) {
    audio_volume_ = std::max(0.0, std::min(2.0, volume));  // Clamp 0-2x
    if (sonification_) {
        sonification_->setVolumeScale(audio_volume_);
    }
    std::cout << "🎵 Volume set to " << static_cast<int>(audio_volume_ * 100) << "%" << std::endl;
}

void BlackHoleRenderer::enableRelativisticAudioEffects(bool enabled) {
    relativistic_audio_effects_ = enabled;
    if (sonification_) {
        sonification_->enableRelativisticEffects(enabled);
    }
    std::cout << "🎵 Relativistic audio effects " << (enabled ? "ENABLED" : "DISABLED") << std::endl;
}

void BlackHoleRenderer::enableSpatialAudio(bool enabled) {
    spatial_audio_enabled_ = enabled;
    if (sonification_) {
        sonification_->setSpatialAudio(enabled);
    }
    std::cout << "🎵 Spatial audio " << (enabled ? "ENABLED" : "DISABLED") << std::endl;
}

std::string BlackHoleRenderer::getPhysicsAudioBreakdown() const {
    if (!sonification_) return "Audio system not initialized";
    
    Vector4 camera_pos = camera_->getPosition();
    Vector4 camera_vel(0, 0, 0, 0);
    SonificationContext context = sonification_->updateContext(camera_pos, camera_vel, last_time_);
    
    return sonification_->getPhysicsBreakdown(context);
}
