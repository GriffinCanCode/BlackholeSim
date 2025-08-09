#include "Renderer.h"
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#else
#include <GL/gl.h>
#endif
#include <iostream>
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
      quad_VAO_(0), quad_VBO_(0), render_mode_(0), black_hole_mass_(1.0),
      frame_time_(0.0), fps_(0.0), last_time_(0.0), frame_count_(0),
      last_mouse_x_(0.0), last_mouse_y_(0.0), first_mouse_(true), mouse_captured_(false) {
    
    std::memset(keys_, false, sizeof(keys_));
}

BlackHoleRenderer::~BlackHoleRenderer() {
    cleanup();
}

bool BlackHoleRenderer::initialize() {
    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return false;
    }
    
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
    
    // On macOS, OpenGL functions are available directly
    // No need to load function pointers
    
    // Setup OpenGL
    if (!setupOpenGL()) {
        return false;
    }
    
    // Initialize physics and ray tracer
    physics_ = std::make_unique<BlackHolePhysics>(black_hole_mass_);
    raytracer_ = std::make_unique<BlackHoleRayTracer>(*physics_, window_width_, window_height_);
    
    // Initialize camera - centered on blackhole with fixed viewing angle
    double initial_distance = 15.0 * physics_->schwarzschildRadius();
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
    
    // Setup callbacks
    setupCallbacks();
    
    // Enable depth testing
    glEnable(GL_DEPTH_TEST);
    
    // Set initial time
    last_time_ = glfwGetTime();
    
    std::cout << "Black Hole Simulator initialized successfully!" << std::endl;
    std::cout << "Controls:" << std::endl;
    std::cout << "  Mouse Wheel: Zoom in/out" << std::endl;
    std::cout << "  ESC: Exit" << std::endl;
    std::cout << "  1-3: Change render mode" << std::endl;
    
    return true;
}

bool BlackHoleRenderer::setupOpenGL() {
    // Set viewport
    glViewport(0, 0, window_width_, window_height_);
    
    // Set clear color to space black
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    
    return true;
}

bool BlackHoleRenderer::loadShaders() {
    blackhole_shader_ = std::make_unique<Shader>();
    
    if (!blackhole_shader_->loadFromFiles("../shaders/blackhole.vert", "../shaders/blackhole.frag")) {
        std::cerr << "Failed to load black hole shaders" << std::endl;
        return false;
    }
    
    // Validate required uniforms exist
    if (!validateShaderUniforms()) {
        std::cerr << "Shader uniform validation failed" << std::endl;
        return false;
    }
    
    return true;
}

bool BlackHoleRenderer::validateShaderUniforms() {
    std::vector<std::string> requiredUniforms = {
        "cameraPos",
        "blackholePos", 
        "schwarzschildRadius",
        "resolution",
        "debugMode"
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
    // No movement controls - camera stays centered on blackhole
    // Only zoom functionality is handled via scroll wheel callback
}

void BlackHoleRenderer::render() {
    updatePerformanceStats();
    processInput();
    
    // Clear the screen
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Use black hole shader
    blackhole_shader_->use();
    
    // Update uniforms
    updateUniforms();
    
    // Render fullscreen quad
    renderFullscreenQuad();
    
    // Render UI (if needed)
    renderUI();
}

void BlackHoleRenderer::updateUniforms() {
    const Vector4& cam_pos = camera_->getPosition();
    
    // Set camera and scene uniforms
    blackhole_shader_->setVec3("cameraPos", 
        static_cast<float>(cam_pos.x), 
        static_cast<float>(cam_pos.y), 
        static_cast<float>(cam_pos.z));
    blackhole_shader_->setVec3("blackholePos", 0.0f, 0.0f, 0.0f);
    blackhole_shader_->setFloat("schwarzschildRadius", 
        static_cast<float>(physics_->schwarzschildRadius()));
    blackhole_shader_->setVec2("resolution", 
        static_cast<float>(window_width_), 
        static_cast<float>(window_height_));
    
    // Set debug mode (0=normal, 1=test pattern, 2=solid color)
    blackhole_shader_->setInt("debugMode", render_mode_);
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
    
    frame_count_++;
    static double fps_timer = 0.0;
    fps_timer += frame_time_;
    
    if (fps_timer >= 1.0) {
        fps_ = frame_count_ / fps_timer;
        frame_count_ = 0;
        fps_timer = 0.0;
    }
}

void BlackHoleRenderer::printDebugInfo() {
    const Vector4& pos = camera_->getPosition();
    double distance = VectorMath::magnitude3D(pos);
    double rs = physics_->schwarzschildRadius();
    
    std::cout << "FPS: " << static_cast<int>(fps_) 
              << " | Frame Time: " << frame_time_ * 1000.0 << "ms"
              << " | Distance: " << distance / rs << " Rs"
              << " | Position: (" << pos.x << ", " << pos.y << ", " << pos.z << ")"
              << std::endl;
}

void BlackHoleRenderer::cleanup() {
    if (quad_VAO_) {
        glDeleteVertexArrays(1, &quad_VAO_);
        glDeleteBuffers(1, &quad_VBO_);
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
    
    // Zoom factor based on scroll
    double zoom_factor = 1.0 + (yoffset * 0.1);  // 10% per scroll step
    double new_distance = current_distance * zoom_factor;
    
    // Clamp distance to reasonable bounds (don't get too close to event horizon or too far away)
    double min_distance = 3.0 * rs;   // Stay well outside event horizon
    double max_distance = 100.0 * rs; // Don't go too far away
    new_distance = std::max(min_distance, std::min(max_distance, new_distance));
    
    // Update camera position while keeping it centered on blackhole
    Vector4 new_pos(0, 0, 0, new_distance);
    renderer->camera_->setPosition(new_pos);
}

void BlackHoleRenderer::framebufferSizeCallback(GLFWwindow* window, int width, int height) {
    BlackHoleRenderer* renderer = static_cast<BlackHoleRenderer*>(glfwGetWindowUserPointer(window));
    
    renderer->window_width_ = width;
    renderer->window_height_ = height;
    renderer->camera_->setAspectRatio(static_cast<float>(width) / height);
    
    glViewport(0, 0, width, height);
}
