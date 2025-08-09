#pragma once
#include "Shader.h"
#include "Physics.h"
#include "RayTracer.h"
#include "LightSystem.h"
#include "RelativisticJets.h"
#include "BlackHoleSonification.h"

// Forward declaration for GLFW (actual header included in .cpp file)
struct GLFWwindow;

#include <memory>

class Camera {
public:
    Camera(const Vector4& position = Vector4(0, 0, 0, 10),
           const Vector4& target = Vector4(0, 0, 0, 0));
    
    // Camera controls
    void setPosition(const Vector4& position) { position_ = position; updateViewMatrix(); }
    void setTarget(const Vector4& target) { target_ = target; updateViewMatrix(); }
    void setFOV(float fov_degrees) { fov_ = fov_degrees; updateProjectionMatrix(); }
    
    // Movement
    void moveForward(float distance);
    void moveRight(float distance);
    void moveUp(float distance);
    void rotate(float yaw_delta, float pitch_delta);
    void orbit(float yaw_delta, float pitch_delta, float distance);
    
    // Getters
    const Vector4& getPosition() const { return position_; }
    const Vector4& getTarget() const { return target_; }
    float getFOV() const { return fov_; }
    const float* getViewMatrix() const { return view_matrix_; }
    const float* getProjectionMatrix() const { return projection_matrix_; }
    
    void setAspectRatio(float aspect) { aspect_ratio_ = aspect; updateProjectionMatrix(); }
    
private:
    Vector4 position_;
    Vector4 target_;
    float fov_;
    float aspect_ratio_;
    float near_plane_;
    float far_plane_;
    
    float view_matrix_[16];
    float projection_matrix_[16];
    
    void updateViewMatrix();
    void updateProjectionMatrix();
    void matrixLookAt(float* result, const Vector4& eye, const Vector4& center, const Vector4& up);
    void matrixPerspective(float* result, float fovy, float aspect, float near, float far);
};

class BlackHoleRenderer {
public:
    BlackHoleRenderer(int window_width = 1200, int window_height = 800);
    ~BlackHoleRenderer();
    
    // Initialization
    bool initialize();
    void cleanup();
    
    // Main render loop
    bool shouldClose() const;
    void processInput();
    void render();
    void swapBuffers();
    void pollEvents();
    
    // Settings
    void setBlackHoleMass(double mass_solar_units);
    void setCameraPosition(const Vector4& position);
    void setCameraTarget(const Vector4& target);
    void setRenderMode(int mode) { render_mode_ = mode; }
    void toggleLightRayVisualization() { show_light_rays_ = !show_light_rays_; }
    bool getLightRayVisualization() const { return show_light_rays_; }
    void toggleJetVisualization() { show_jets_ = !show_jets_; }
    bool getJetVisualization() const { return show_jets_; }
    void setJetLorentzFactor(float gamma) { jet_lorentz_factor_ = gamma; }
    void setJetOpeningAngle(float angle) { jet_opening_angle_ = angle; }
    
    // Audio sonification controls
    void toggleSonification() { audio_enabled_ = !audio_enabled_; }
    bool isSonificationEnabled() const { return audio_enabled_; }
    void setSonificationVolume(double volume);
    void enableRelativisticAudioEffects(bool enabled);
    void enableSpatialAudio(bool enabled);
    std::string getPhysicsAudioBreakdown() const;
    
    // Performance monitoring
    double getFrameTime() const { return frame_time_; }
    double getFPS() const { return fps_; }
    
    // Window management
    GLFWwindow* getWindow() const { return window_; }
    int getWindowWidth() const { return window_width_; }
    int getWindowHeight() const { return window_height_; }
    
private:
    // Window and OpenGL context
    GLFWwindow* window_;
    int window_width_, window_height_;
    
    // Rendering components
    std::unique_ptr<Shader> blackhole_shader_;
    std::unique_ptr<Shader> blackhole_optimized_shader_;
    std::unique_ptr<BlackHolePhysics> physics_;
    std::unique_ptr<BlackHoleRayTracer> raytracer_;

    std::unique_ptr<SimpleLightSystem> light_system_;
    std::unique_ptr<RelativisticJets> jets_;
    std::unique_ptr<BlackHoleSonification> sonification_;
    std::unique_ptr<BlackHoleAudioBuffer> audio_buffer_;
    std::unique_ptr<Camera> camera_;
    
    // Render data
    unsigned int quad_VAO_, quad_VBO_;
    
    // Settings
    int render_mode_;  // 0: Ray tracing, 1: Wireframe, 2: Debug
    double black_hole_mass_;
    bool show_light_rays_;
    bool show_jets_;
    float jet_lorentz_factor_;
    float jet_opening_angle_;
    float jet_magnetic_field_;
    Vector4 jet_axis_;
    
    // Audio sonification settings
    bool audio_enabled_;
    double audio_volume_;
    bool relativistic_audio_effects_;
    bool spatial_audio_enabled_;
    
    // Auto-movement settings
    bool auto_movement_enabled_;
    
    // Performance optimization settings
    bool use_optimized_rendering_;
    bool dramatic_effects_enabled_; // Toggle for enhanced visual effects
    float target_fps_;
    float performance_scale_; // 0.1-1.0 for shader performance scaling
    int current_lod_level_;   // 0-3 level of detail
    
    // Performance tracking
    double frame_time_;
    double fps_;
    double last_time_;
    int frame_count_;
    
    // Input handling
    bool keys_[1024];
    double last_mouse_x_, last_mouse_y_;
    bool first_mouse_;
    bool mouse_captured_;
    
    // Control panel
    bool show_control_panel_;
    unsigned int ui_quad_VAO_, ui_quad_VBO_;
    std::unique_ptr<Shader> ui_shader_;
    
    // Setup functions
    bool setupOpenGL();
    bool loadShaders();
    bool validateShaderUniforms();
    void setupGeometry();
    void setupCallbacks();
    
    // Rendering functions
    void renderFullscreenQuad();
    void updateUniforms();
    void renderUI();
    
    // Control panel functions
    void setupControlPanelGeometry();
    void renderControlPanel();
    void renderControlPanelContent(float x, float y, float width, float height);
    void toggleControlPanel() { show_control_panel_ = !show_control_panel_; }
    
    // Input callbacks (static functions for GLFW)
    static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
    static void mouseCallback(GLFWwindow* window, double xpos, double ypos);
    static void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);
    static void framebufferSizeCallback(GLFWwindow* window, int width, int height);
    
    // Utility functions
    void updatePerformanceStats();
    void printDebugInfo();
    
    // Audio sonification functions
    void updateAudioSonification();
    void initializeAudioSystem();
    
    // Performance optimization functions
    void updateAdaptivePerformance();
    void adjustPerformanceSettings();
    void switchToOptimizedRendering();
    void switchToStandardRendering();
};
