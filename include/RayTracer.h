#pragma once
#include "Physics.h"
#include "LightSystem.h"
#include <vector>
#include <array>

struct Ray {
    Vector4 origin;
    Vector4 direction;
    Vector4 position;
    Vector4 momentum;
    bool active;
    double affine_parameter;
    std::vector<Vector4> path;  // For visualization
    
    Ray(const Vector4& orig, const Vector4& dir) 
        : origin(orig), direction(dir), position(orig), active(true), affine_parameter(0.0) {
        // Initialize momentum for null geodesic (photon)
        momentum = Vector4(1.0, dir.x, dir.y, dir.z);  // E=1 for photons
        path.push_back(position);
    }
};

struct RayTracingResult {
    std::vector<std::vector<Vector4>> ray_paths;
    std::vector<bool> escaped_rays;
    std::vector<bool> captured_rays;
    std::vector<double> final_distances;
    std::vector<SpectralColor> final_colors;  // Enhanced with spectral color information
};

class BlackHoleRayTracer {
public:
    BlackHoleRayTracer(const BlackHolePhysics& physics, int width = 800, int height = 600);
    
    // Main ray tracing function
    RayTracingResult traceRays(const Vector4& camera_position, 
                              const Vector4& camera_target,
                              double fov_degrees = 45.0);
    
    // Trace a single ray through spacetime
    bool traceSingleRay(Ray& ray, double max_affine_parameter = 100.0);
    
    // Generate camera rays
    std::vector<Ray> generateCameraRays(const Vector4& camera_pos,
                                       const Vector4& camera_target,
                                       double fov_degrees);
    
    // Enhanced color sampling with light simulation
    SpectralColor sampleBackground(const Vector4& direction, const Vector4& observer) const;
    SpectralColor sampleAccretionDisk(const Vector4& position, const Vector4& observer) const;
    
    // Light simulation integration
    SpectralColor calculateRayColor(const Ray& ray, const Vector4& observer_position) const;
    
    // Settings
    void setMaxIterations(int max_iter) { max_iterations_ = max_iter; }
    void setMinStepSize(double min_step) { min_step_size_ = min_step; }
    void setMaxStepSize(double max_step) { max_step_size_ = max_step; }
    void setAccretionDiskEnabled(bool enabled) { accretion_disk_enabled_ = enabled; }
    void enableSpectralRendering(bool enabled) { spectral_rendering_ = enabled; }
    
    // Getters
    int getWidth() const { return width_; }
    int getHeight() const { return height_; }
    
private:
    const BlackHolePhysics& physics_;
    std::unique_ptr<SimpleLightSystem> light_simulation_;
    int width_, height_;
    int max_iterations_;
    double min_step_size_;
    double max_step_size_;
    bool accretion_disk_enabled_;
    bool spectral_rendering_;
    
    // Accretion disk parameters
    double disk_inner_radius_;
    double disk_outer_radius_;
    double disk_height_;
    
    // Helper functions
    bool checkTerminationConditions(const Ray& ray) const;
    Vector4 calculateInitialMomentum(const Vector4& position, const Vector4& direction) const;
    double calculateProperTime(const Vector4& position, const Vector4& momentum) const;
    
    // Background generation
    void generateStarField();
    std::vector<Vector4> star_positions_;
    std::vector<double> star_brightnesses_;
};

// Utility functions for camera setup
namespace CameraUtils {
    std::array<Vector4, 3> calculateCameraBasis(const Vector4& position, const Vector4& target);
    Vector4 screenToWorldDirection(double screen_x, double screen_y, 
                                  const std::array<Vector4, 3>& camera_basis,
                                  double fov_radians, double aspect_ratio);
}
