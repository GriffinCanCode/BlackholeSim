#include "RayTracer.h"
#include <cmath>
#include <random>
#include <algorithm>

BlackHoleRayTracer::BlackHoleRayTracer(const BlackHolePhysics& physics, int width, int height)
    : physics_(physics), width_(width), height_(height), max_iterations_(10000),
      min_step_size_(0.001), max_step_size_(0.1), accretion_disk_enabled_(true) {
    
    // Accretion disk parameters (in units of Schwarzschild radius)
    disk_inner_radius_ = 3.0 * physics_.schwarzschildRadius();  // ISCO
    disk_outer_radius_ = 20.0 * physics_.schwarzschildRadius();
    disk_height_ = 0.1 * physics_.schwarzschildRadius();
    
    generateStarField();
}

RayTracingResult BlackHoleRayTracer::traceRays(const Vector4& camera_position,
                                              const Vector4& camera_target,
                                              double fov_degrees) {
    RayTracingResult result;
    
    // Generate rays from camera
    std::vector<Ray> rays = generateCameraRays(camera_position, camera_target, fov_degrees);
    
    result.ray_paths.reserve(rays.size());
    result.escaped_rays.reserve(rays.size());
    result.captured_rays.reserve(rays.size());
    result.final_distances.reserve(rays.size());
    
    // Trace each ray
    for (auto& ray : rays) {
        bool escaped = traceSingleRay(ray);
        
        result.ray_paths.push_back(ray.path);
        result.escaped_rays.push_back(escaped);
        result.captured_rays.push_back(!escaped && physics_.isInsideEventHorizon(
            VectorMath::magnitude3D(ray.position)));
        result.final_distances.push_back(VectorMath::magnitude3D(ray.position));
    }
    
    return result;
}

bool BlackHoleRayTracer::traceSingleRay(Ray& ray, double max_affine_parameter) {
    // Convert to spherical coordinates for integration
    Vector4 spherical_pos = physics_.cartesianToSpherical(ray.position);
    
    // Initialize momentum in spherical coordinates
    Vector4 spherical_mom = calculateInitialMomentum(spherical_pos, ray.direction);
    
    int iterations = 0;
    double affine_param = 0.0;
    
    while (iterations < max_iterations_ && affine_param < max_affine_parameter) {
        // Check termination conditions
        if (checkTerminationConditions(ray)) {
            break;
        }
        
        double r = spherical_pos.x;
        
        // Adaptive step size based on distance from black hole
        double step_size = physics_.adaptiveStepSize(r);
        step_size = std::max(min_step_size_, std::min(max_step_size_, step_size));
        
        // Integrate geodesic
        auto [new_pos, new_mom] = physics_.integrateGeodesic(spherical_pos, spherical_mom, step_size);
        
        spherical_pos = new_pos;
        spherical_mom = new_mom;
        affine_param += step_size;
        
        // Convert back to Cartesian for path storage
        Vector4 cartesian_pos = physics_.sphericalToCartesian(spherical_pos);
        ray.position = cartesian_pos;
        ray.affine_parameter = affine_param;
        
        // Store path point for visualization (subsample to avoid too many points)
        if (iterations % 10 == 0) {
            ray.path.push_back(cartesian_pos);
        }
        
        iterations++;
    }
    
    // Ray escaped if it's far from the black hole
    double final_distance = VectorMath::magnitude3D(ray.position);
    return final_distance > 50.0 * physics_.schwarzschildRadius();
}

std::vector<Ray> BlackHoleRayTracer::generateCameraRays(const Vector4& camera_pos,
                                                       const Vector4& camera_target,
                                                       double fov_degrees) {
    std::vector<Ray> rays;
    rays.reserve(width_ * height_);
    
    // Calculate camera basis vectors
    auto camera_basis = CameraUtils::calculateCameraBasis(camera_pos, camera_target);
    double fov_radians = fov_degrees * Constants::PI / 180.0;
    double aspect_ratio = static_cast<double>(width_) / height_;
    
    // Generate ray for each pixel
    for (int y = 0; y < height_; ++y) {
        for (int x = 0; x < width_; ++x) {
            // Convert pixel coordinates to normalized screen coordinates [-1, 1]
            double screen_x = (2.0 * x / width_) - 1.0;
            double screen_y = 1.0 - (2.0 * y / height_);
            
            // Calculate ray direction
            Vector4 ray_direction = CameraUtils::screenToWorldDirection(
                screen_x, screen_y, camera_basis, fov_radians, aspect_ratio);
            
            rays.emplace_back(camera_pos, ray_direction);
        }
    }
    
    return rays;
}

Vector4 BlackHoleRayTracer::sampleBackground(const Vector4& direction) const {
    // Simple procedural star field
    Vector4 normalized_dir = VectorMath::normalize3D(direction);
    
    // Create pseudo-random stars based on direction
    double hash = std::sin(normalized_dir.x * 12.9898 + normalized_dir.y * 78.233 + normalized_dir.z * 37.719) * 43758.5453;
    hash = hash - std::floor(hash);
    
    if (hash > 0.995) {
        // Bright star
        return Vector4(0, 1.0, 1.0, 0.8);  // Bright white-blue
    } else if (hash > 0.99) {
        // Medium star
        return Vector4(0, 0.8, 0.7, 0.5);  // Yellow-white
    } else if (hash > 0.985) {
        // Dim star
        return Vector4(0, 0.6, 0.4, 0.3);  // Orange-red
    } else {
        // Dark space
        return Vector4(0, 0.05, 0.05, 0.1);  // Very dark blue
    }
}

Vector4 BlackHoleRayTracer::sampleAccretionDisk(const Vector4& position) const {
    if (!accretion_disk_enabled_) {
        return Vector4(0, 0, 0, 0);
    }
    
    double r = std::sqrt(position.x * position.x + position.y * position.y + position.z * position.z);
    double height = std::abs(position.z);
    
    // Check if position is within accretion disk
    if (r >= disk_inner_radius_ && r <= disk_outer_radius_ && height <= disk_height_) {
        // Temperature gradient: hotter closer to black hole
        double temp_factor = disk_outer_radius_ / r;
        
        // Blackbody radiation approximation
        double red = std::min(1.0, temp_factor * 0.3);
        double green = std::min(1.0, temp_factor * 0.8);
        double blue = std::min(1.0, temp_factor * 1.0);
        
        // Add some turbulence
        double turbulence = 0.1 * std::sin(position.x * 10) * std::cos(position.y * 8);
        
        return Vector4(0, red + turbulence, green + turbulence, blue);
    }
    
    return Vector4(0, 0, 0, 0);
}

bool BlackHoleRayTracer::checkTerminationConditions(const Ray& ray) const {
    double distance = VectorMath::magnitude3D(ray.position);
    
    // Ray captured by black hole
    if (physics_.isInsideEventHorizon(distance)) {
        return true;
    }
    
    // Ray escaped to infinity
    if (distance > 100.0 * physics_.schwarzschildRadius()) {
        return true;
    }
    
    return false;
}

Vector4 BlackHoleRayTracer::calculateInitialMomentum(const Vector4& position, const Vector4& direction) const {
    // For photons, we need to ensure the momentum satisfies the null condition
    // in Schwarzschild coordinates
    double r = position.x;
    
    if (r <= physics_.schwarzschildRadius()) {
        return Vector4(0, 0, 0, 0);
    }
    
    // Normalize direction
    Vector4 norm_dir = VectorMath::normalize3D(direction);
    
    // Energy component (conserved)
    double E = 1.0;  // Photon energy
    
    // Calculate momentum components ensuring null geodesic condition
    double gtt = physics_.metricGtt(r);
    double pt = E / (-gtt);  // dt/dλ
    
    return Vector4(pt, norm_dir.x, norm_dir.y, norm_dir.z);
}

void BlackHoleRayTracer::generateStarField() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    std::uniform_real_distribution<> brightness_dis(0.1, 1.0);
    
    int num_stars = 10000;
    star_positions_.reserve(num_stars);
    star_brightnesses_.reserve(num_stars);
    
    for (int i = 0; i < num_stars; ++i) {
        // Generate random direction on unit sphere
        Vector4 dir;
        do {
            dir = Vector4(0, dis(gen), dis(gen), dis(gen));
        } while (VectorMath::magnitude3D(dir) > 1.0);
        
        dir = VectorMath::normalize3D(dir);
        star_positions_.push_back(dir);
        star_brightnesses_.push_back(brightness_dis(gen));
    }
}

namespace CameraUtils {
    std::array<Vector4, 3> calculateCameraBasis(const Vector4& position, const Vector4& target) {
        Vector4 forward = VectorMath::normalize3D(Vector4(0, 
            target.x - position.x, 
            target.y - position.y, 
            target.z - position.z));
        
        // Assume up vector is along z-axis
        Vector4 world_up(0, 0, 0, 1);
        Vector4 right = VectorMath::normalize3D(VectorMath::cross3D(forward, world_up));
        Vector4 up = VectorMath::normalize3D(VectorMath::cross3D(right, forward));
        
        return {forward, right, up};
    }
    
    Vector4 screenToWorldDirection(double screen_x, double screen_y,
                                  const std::array<Vector4, 3>& camera_basis,
                                  double fov_radians, double aspect_ratio) {
        double tan_half_fov = std::tan(fov_radians * 0.5);
        
        Vector4 direction = camera_basis[0];  // Forward
        direction = direction + camera_basis[1] * (screen_x * tan_half_fov * aspect_ratio);  // Right
        direction = direction + camera_basis[2] * (screen_y * tan_half_fov);  // Up
        
        return VectorMath::normalize3D(direction);
    }
}
