#include "RayTracer.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <memory>

BlackHoleRayTracer::BlackHoleRayTracer(const BlackHolePhysics& physics, int width, int height)
    : physics_(physics), width_(width), height_(height), max_iterations_(10000),
      min_step_size_(0.001), max_step_size_(0.1), accretion_disk_enabled_(true),
      spectral_rendering_(true) {
    
    // Create light simulation system
    light_simulation_ = std::make_unique<SimpleLightSystem>(physics_);
    
    // Accretion disk parameters - use ISCO for inner radius
    disk_inner_radius_ = physics_.innerStableCircularOrbit();    // ISCO (Kerr-aware)
    disk_outer_radius_ = 20.0 * physics_.outerHorizonRadius();   // Use outer horizon for Kerr
    disk_height_ = 0.1 * physics_.outerHorizonRadius();
    
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
    result.final_colors.reserve(rays.size());
    
    // Trace each ray
    for (auto& ray : rays) {
        bool escaped = traceSingleRay(ray);
        SpectralColor ray_color = calculateRayColor(ray, camera_position);
        
        result.ray_paths.push_back(ray.path);
        result.escaped_rays.push_back(escaped);
        result.captured_rays.push_back(!escaped && physics_.isInsideEventHorizon(
            VectorMath::magnitude3D(ray.position)));
        result.final_distances.push_back(VectorMath::magnitude3D(ray.position));
        result.final_colors.push_back(ray_color);
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
    return final_distance > 50.0 * physics_.outerHorizonRadius();
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

SpectralColor BlackHoleRayTracer::sampleBackground(const Vector4& direction, const Vector4& observer) const {
    if (!spectral_rendering_) {
        // Fallback to simple star field
        Vector4 normalized_dir = VectorMath::normalize3D(direction);
        double hash = std::sin(normalized_dir.x * 12.9898 + normalized_dir.y * 78.233 + normalized_dir.z * 37.719) * 43758.5453;
        hash = hash - std::floor(hash);
        
        if (hash > 0.995) {
            return SpectralColor(1.0, 1.0, 0.8);  // Bright white-blue
        } else if (hash > 0.99) {
            return SpectralColor(0.8, 0.7, 0.5);  // Yellow-white
        } else if (hash > 0.985) {
            return SpectralColor(0.6, 0.4, 0.3);  // Orange-red
        } else {
            return SpectralColor(0.05, 0.05, 0.1);  // Very dark blue
        }
    }
    
    Vector4 normalized_dir = VectorMath::normalize3D(direction);
    double hash = std::sin(normalized_dir.x * 12.9898 + normalized_dir.y * 78.233 + normalized_dir.z * 37.719) * 43758.5453;
    hash = hash - std::floor(hash);
    
    if (hash > 0.995) {
        // Bright star with realistic stellar spectrum
        double star_temp = 8000.0 + hash * 10000.0;  // Hot blue-white stars
        // Simple starlight color based on temperature
        return light_simulation_->temperatureToRGB(star_temp);
    } else if (hash > 0.99) {
        // Medium star - sun-like
        double star_temp = 5000.0 + hash * 3000.0;  // Yellow-white stars
        // Simple starlight color based on temperature
        return light_simulation_->temperatureToRGB(star_temp);
    } else if (hash > 0.985) {
        // Dim star - red dwarf
        double star_temp = 3000.0 + hash * 2000.0;  // Red stars
        // Simple starlight color based on temperature
        return light_simulation_->temperatureToRGB(star_temp);
    } else {
        // Deep space background
        return SpectralColor(0.02, 0.02, 0.05);  // Very dark blue
    }
}

SpectralColor BlackHoleRayTracer::sampleAccretionDisk(const Vector4& position, const Vector4& observer) const {
    if (!accretion_disk_enabled_) {
        return SpectralColor(0, 0, 0, 0);
    }
    
    if (spectral_rendering_ && light_simulation_) {
        // Use enhanced accretion disk model with proper physics
        return light_simulation_->accretionDiskColor(position, observer);
    }
    
    // Fallback to simple temperature gradient
    double r = std::sqrt(position.x * position.x + position.y * position.y + position.z * position.z);
    double height = std::abs(position.z);
    
    if (r >= disk_inner_radius_ && r <= disk_outer_radius_ && height <= disk_height_) {
        double temp_factor = disk_outer_radius_ / r;
        double red = std::min(1.0, temp_factor * 0.3);
        double green = std::min(1.0, temp_factor * 0.8);
        double blue = std::min(1.0, temp_factor * 1.0);
        
        double turbulence = 0.1 * std::sin(position.x * 10) * std::cos(position.y * 8);
        return SpectralColor(red + turbulence, green + turbulence, blue);
    }
    
    return SpectralColor(0, 0, 0, 0);
}

SpectralColor BlackHoleRayTracer::calculateRayColor(const Ray& ray, const Vector4& observer_position) const {
    double final_distance = VectorMath::magnitude3D(ray.position);
    
    // Ray captured by black hole - pure black
    if (physics_.isInsideEventHorizon(final_distance)) {
        return SpectralColor(0, 0, 0, 0);
    }
    
    // Check if ray intersected accretion disk along its path
    SpectralColor total_color(0, 0, 0, 0);
    bool found_disk_intersection = false;
    
    // Sample along the ray path for disk intersections
    if (ray.path.size() > 1) {
        for (size_t i = 1; i < ray.path.size(); ++i) {
            SpectralColor disk_color = sampleAccretionDisk(ray.path[i], observer_position);
            if (disk_color.r > 0.01 || disk_color.g > 0.01 || disk_color.b > 0.01) {
                total_color = total_color + disk_color;
                found_disk_intersection = true;
            }
        }
    }
    
    // If no disk intersection, sample background in final ray direction
    if (!found_disk_intersection) {
        Vector4 final_direction = VectorMath::normalize3D(ray.direction);
        total_color = sampleBackground(final_direction, observer_position);
    }
    
    // Apply gravitational lensing brightness enhancement
    double lensing_brightness = 1.0;
    if (final_distance > 2.0 * physics_.outerHorizonRadius()) {
        double horizon_radius = physics_.outerHorizonRadius();
        lensing_brightness = 1.0 + 0.1 * horizon_radius / (final_distance * final_distance);
    }
    
    return total_color * lensing_brightness;
}

bool BlackHoleRayTracer::checkTerminationConditions(const Ray& ray) const {
    double distance = VectorMath::magnitude3D(ray.position);
    
    // Ray captured by black hole
    if (physics_.isInsideEventHorizon(distance)) {
        return true;
    }
    
    // Ray escaped to infinity
    if (distance > 100.0 * physics_.outerHorizonRadius()) {
        return true;
    }
    
    return false;
}

Vector4 BlackHoleRayTracer::calculateInitialMomentum(const Vector4& position, const Vector4& direction) const {
    // For photons, we need to ensure the momentum satisfies the null geodesic condition:
    // g_μν p^μ p^ν = 0 in curved spacetime
    double r = position.x;
    
    if (physics_.isInsideEventHorizon(r)) {
        return Vector4(0, 0, 0, 0);
    }
    
    // Normalize direction vector
    Vector4 norm_dir = VectorMath::normalize3D(direction);
    
    // Get metric components at this position
    double gtt = physics_.metricGtt(r);
    double grr = physics_.metricGrr(r);
    double gtheta = physics_.metricGthetatheta(r);
    
    // Extract theta from position (assuming spherical coordinates in Vector4)
    // For proper coordinate handling, convert cartesian to spherical if needed
    double theta = (position.y != 0) ? position.y : Constants::PI / 2.0;  // Default to equatorial plane
    double gphi = physics_.metricGphiphi(r, theta);
    
    // Set spatial momentum components from normalized direction
    double pr = norm_dir.x;
    double ptheta = norm_dir.y;  
    double pphi = norm_dir.z;
    
    // Calculate spatial momentum squared: g_ij p^i p^j
    double spatial_momentum_sq = grr * pr * pr + 
                                gtheta * ptheta * ptheta + 
                                gphi * pphi * pphi;
    
    // From null condition: g_tt (p^t)² + spatial_momentum_sq = 0
    // Therefore: p^t = sqrt(-spatial_momentum_sq / g_tt)
    double pt_squared = -spatial_momentum_sq / gtt;
    
    // Ensure we don't take sqrt of negative number (numerical safety)
    if (pt_squared < 0.0) {
        pt_squared = 0.0;
    }
    
    double pt = std::sqrt(pt_squared);
    
    // Ensure reasonable energy scale (photon energy E ~ 1)
    double energy_scale = 1.0;
    if (pt > 0.0) {
        energy_scale = 1.0 / pt;  // Normalize to unit energy
        pt *= energy_scale;
        pr *= energy_scale;
        ptheta *= energy_scale;
        pphi *= energy_scale;
    }
    
    return Vector4(pt, pr, ptheta, pphi);
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
