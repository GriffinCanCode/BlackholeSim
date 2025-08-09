#include "LightSystem.h"
#include <cmath>
#include <algorithm>

SimpleLightSystem::SimpleLightSystem(const BlackHolePhysics& physics)
    : physics_(physics) {
}

SpectralColor SimpleLightSystem::temperatureToRGB(double temperature) const {
    // Simple blackbody to RGB conversion
    if (temperature < 3000) {
        return SpectralColor(1.0, 0.3, 0.0);  // Red
    } else if (temperature < 5000) {
        double t = (temperature - 3000) / 2000.0;
        return SpectralColor(1.0, 0.3 + 0.5 * t, 0.2 * t);  // Red-yellow
    } else if (temperature < 6500) {
        double t = (temperature - 5000) / 1500.0;
        return SpectralColor(1.0, 0.8 + 0.2 * t, 0.2 + 0.6 * t);  // Yellow-white
    } else {
        double t = std::min(1.0, (temperature - 6500) / 3500.0);
        return SpectralColor(1.0 - 0.2 * t, 1.0, 1.0);  // White-blue
    }
}

SpectralColor SimpleLightSystem::wavelengthToRGB(double wavelength) const {
    // Clamp wavelength to visible spectrum
    wavelength = std::max(SpectralConstants::WAVELENGTH_MIN, 
                         std::min(SpectralConstants::WAVELENGTH_MAX, wavelength));
    
    double r = 0.0, g = 0.0, b = 0.0;
    
    // Simple wavelength to RGB mapping
    if (wavelength >= 580 && wavelength < 645) {
        r = (wavelength - 580) / (645 - 580);
    } else if (wavelength >= 645 && wavelength <= 750) {
        r = 1.0;
    }
    
    if (wavelength >= 490 && wavelength < 510) {
        g = (wavelength - 490) / (510 - 490);
    } else if (wavelength >= 510 && wavelength < 580) {
        g = 1.0;
    } else if (wavelength >= 580 && wavelength <= 645) {
        g = 1.0 - (wavelength - 580) / (645 - 580);
    }
    
    if (wavelength >= 380 && wavelength < 440) {
        b = 0.3 + 0.7 * (wavelength - 380) / (440 - 380);
    } else if (wavelength >= 440 && wavelength <= 490) {
        b = 1.0;
    } else if (wavelength > 490 && wavelength < 510) {
        b = 1.0 - (wavelength - 490) / (510 - 490);
    }
    
    // Intensity falloff at spectrum edges
    double intensity = 1.0;
    if (wavelength < 420) {
        intensity = 0.3 + 0.7 * (wavelength - 380) / (420 - 380);
    } else if (wavelength > 700) {
        intensity = 0.3 + 0.7 * (750 - wavelength) / (750 - 700);
    }
    
    return SpectralColor(r * intensity, g * intensity, b * intensity, intensity);
}

double SimpleLightSystem::gravitationalRedshift(double r_emission, double r_observer) const {
    double rs = physics_.schwarzschildRadius();
    
    // Avoid division by zero near event horizon
    if (r_emission <= rs || r_observer <= rs) return 0.0;
    
    // Simplified gravitational redshift
    double factor_emission = 1.0 - rs / r_emission;
    double factor_observer = 1.0 - rs / r_observer;
    
    if (factor_emission <= 0.0 || factor_observer <= 0.0) return 0.0;
    
    return std::sqrt(factor_observer / factor_emission);
}

double SimpleLightSystem::dopplerShift(const Vector4& velocity, const Vector4& photon_direction) const {
    // Simplified Doppler shift calculation
    double beta_magnitude = VectorMath::magnitude3D(velocity) / SpectralConstants::SPEED_OF_LIGHT;
    beta_magnitude = std::min(beta_magnitude, 0.99);  // Avoid superluminal velocities
    
    double gamma = 1.0 / std::sqrt(1.0 - beta_magnitude * beta_magnitude);
    
    Vector4 normalized_velocity = VectorMath::normalize3D(velocity);
    Vector4 normalized_direction = VectorMath::normalize3D(photon_direction);
    double beta_dot_n = VectorMath::dot3D(normalized_velocity, normalized_direction);
    
    return gamma * (1.0 + beta_dot_n);
}

SpectralColor SimpleLightSystem::accretionDiskColor(const Vector4& position, const Vector4& observer) const {
    double r = std::sqrt(position.x * position.x + position.y * position.y);
    double height = std::abs(position.z);
    
    double rs = physics_.schwarzschildRadius();
    double inner_radius = 3.0 * rs;  // ISCO
    double outer_radius = 20.0 * rs;
    double disk_height = 0.1 * rs;
    
    // Check if position is within accretion disk bounds
    if (r < inner_radius || r > outer_radius || height > disk_height) {
        return SpectralColor(0, 0, 0, 0);
    }
    
    // Calculate local temperature
    double temperature = diskTemperature(r);
    SpectralColor thermal_color = temperatureToRGB(temperature);
    
    // Height falloff
    double height_factor = std::exp(-height * height / (disk_height * disk_height));
    thermal_color = thermal_color * height_factor;
    
    return thermal_color;
}

double SimpleLightSystem::planckFunction(double temperature, double wavelength) const {
    // Simplified Planck function
    double h = SpectralConstants::PLANCK_CONSTANT;
    double c = SpectralConstants::SPEED_OF_LIGHT;
    double k = SpectralConstants::BOLTZMANN_CONSTANT;
    
    double lambda5 = std::pow(wavelength, 5);
    double exp_term = std::exp(h * c / (wavelength * k * temperature));
    
    if (exp_term > 1e10) return 0.0;  // Avoid numerical overflow
    
    return (2.0 * h * c * c / lambda5) / (exp_term - 1.0);
}

double SimpleLightSystem::diskTemperature(double radius) const {
    double rs = physics_.schwarzschildRadius();
    double inner_radius = 3.0 * rs;
    double central_temperature = 3e7;  // 30 million Kelvin
    
    if (radius <= inner_radius) return central_temperature;
    
    // Temperature profile: T ∝ r^(-3/4)
    return central_temperature * std::pow(inner_radius / radius, 0.75);
}
