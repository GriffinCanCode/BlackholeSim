#pragma once
#include "Physics.h"

// Simple spectral constants for basic light physics
namespace SpectralConstants {
    constexpr double PLANCK_CONSTANT = 6.62607015e-34;    // J⋅s
    constexpr double BOLTZMANN_CONSTANT = 1.380649e-23;   // J/K
    constexpr double SPEED_OF_LIGHT = Constants::c;       // m/s
    constexpr double WIEN_CONSTANT = 2.897771955e-3;      // m⋅K
    
    // Visible spectrum wavelengths (in nanometers)
    constexpr double WAVELENGTH_RED = 700.0;
    constexpr double WAVELENGTH_GREEN = 550.0;
    constexpr double WAVELENGTH_BLUE = 450.0;
    constexpr double WAVELENGTH_MIN = 380.0;
    constexpr double WAVELENGTH_MAX = 750.0;
}

// Simple RGB color structure for light simulation
struct SpectralColor {
    double r, g, b, intensity;
    
    SpectralColor(double red = 0.0, double green = 0.0, double blue = 0.0, double i = 1.0)
        : r(red), g(green), b(blue), intensity(i) {}
    
    SpectralColor operator+(const SpectralColor& other) const {
        return SpectralColor(r + other.r, g + other.g, b + other.b, 
                           std::max(intensity, other.intensity));
    }
    
    SpectralColor operator*(double scalar) const {
        return SpectralColor(r * scalar, g * scalar, b * scalar, intensity * scalar);
    }
};

// Simple light system for basic physics calculations
class SimpleLightSystem {
public:
    SimpleLightSystem(const BlackHolePhysics& physics);
    
    // Basic color calculations
    SpectralColor temperatureToRGB(double temperature) const;
    SpectralColor wavelengthToRGB(double wavelength) const;
    
    // Basic physics calculations
    double gravitationalRedshift(double r_emission, double r_observer) const;
    double dopplerShift(const Vector4& velocity, const Vector4& photon_direction) const;
    
    // Simple accretion disk color
    SpectralColor accretionDiskColor(const Vector4& position, const Vector4& observer) const;
    
private:
    const BlackHolePhysics& physics_;
    
    // Helper functions
    double planckFunction(double temperature, double wavelength) const;
    double diskTemperature(double radius) const;
};
