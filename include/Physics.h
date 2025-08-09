#pragma once
#include <array>
#include <vector>
#include <cmath>

// Physical constants (in geometric units where c = G = 1)
namespace Constants {
    constexpr double G = 6.67430e-11;           // Gravitational constant (m³/kg⋅s²)
    constexpr double c = 299792458.0;           // Speed of light (m/s)
    constexpr double SOLAR_MASS = 1.989e30;     // Solar mass (kg)
    constexpr double PI = 3.14159265358979323846;
}

// 4-vector for spacetime coordinates and momentum
struct Vector4 {
    double t, x, y, z;
    
    Vector4(double t = 0, double x = 0, double y = 0, double z = 0) 
        : t(t), x(x), y(y), z(z) {}
    
    Vector4 operator+(const Vector4& other) const {
        return Vector4(t + other.t, x + other.x, y + other.y, z + other.z);
    }
    
    Vector4 operator*(double scalar) const {
        return Vector4(t * scalar, x * scalar, y * scalar, z * scalar);
    }
};

// Spherical coordinates for easier black hole physics
struct SphericalCoords {
    double r, theta, phi;
    
    SphericalCoords(double r = 0, double theta = 0, double phi = 0)
        : r(r), theta(theta), phi(phi) {}
};

class BlackHolePhysics {
public:
    BlackHolePhysics(double mass_solar_units = 1.0);
    
    // Core physics calculations
    double schwarzschildRadius() const { return schwarzschild_radius_; }
    double photonSphere() const { return 1.5 * schwarzschild_radius_; }
    double innerStableCircularOrbit() const { return 3.0 * schwarzschild_radius_; }
    
    // Schwarzschild metric components
    double metricGtt(double r) const;
    double metricGrr(double r) const;
    double metricGthetatheta(double r) const;
    double metricGphiphi(double r, double theta) const;
    
    // Geodesic equations for photon trajectories
    Vector4 geodesicDerivative(const Vector4& position, const Vector4& momentum) const;
    
    // Runge-Kutta 4th order integration for ray tracing
    std::pair<Vector4, Vector4> integrateGeodesic(
        const Vector4& pos, const Vector4& mom, double step_size) const;
    
    // Coordinate transformations
    Vector4 cartesianToSpherical(const Vector4& cartesian) const;
    Vector4 sphericalToCartesian(const Vector4& spherical) const;
    
    // Ray tracing utilities
    double adaptiveStepSize(double r) const;
    bool isInsideEventHorizon(double r) const { return r < schwarzschild_radius_; }
    
    // Gravitational lensing calculations
    double deflectionAngle(double impact_parameter) const;
    
private:
    double mass_;                    // Black hole mass in solar masses
    double schwarzschild_radius_;   // Event horizon radius
    
    // Helper functions for geodesic calculations
    std::array<double, 4> christoffelSymbols(const Vector4& pos, int mu, int nu, int lambda) const;
    double christoffel(const Vector4& pos, int mu, int nu, int lambda) const;
};

// Utility functions for vector operations
namespace VectorMath {
    double dot3D(const Vector4& a, const Vector4& b);
    Vector4 cross3D(const Vector4& a, const Vector4& b);
    double magnitude3D(const Vector4& v);
    Vector4 normalize3D(const Vector4& v);
}
