#pragma once
#include <array>
#include <vector>
#include <cmath>
#include <memory>

// Forward declarations
class KerrPhysics;
class IBlackHolePhysics;

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

// Facade class providing backward compatibility while using the factory pattern internally
class BlackHolePhysics {
public:
    BlackHolePhysics(double mass_solar_units = 1.0, double spin_parameter = 0.0);
    ~BlackHolePhysics();  // Defined in cpp to ensure IBlackHolePhysics is complete
    
    // Core physics calculations - delegated to internal physics engine
    double schwarzschildRadius() const;
    double photonSphere() const;
    double innerStableCircularOrbit() const;
    
    // Kerr metric support
    bool isRotating() const;
    double getSpinParameter() const;
    void setSpinParameter(double a);
    
    // Mass management
    double getMass() const;
    void setMass(double mass_solar_units);
    
    // Enhanced properties for rotating black holes
    double outerHorizonRadius() const;
    double innerHorizonRadius() const;
    double ergosphereRadius(double theta = 0.0) const;
    bool isInErgosphere(double r, double theta) const;
    
    // Metric components
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
    bool isInsideEventHorizon(double r) const;
    
    // Gravitational lensing calculations
    double deflectionAngle(double impact_parameter) const;
    
    // Factory access for advanced usage
    const IBlackHolePhysics* getPhysicsEngine() const { return physics_engine_.get(); }
    std::string getPhysicsType() const;
    
private:
    std::unique_ptr<IBlackHolePhysics> physics_engine_;  // Factory-created physics implementation
    double mass_;                                        // Cached for quick access
    double spin_parameter_;                              // Cached for quick access
    
    // Internal helpers
    void updatePhysicsEngine();
    void syncParametersFromEngine();
};

// Utility functions for vector operations
namespace VectorMath {
    double dot3D(const Vector4& a, const Vector4& b);
    Vector4 cross3D(const Vector4& a, const Vector4& b);
    double magnitude3D(const Vector4& v);
    Vector4 normalize3D(const Vector4& v);
}
