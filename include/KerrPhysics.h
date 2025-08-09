#pragma once
#include "Physics.h"
#include <array>
#include <cmath>

// Forward declarations
struct Vector4;

namespace KerrConstants {
    // Critical radii relationships for Kerr black holes
    constexpr double MAX_SPIN_PARAMETER = 0.998;  // Maximum stable a/M ratio
    constexpr double ERGOSPHERE_SAFETY_FACTOR = 1.001;
}

// Kerr spacetime coordinates and parameters
struct KerrParameters {
    double mass;                    // Black hole mass (solar masses)
    double angular_momentum;        // Dimensionless spin parameter a = J/(Mc)
    double outer_horizon_radius;    // r+ = M + √(M² - a²)
    double inner_horizon_radius;    // r- = M - √(M² - a²)
    double ergosphere_radius;       // r_ergo = M + √(M² - a²cos²θ)
    double static_limit;            // g_tt = 0 surface
    
    KerrParameters(double m = 1.0, double a = 0.0) 
        : mass(m), angular_momentum(std::max(-KerrConstants::MAX_SPIN_PARAMETER, 
                                           std::min(KerrConstants::MAX_SPIN_PARAMETER, a))) {
        updateDerivedParameters();
    }
    
    void updateDerivedParameters();
};

// Enhanced 4-vector operations for Kerr spacetime
struct KerrVector4 : public Vector4 {
    KerrVector4(double t = 0, double x = 0, double y = 0, double z = 0) 
        : Vector4(t, x, y, z) {}
    
    KerrVector4(const Vector4& v) : Vector4(v) {}
    
    // Boyer-Lindquist specific operations
    double getR() const { return x; }           // Radial coordinate
    double getTheta() const { return y; }       // Polar angle
    double getPhi() const { return z; }         // Azimuthal angle
    double getT() const { return t; }           // Time coordinate
    
    void setR(double r) { x = r; }
    void setTheta(double theta) { y = theta; }
    void setPhi(double phi) { z = phi; }
    void setT(double time) { t = time; }
    
    // Arithmetic operations
    KerrVector4 operator+(const KerrVector4& other) const {
        return KerrVector4(t + other.t, x + other.x, y + other.y, z + other.z);
    }
    
    KerrVector4 operator*(double scalar) const {
        return KerrVector4(t * scalar, x * scalar, y * scalar, z * scalar);
    }
};

class KerrPhysics {
public:
    KerrPhysics(double mass_solar_units = 1.0, double spin_parameter = 0.0);
    
    // Core Kerr metric properties
    const KerrParameters& getParameters() const { return params_; }
    void setSpinParameter(double a);
    void setMass(double mass_solar_units);
    
    // Kerr metric components in Boyer-Lindquist coordinates
    double metricGtt(double r, double theta) const;
    double metricGrr(double r) const;
    double metricGthetatheta(double r) const;
    double metricGphiphi(double r, double theta) const;
    double metricGtphi(double r, double theta) const;  // Off-diagonal term for frame-dragging
    
    // Derived metric quantities
    double metricDeterminant(double r, double theta) const;
    double deltaFunction(double r) const;              // Δ = r² - 2Mr + a²
    double sigmaFunction(double r, double theta) const; // Σ = r² + a²cos²θ
    double aFunction(double r, double theta) const;     // A = (r² + a²)² - a²Δsin²θ
    
    // Critical surfaces and radii
    double outerHorizonRadius() const { return params_.outer_horizon_radius; }
    double innerHorizonRadius() const { return params_.inner_horizon_radius; }
    double ergosphereRadius(double theta) const;
    double staticLimit(double theta) const;
    bool isInErgosphere(double r, double theta) const;
    bool isInsideOuterHorizon(double r) const { return r < params_.outer_horizon_radius; }
    
    // Frame-dragging and angular momentum effects
    double frameReferenceVelocity(double r, double theta) const;  // ω = -g_tφ/g_φφ
    double lenseThirringPrecession(double r, double theta) const;
    Vector4 applyFrameDragging(const Vector4& momentum, const Vector4& position) const;
    
    // Geodesic equations for Kerr spacetime
    KerrVector4 geodesicDerivativeKerr(const KerrVector4& position, const KerrVector4& momentum) const;
    std::pair<KerrVector4, KerrVector4> integrateKerrGeodesic(
        const KerrVector4& pos, const KerrVector4& mom, double step_size) const;
    
    // Conserved quantities in Kerr spacetime
    double specificEnergy(const KerrVector4& position, const KerrVector4& momentum) const;
    double specificAngularMomentum(const KerrVector4& position, const KerrVector4& momentum) const;
    double carterConstant(const KerrVector4& position, const KerrVector4& momentum) const;
    
    // Orbital mechanics for Kerr spacetime
    double innerMostStableCircularOrbit(bool prograde = true) const;
    double photonSphereRadius(bool prograde = true) const;
    double circularOrbitFrequency(double r, bool prograde = true) const;
    
    // Enhanced coordinate transformations
    KerrVector4 cartesianToBoyer(const Vector4& cartesian) const;
    Vector4 boyerToCartesian(const KerrVector4& boyer) const;
    KerrVector4 schwarzschildToBoyer(const Vector4& schwarzschild) const;
    
    // Ray tracing utilities specific to Kerr
    double adaptiveKerrStepSize(double r, double theta) const;
    bool checkKerrTermination(const KerrVector4& position) const;
    
    // Gravitational effects enhanced by rotation
    double deflectionAngleKerr(double impact_parameter, double inclination) const;
    double gravitationalRedshift(double r, double theta, const Vector4& observer_velocity) const;
    
    // Integration with existing BlackHolePhysics for compatibility
    BlackHolePhysics createEquivalentSchwarzschildPhysics() const;
    bool isEffectivelyNonRotating(double tolerance = 1e-6) const;
    
private:
    KerrParameters params_;
    
    // Internal computation helpers
    std::array<double, 16> christoffelSymbolsKerr(const KerrVector4& pos) const;
    double christoffelKerr(const KerrVector4& pos, int mu, int nu, int lambda) const;
    
    // Geodesic equation components
    double dtdlambda_derivative(const KerrVector4& pos, const KerrVector4& mom) const;
    double drdlambda_derivative(const KerrVector4& pos, const KerrVector4& mom) const;
    double dthetadlambda_derivative(const KerrVector4& pos, const KerrVector4& mom) const;
    double dphidlambda_derivative(const KerrVector4& pos, const KerrVector4& mom) const;
    
    // Validation helpers
    void validateParameters() const;
    void updateCachedQuantities();
    
    // Cached calculations for performance
    mutable double cached_delta_;
    mutable double cached_sigma_;
    mutable double cached_r_;
    mutable double cached_theta_;
    mutable bool cache_valid_;
    
    void invalidateCache() const { cache_valid_ = false; }
    void updateCache(double r, double theta) const;
};

// Utility functions for Kerr spacetime calculations
namespace KerrUtils {
    // Physical interpretations
    double extractionEfficiency(double spin_parameter);
    double maximallyExtractableEnergy(double mass, double spin_parameter);
    Vector4 ergosphereNormalVector(double r, double theta, double spin_parameter);
    
    // Coordinate system utilities
    bool isValidBoyerLindquist(double r, double theta, double phi);
    double coordinateTransformationJacobian(const KerrVector4& boyer_coords);
    
    // Specialized calculations
    double effectiveGravitationalRadius(double r, double theta, double spin_parameter);
    double rotationalEnergyDensity(const KerrParameters& params);
    double totalAngularMomentum(const KerrParameters& params);
    
    // Compatibility with existing Vector4 operations
    KerrVector4 convert(const Vector4& v);
    Vector4 convert(const KerrVector4& kv);
}
