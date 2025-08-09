#pragma once
#include "Physics.h"
#include <memory>
#include <string>
#include <functional>
#include <unordered_map>
#include <vector>

// Abstract interface for all black hole physics implementations
class IBlackHolePhysics {
public:
    virtual ~IBlackHolePhysics() = default;
    
    // Core physics calculations - pure virtual interface
    virtual double schwarzschildRadius() const = 0;
    virtual double photonSphere() const = 0;
    virtual double innerStableCircularOrbit() const = 0;
    
    // Enhanced properties for rotating black holes
    virtual double outerHorizonRadius() const = 0;
    virtual double innerHorizonRadius() const = 0;
    virtual double ergosphereRadius(double theta = 0.0) const = 0;
    virtual bool isInErgosphere(double r, double theta) const = 0;
    
    // Metric components
    virtual double metricGtt(double r) const = 0;
    virtual double metricGrr(double r) const = 0;
    virtual double metricGthetatheta(double r) const = 0;
    virtual double metricGphiphi(double r, double theta) const = 0;
    
    // Geodesic equations for photon trajectories
    virtual Vector4 geodesicDerivative(const Vector4& position, const Vector4& momentum) const = 0;
    
    // Runge-Kutta 4th order integration for ray tracing
    virtual std::pair<Vector4, Vector4> integrateGeodesic(
        const Vector4& pos, const Vector4& mom, double step_size) const = 0;
    
    // Coordinate transformations
    virtual Vector4 cartesianToSpherical(const Vector4& cartesian) const = 0;
    virtual Vector4 sphericalToCartesian(const Vector4& spherical) const = 0;
    
    // Ray tracing utilities
    virtual double adaptiveStepSize(double r) const = 0;
    virtual bool isInsideEventHorizon(double r) const = 0;
    
    // Gravitational lensing calculations
    virtual double deflectionAngle(double impact_parameter) const = 0;
    
    // Physics properties
    virtual double getMass() const = 0;
    virtual double getSpinParameter() const = 0;
    virtual bool isRotating() const = 0;
    
    // Type identification
    virtual std::string getPhysicsType() const = 0;
};

// Schwarzschild (non-rotating) black hole implementation
class SchwarzchildPhysics : public IBlackHolePhysics {
public:
    SchwarzchildPhysics(double mass_solar_units);
    virtual ~SchwarzchildPhysics() = default;
    
    // IBlackHolePhysics interface implementation
    double schwarzschildRadius() const override { return schwarzschild_radius_; }
    double photonSphere() const override;
    double innerStableCircularOrbit() const override;
    
    double outerHorizonRadius() const override { return schwarzschild_radius_; }
    double innerHorizonRadius() const override { return 0.0; }
    double ergosphereRadius(double theta = 0.0) const override { return schwarzschild_radius_; }
    bool isInErgosphere(double r, double theta) const override { return false; }
    
    double metricGtt(double r) const override;
    double metricGrr(double r) const override;
    double metricGthetatheta(double r) const override;
    double metricGphiphi(double r, double theta) const override;
    
    Vector4 geodesicDerivative(const Vector4& position, const Vector4& momentum) const override;
    std::pair<Vector4, Vector4> integrateGeodesic(
        const Vector4& pos, const Vector4& mom, double step_size) const override;
    
    Vector4 cartesianToSpherical(const Vector4& cartesian) const override;
    Vector4 sphericalToCartesian(const Vector4& spherical) const override;
    
    double adaptiveStepSize(double r) const override;
    bool isInsideEventHorizon(double r) const override;
    double deflectionAngle(double impact_parameter) const override;
    
    double getMass() const override { return mass_; }
    double getSpinParameter() const override { return 0.0; }
    bool isRotating() const override { return false; }
    std::string getPhysicsType() const override { return "Schwarzschild"; }

private:
    double mass_;
    double schwarzschild_radius_;
    
    // Helper functions for geodesic calculations
    std::array<double, 4> christoffelSymbols(const Vector4& pos, int mu, int nu, int lambda) const;
    double christoffel(const Vector4& pos, int mu, int nu, int lambda) const;
};

// Kerr (rotating) black hole implementation
class KerrPhysicsEngine : public IBlackHolePhysics {
public:
    KerrPhysicsEngine(double mass_solar_units, double spin_parameter);
    virtual ~KerrPhysicsEngine() = default;
    
    // IBlackHolePhysics interface implementation
    double schwarzschildRadius() const override { return 2.0 * mass_; }
    double photonSphere() const override;
    double innerStableCircularOrbit() const override;
    
    double outerHorizonRadius() const override;
    double innerHorizonRadius() const override;
    double ergosphereRadius(double theta = 0.0) const override;
    bool isInErgosphere(double r, double theta) const override;
    
    double metricGtt(double r) const override;
    double metricGrr(double r) const override;
    double metricGthetatheta(double r) const override;
    double metricGphiphi(double r, double theta) const override;
    
    Vector4 geodesicDerivative(const Vector4& position, const Vector4& momentum) const override;
    std::pair<Vector4, Vector4> integrateGeodesic(
        const Vector4& pos, const Vector4& mom, double step_size) const override;
    
    Vector4 cartesianToSpherical(const Vector4& cartesian) const override;
    Vector4 sphericalToCartesian(const Vector4& spherical) const override;
    
    double adaptiveStepSize(double r) const override;
    bool isInsideEventHorizon(double r) const override;
    double deflectionAngle(double impact_parameter) const override;
    
    double getMass() const override { return mass_; }
    double getSpinParameter() const override { return spin_parameter_; }
    bool isRotating() const override { return std::abs(spin_parameter_) > 1e-6; }
    std::string getPhysicsType() const override { return "Kerr"; }
    
    // Kerr-specific methods
    void setSpinParameter(double a);

private:
    double mass_;
    double spin_parameter_;
    std::unique_ptr<KerrPhysics> kerr_physics_;  // Delegate to existing Kerr implementation
    
    void updateKerrPhysics();
};

// Factory class for creating physics engines
class BlackHolePhysicsFactory {
public:
    enum class PhysicsType {
        Schwarzschild,
        Kerr,
        Auto  // Automatically choose based on spin parameter
    };
    
    // Factory methods
    static std::unique_ptr<IBlackHolePhysics> create(
        double mass_solar_units, 
        double spin_parameter = 0.0, 
        PhysicsType type = PhysicsType::Auto
    );
    
    static std::unique_ptr<IBlackHolePhysics> createSchwarzschild(double mass_solar_units);
    static std::unique_ptr<IBlackHolePhysics> createKerr(double mass_solar_units, double spin_parameter);
    
    // Utility functions
    static PhysicsType determinePhysicsType(double spin_parameter);
    static std::string getPhysicsTypeName(PhysicsType type);
    
private:
    static constexpr double SPIN_THRESHOLD = 1e-6;
};

// Extension point for future physics models
class BlackHolePhysicsRegistry {
public:
    using PhysicsCreator = std::function<std::unique_ptr<IBlackHolePhysics>(double, double)>;
    
    static void registerPhysicsType(const std::string& name, PhysicsCreator creator);
    static std::unique_ptr<IBlackHolePhysics> createByName(
        const std::string& name, double mass, double spin = 0.0);
    static std::vector<std::string> getAvailableTypes();
    
private:
    static std::unordered_map<std::string, PhysicsCreator> creators_;
};
