#include "Physics.h"
#include "PhysicsFactory.h"
#include <cmath>
#include <algorithm>

BlackHolePhysics::BlackHolePhysics(double mass_solar_units, double spin_parameter) 
    : mass_(mass_solar_units), spin_parameter_(spin_parameter) {
    updatePhysicsEngine();
}

BlackHolePhysics::~BlackHolePhysics() = default;

void BlackHolePhysics::updatePhysicsEngine() {
    // Use factory to create appropriate physics engine
    physics_engine_ = BlackHolePhysicsFactory::create(mass_, spin_parameter_);
    syncParametersFromEngine();
}

void BlackHolePhysics::syncParametersFromEngine() {
    if (physics_engine_) {
        mass_ = physics_engine_->getMass();
        spin_parameter_ = physics_engine_->getSpinParameter();
    }
}

// Delegated methods - all forward to the internal physics engine

double BlackHolePhysics::schwarzschildRadius() const {
    return physics_engine_ ? physics_engine_->schwarzschildRadius() : 0.0;
}

double BlackHolePhysics::photonSphere() const {
    return physics_engine_ ? physics_engine_->photonSphere() : 0.0;
}

double BlackHolePhysics::innerStableCircularOrbit() const {
    return physics_engine_ ? physics_engine_->innerStableCircularOrbit() : 0.0;
}

bool BlackHolePhysics::isRotating() const {
    return physics_engine_ ? physics_engine_->isRotating() : false;
}

double BlackHolePhysics::getSpinParameter() const {
    return spin_parameter_;
}

void BlackHolePhysics::setSpinParameter(double a) {
    spin_parameter_ = a;
    updatePhysicsEngine();
}

double BlackHolePhysics::getMass() const {
    return mass_;
}

void BlackHolePhysics::setMass(double mass_solar_units) {
    mass_ = mass_solar_units;
    updatePhysicsEngine();
}

double BlackHolePhysics::outerHorizonRadius() const {
    return physics_engine_ ? physics_engine_->outerHorizonRadius() : 0.0;
}

double BlackHolePhysics::innerHorizonRadius() const {
    return physics_engine_ ? physics_engine_->innerHorizonRadius() : 0.0;
}

double BlackHolePhysics::ergosphereRadius(double theta) const {
    return physics_engine_ ? physics_engine_->ergosphereRadius(theta) : 0.0;
}

bool BlackHolePhysics::isInErgosphere(double r, double theta) const {
    return physics_engine_ ? physics_engine_->isInErgosphere(r, theta) : false;
}

double BlackHolePhysics::metricGtt(double r) const {
    return physics_engine_ ? physics_engine_->metricGtt(r) : 0.0;
}

double BlackHolePhysics::metricGrr(double r) const {
    return physics_engine_ ? physics_engine_->metricGrr(r) : 0.0;
}

double BlackHolePhysics::metricGthetatheta(double r) const {
    return physics_engine_ ? physics_engine_->metricGthetatheta(r) : 0.0;
}

double BlackHolePhysics::metricGphiphi(double r, double theta) const {
    return physics_engine_ ? physics_engine_->metricGphiphi(r, theta) : 0.0;
}

Vector4 BlackHolePhysics::geodesicDerivative(const Vector4& position, const Vector4& momentum) const {
    return physics_engine_ ? physics_engine_->geodesicDerivative(position, momentum) : Vector4();
}

std::pair<Vector4, Vector4> BlackHolePhysics::integrateGeodesic(
    const Vector4& pos, const Vector4& mom, double step_size) const {
    return physics_engine_ ? 
        physics_engine_->integrateGeodesic(pos, mom, step_size) : 
        std::make_pair(pos, mom);
}

Vector4 BlackHolePhysics::cartesianToSpherical(const Vector4& cartesian) const {
    return physics_engine_ ? physics_engine_->cartesianToSpherical(cartesian) : cartesian;
}

Vector4 BlackHolePhysics::sphericalToCartesian(const Vector4& spherical) const {
    return physics_engine_ ? physics_engine_->sphericalToCartesian(spherical) : spherical;
}

double BlackHolePhysics::adaptiveStepSize(double r) const {
    return physics_engine_ ? physics_engine_->adaptiveStepSize(r) : 0.1;
}

bool BlackHolePhysics::isInsideEventHorizon(double r) const {
    return physics_engine_ ? physics_engine_->isInsideEventHorizon(r) : false;
}

double BlackHolePhysics::deflectionAngle(double impact_parameter) const {
    return physics_engine_ ? physics_engine_->deflectionAngle(impact_parameter) : 0.0;
}

std::string BlackHolePhysics::getPhysicsType() const {
    return physics_engine_ ? physics_engine_->getPhysicsType() : "None";
}

// VectorMath namespace implementation remains unchanged
namespace VectorMath {
    double dot3D(const Vector4& a, const Vector4& b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }
    
    Vector4 cross3D(const Vector4& a, const Vector4& b) {
        return Vector4(0,
                      a.y * b.z - a.z * b.y,
                      a.z * b.x - a.x * b.z,
                      a.x * b.y - a.y * b.x);
    }
    
    double magnitude3D(const Vector4& v) {
        return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }
    
    Vector4 normalize3D(const Vector4& v) {
        double mag = magnitude3D(v);
        if (mag > 0) {
            return Vector4(v.t, v.x / mag, v.y / mag, v.z / mag);
        }
        return v;
    }
}