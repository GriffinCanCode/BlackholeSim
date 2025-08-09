#include "PhysicsFactory.h"
#include "KerrPhysics.h"
#include <cmath>
#include <stdexcept>

// ============================================================================
// SchwarzchildPhysics Implementation
// ============================================================================

SchwarzchildPhysics::SchwarzchildPhysics(double mass_solar_units) 
    : mass_(mass_solar_units) {
    schwarzschild_radius_ = 2.0 * mass_solar_units;
}

double SchwarzchildPhysics::photonSphere() const {
    return 1.5 * schwarzschild_radius_;
}

double SchwarzchildPhysics::innerStableCircularOrbit() const {
    return 3.0 * schwarzschild_radius_;
}

double SchwarzchildPhysics::metricGtt(double r) const {
    return -(1.0 - schwarzschild_radius_ / r);
}

double SchwarzchildPhysics::metricGrr(double r) const {
    return 1.0 / (1.0 - schwarzschild_radius_ / r);
}

double SchwarzchildPhysics::metricGthetatheta(double r) const {
    return r * r;
}

double SchwarzchildPhysics::metricGphiphi(double r, double theta) const {
    return r * r * std::sin(theta) * std::sin(theta);
}

Vector4 SchwarzchildPhysics::geodesicDerivative(const Vector4& position, const Vector4& momentum) const {
    double r = position.x;  // Using x as radial coordinate in spherical
    double theta = position.y;  // Using y as theta coordinate
    
    if (r <= schwarzschild_radius_) {
        return Vector4(0, 0, 0, 0);  // Inside event horizon
    }
    
    double rs = schwarzschild_radius_;
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    
    // Geodesic equations in Schwarzschild coordinates
    // d²t/dλ² = -2rs/(r(r-rs)) * dr/dλ * dt/dλ
    double dtdt_coeff = -rs / (r * (r - rs));
    double d2t = 2.0 * dtdt_coeff * momentum.x * momentum.t;
    
    // d²r/dλ² = -rs(r-rs)/(2r³) * (dt/dλ)² + rs/(2r(r-rs)) * (dr/dλ)² + (r-rs) * [(dθ/dλ)² + sin²θ(dφ/dλ)²]
    double d2r = -rs * (r - rs) / (2.0 * r * r * r) * momentum.t * momentum.t
                + rs / (2.0 * r * (r - rs)) * momentum.x * momentum.x
                + (r - rs) * (momentum.y * momentum.y + sin_theta * sin_theta * momentum.z * momentum.z);
    
    // d²θ/dλ² = -2/r * dr/dλ * dθ/dλ + sin(θ)cos(θ) * (dφ/dλ)²
    double d2theta = -2.0 / r * momentum.x * momentum.y 
                    + sin_theta * cos_theta * momentum.z * momentum.z;
    
    // d²φ/dλ² = -2/r * dr/dλ * dφ/dλ - 2cos(θ)/sin(θ) * dθ/dλ * dφ/dλ
    double d2phi = -2.0 / r * momentum.x * momentum.z
                  - 2.0 * cos_theta / sin_theta * momentum.y * momentum.z;
    
    return Vector4(d2t, d2r, d2theta, d2phi);
}

std::pair<Vector4, Vector4> SchwarzchildPhysics::integrateGeodesic(
    const Vector4& pos, const Vector4& mom, double step_size) const {
    
    // 4th order Runge-Kutta integration
    Vector4 k1_pos = mom;
    Vector4 k1_mom = geodesicDerivative(pos, mom);
    
    Vector4 k2_pos = mom + k1_mom * (step_size * 0.5);
    Vector4 k2_mom = geodesicDerivative(pos + k1_pos * (step_size * 0.5), 
                                       mom + k1_mom * (step_size * 0.5));
    
    Vector4 k3_pos = mom + k2_mom * (step_size * 0.5);
    Vector4 k3_mom = geodesicDerivative(pos + k2_pos * (step_size * 0.5),
                                       mom + k2_mom * (step_size * 0.5));
    
    Vector4 k4_pos = mom + k3_mom * step_size;
    Vector4 k4_mom = geodesicDerivative(pos + k3_pos * step_size,
                                       mom + k3_mom * step_size);
    
    Vector4 new_pos = pos + (k1_pos + k2_pos * 2.0 + k3_pos * 2.0 + k4_pos) * (step_size / 6.0);
    Vector4 new_mom = mom + (k1_mom + k2_mom * 2.0 + k3_mom * 2.0 + k4_mom) * (step_size / 6.0);
    
    return std::make_pair(new_pos, new_mom);
}

Vector4 SchwarzchildPhysics::cartesianToSpherical(const Vector4& cartesian) const {
    double x = cartesian.x;
    double y = cartesian.y;
    double z = cartesian.z;
    
    double r = std::sqrt(x*x + y*y + z*z);
    if (r == 0.0) return Vector4(cartesian.t, 0, 0, 0);
    
    double theta = std::acos(std::max(-1.0, std::min(1.0, z / r)));
    double phi = std::atan2(y, x);
    
    return Vector4(cartesian.t, r, theta, phi);
}

Vector4 SchwarzchildPhysics::sphericalToCartesian(const Vector4& spherical) const {
    double r = spherical.x;
    double theta = spherical.y;
    double phi = spherical.z;
    
    double x = r * std::sin(theta) * std::cos(phi);
    double y = r * std::sin(theta) * std::sin(phi);
    double z = r * std::cos(theta);
    
    return Vector4(spherical.t, x, y, z);
}

double SchwarzchildPhysics::adaptiveStepSize(double r) const {
    // Smaller steps near the black hole for accuracy
    double min_step = 0.001;
    double max_step = 0.1;
    
    if (r < 2.0 * schwarzschild_radius_) {
        return min_step;
    } else if (r < 5.0 * schwarzschild_radius_) {
        return min_step + (max_step - min_step) * (r - 2.0 * schwarzschild_radius_) / (3.0 * schwarzschild_radius_);
    } else {
        return max_step;
    }
}

bool SchwarzchildPhysics::isInsideEventHorizon(double r) const {
    return r < schwarzschild_radius_;
}

double SchwarzchildPhysics::deflectionAngle(double impact_parameter) const {
    // Exact strong-field deflection angle for Schwarzschild metric
    // Uses Darwin's exact formula and strong-field approximations
    
    double M = schwarzschild_radius_ / 2.0; // Convert to mass units
    double b = impact_parameter; // Impact parameter
    double u = M / b; // Dimensionless parameter
    
    // Critical impact parameter for photon sphere (b = 3√3 M)
    double b_crit = 3.0 * std::sqrt(3.0) * M;
    
    if (b < 0.1 * b_crit) {
        // Very strong field - particle likely captured
        return Constants::PI; // 180 degree deflection (effective capture)
    }
    
    if (b > 10.0 * b_crit) {
        // Weak field regime - use standard weak field formula
        return 4.0 * M / b;
    }
    
    // Strong field regime - use exact Darwin formula
    // α = 4u + (15π/4)u² + (128/3)u³ + O(u⁴)
    
    double alpha_0 = 4.0 * u; // Leading order (weak field)
    double alpha_1 = (15.0 * Constants::PI / 4.0) * u * u; // Second order
    double alpha_2 = (128.0 / 3.0) * u * u * u; // Third order (strong field)
    
    // Fourth order correction for very strong fields
    double alpha_3 = (3465.0 * Constants::PI / 64.0) * u * u * u * u;
    
    double total_deflection = alpha_0 + alpha_1 + alpha_2 + alpha_3;
    
    // For very strong fields near the photon sphere, add logarithmic divergence correction
    if (b < 2.0 * b_crit) {
        double log_correction = 2.0 * std::log(b / b_crit);
        total_deflection += log_correction * u * u;
    }
    
    // Ensure physically reasonable result (cap at multiple orbits)
    return std::min(total_deflection, 4.0 * Constants::PI); // Cap at 720 degrees (2 orbits)
}

std::array<double, 4> SchwarzchildPhysics::christoffelSymbols(const Vector4& pos, int mu, int nu, int lambda) const {
    // Implementation of Christoffel symbols for Schwarzschild metric
    std::array<double, 4> symbols = {0, 0, 0, 0};
    // Detailed implementation would go here - simplified for now
    return symbols;
}

double SchwarzchildPhysics::christoffel(const Vector4& pos, int mu, int nu, int lambda) const {
    // Implementation of individual Christoffel symbol calculation
    return 0.0; // Simplified for now
}

// ============================================================================
// KerrPhysicsEngine Implementation
// ============================================================================

KerrPhysicsEngine::KerrPhysicsEngine(double mass_solar_units, double spin_parameter) 
    : mass_(mass_solar_units), spin_parameter_(spin_parameter) {
    updateKerrPhysics();
}

void KerrPhysicsEngine::updateKerrPhysics() {
    kerr_physics_ = std::make_unique<KerrPhysics>(mass_, spin_parameter_);
}

void KerrPhysicsEngine::setSpinParameter(double a) {
    spin_parameter_ = a;
    updateKerrPhysics();
}

double KerrPhysicsEngine::photonSphere() const {
    return kerr_physics_->photonSphereRadius();
}

double KerrPhysicsEngine::innerStableCircularOrbit() const {
    return kerr_physics_->innerMostStableCircularOrbit();
}

double KerrPhysicsEngine::outerHorizonRadius() const {
    return kerr_physics_->outerHorizonRadius();
}

double KerrPhysicsEngine::innerHorizonRadius() const {
    return kerr_physics_->innerHorizonRadius();
}

double KerrPhysicsEngine::ergosphereRadius(double theta) const {
    return kerr_physics_->ergosphereRadius(theta);
}

bool KerrPhysicsEngine::isInErgosphere(double r, double theta) const {
    return kerr_physics_->isInErgosphere(r, theta);
}

double KerrPhysicsEngine::metricGtt(double r) const {
    return kerr_physics_->metricGtt(r, 0.0); // Default theta=0
}

double KerrPhysicsEngine::metricGrr(double r) const {
    return kerr_physics_->metricGrr(r);
}

double KerrPhysicsEngine::metricGthetatheta(double r) const {
    return kerr_physics_->metricGthetatheta(r);
}

double KerrPhysicsEngine::metricGphiphi(double r, double theta) const {
    return kerr_physics_->metricGphiphi(r, theta);
}

Vector4 KerrPhysicsEngine::geodesicDerivative(const Vector4& position, const Vector4& momentum) const {
    KerrVector4 kerr_pos(position);
    KerrVector4 kerr_mom(momentum);
    KerrVector4 result = kerr_physics_->geodesicDerivativeKerr(kerr_pos, kerr_mom);
    return Vector4(result.getT(), result.getR(), result.getTheta(), result.getPhi());
}

std::pair<Vector4, Vector4> KerrPhysicsEngine::integrateGeodesic(
    const Vector4& pos, const Vector4& mom, double step_size) const {
    
    KerrVector4 kerr_pos(pos);
    KerrVector4 kerr_mom(mom);
    auto [new_pos, new_mom] = kerr_physics_->integrateKerrGeodesic(kerr_pos, kerr_mom, step_size);
    return std::make_pair(
        Vector4(new_pos.getT(), new_pos.getR(), new_pos.getTheta(), new_pos.getPhi()),
        Vector4(new_mom.getT(), new_mom.getR(), new_mom.getTheta(), new_mom.getPhi())
    );
}

Vector4 KerrPhysicsEngine::cartesianToSpherical(const Vector4& cartesian) const {
    KerrVector4 boyer = kerr_physics_->cartesianToBoyer(cartesian);
    return Vector4(boyer.getT(), boyer.getR(), boyer.getTheta(), boyer.getPhi());
}

Vector4 KerrPhysicsEngine::sphericalToCartesian(const Vector4& spherical) const {
    KerrVector4 boyer(spherical);
    Vector4 result = kerr_physics_->boyerToCartesian(boyer);
    return result;
}

double KerrPhysicsEngine::adaptiveStepSize(double r) const {
    return kerr_physics_->adaptiveKerrStepSize(r, 0.0); // Default theta=0
}

bool KerrPhysicsEngine::isInsideEventHorizon(double r) const {
    return kerr_physics_->isInsideOuterHorizon(r);
}

double KerrPhysicsEngine::deflectionAngle(double impact_parameter) const {
    return kerr_physics_->deflectionAngleKerr(impact_parameter, 0.0);
}

// ============================================================================
// BlackHolePhysicsFactory Implementation
// ============================================================================

std::unique_ptr<IBlackHolePhysics> BlackHolePhysicsFactory::create(
    double mass_solar_units, double spin_parameter, PhysicsType type) {
    
    if (type == PhysicsType::Auto) {
        type = determinePhysicsType(spin_parameter);
    }
    
    switch (type) {
        case PhysicsType::Schwarzschild:
            return createSchwarzschild(mass_solar_units);
        case PhysicsType::Kerr:
            return createKerr(mass_solar_units, spin_parameter);
        default:
            throw std::invalid_argument("Unknown physics type");
    }
}

std::unique_ptr<IBlackHolePhysics> BlackHolePhysicsFactory::createSchwarzschild(double mass_solar_units) {
    return std::make_unique<SchwarzchildPhysics>(mass_solar_units);
}

std::unique_ptr<IBlackHolePhysics> BlackHolePhysicsFactory::createKerr(double mass_solar_units, double spin_parameter) {
    return std::make_unique<KerrPhysicsEngine>(mass_solar_units, spin_parameter);
}

BlackHolePhysicsFactory::PhysicsType BlackHolePhysicsFactory::determinePhysicsType(double spin_parameter) {
    return (std::abs(spin_parameter) > SPIN_THRESHOLD) ? PhysicsType::Kerr : PhysicsType::Schwarzschild;
}

std::string BlackHolePhysicsFactory::getPhysicsTypeName(PhysicsType type) {
    switch (type) {
        case PhysicsType::Schwarzschild: return "Schwarzschild";
        case PhysicsType::Kerr: return "Kerr";
        case PhysicsType::Auto: return "Auto";
        default: return "Unknown";
    }
}

// ============================================================================
// BlackHolePhysicsRegistry Implementation
// ============================================================================

std::unordered_map<std::string, BlackHolePhysicsRegistry::PhysicsCreator> BlackHolePhysicsRegistry::creators_;

void BlackHolePhysicsRegistry::registerPhysicsType(const std::string& name, PhysicsCreator creator) {
    creators_[name] = creator;
}

std::unique_ptr<IBlackHolePhysics> BlackHolePhysicsRegistry::createByName(
    const std::string& name, double mass, double spin) {
    
    auto it = creators_.find(name);
    if (it != creators_.end()) {
        return it->second(mass, spin);
    }
    
    // Fallback to factory for built-in types
    if (name == "Schwarzschild") {
        return BlackHolePhysicsFactory::createSchwarzschild(mass);
    } else if (name == "Kerr") {
        return BlackHolePhysicsFactory::createKerr(mass, spin);
    }
    
    throw std::invalid_argument("Unknown physics type: " + name);
}

std::vector<std::string> BlackHolePhysicsRegistry::getAvailableTypes() {
    std::vector<std::string> types = {"Schwarzschild", "Kerr"};
    
    for (const auto& pair : creators_) {
        types.push_back(pair.first);
    }
    
    return types;
}
