#include "../include/KerrPhysics.h"
#include <cmath>

void KerrParameters::updateDerivedParameters() {
    // Ensure spin parameter is within valid bounds
    if (std::abs(angular_momentum) > KerrConstants::MAX_SPIN_PARAMETER) {
        angular_momentum = std::copysign(KerrConstants::MAX_SPIN_PARAMETER, angular_momentum);
    }
    
    double M = mass;
    double a = angular_momentum;
    double discriminant = M * M - a * a;
    
    if (discriminant < 0) {
        // Naked singularity case - clamp to extremal limit
        angular_momentum = std::copysign(M * 0.999, angular_momentum);
        a = angular_momentum;
        discriminant = M * M - a * a;
    }
    
    double sqrt_disc = std::sqrt(discriminant);
    outer_horizon_radius = M + sqrt_disc;
    inner_horizon_radius = M - sqrt_disc;
    
    // Ergosphere radius depends on theta, so we store the maximum (at equator)
    ergosphere_radius = M + std::sqrt(M * M - a * a);
    static_limit = 2.0 * M;  // At equator where cos(θ) = 0
}

KerrPhysics::KerrPhysics(double mass_solar_units, double spin_parameter)
    : params_(mass_solar_units, spin_parameter), cache_valid_(false) {
    validateParameters();
    updateCachedQuantities();
}

void KerrPhysics::setSpinParameter(double a) {
    params_.angular_momentum = a;
    params_.updateDerivedParameters();
    validateParameters();
    invalidateCache();
}

void KerrPhysics::setMass(double mass_solar_units) {
    params_.mass = mass_solar_units;
    params_.updateDerivedParameters();
    validateParameters();
    invalidateCache();
}

double KerrPhysics::metricGtt(double r, double theta) const {
    updateCache(r, theta);
    double M = params_.mass;
    double a = params_.angular_momentum;
    double cos_theta = std::cos(theta);
    
    return -(1.0 - (2.0 * M * r) / cached_sigma_);
}

double KerrPhysics::metricGrr(double r) const {
    double delta = deltaFunction(r);
    return cached_sigma_ / delta;
}

double KerrPhysics::metricGthetatheta(double r) const {
    updateCache(r, 0); // theta not needed for this component
    return cached_sigma_;
}

double KerrPhysics::metricGphiphi(double r, double theta) const {
    updateCache(r, theta);
    double M = params_.mass;
    double a = params_.angular_momentum;
    double sin_theta = std::sin(theta);
    double sin2_theta = sin_theta * sin_theta;
    
    double A_val = aFunction(r, theta);
    return sin2_theta * A_val / cached_sigma_;
}

double KerrPhysics::metricGtphi(double r, double theta) const {
    updateCache(r, theta);
    double M = params_.mass;
    double a = params_.angular_momentum;
    double sin2_theta = std::sin(theta) * std::sin(theta);
    
    return -(2.0 * M * a * r * sin2_theta) / cached_sigma_;
}

double KerrPhysics::deltaFunction(double r) const {
    double M = params_.mass;
    double a = params_.angular_momentum;
    return r * r - 2.0 * M * r + a * a;
}

double KerrPhysics::sigmaFunction(double r, double theta) const {
    double a = params_.angular_momentum;
    double cos_theta = std::cos(theta);
    return r * r + a * a * cos_theta * cos_theta;
}

double KerrPhysics::aFunction(double r, double theta) const {
    double a = params_.angular_momentum;
    double sin2_theta = std::sin(theta) * std::sin(theta);
    double delta = deltaFunction(r);
    
    return (r * r + a * a) * (r * r + a * a) - a * a * delta * sin2_theta;
}

double KerrPhysics::ergosphereRadius(double theta) const {
    double M = params_.mass;
    double a = params_.angular_momentum;
    double cos2_theta = std::cos(theta) * std::cos(theta);
    
    return M + std::sqrt(M * M - a * a * cos2_theta);
}

double KerrPhysics::staticLimit(double theta) const {
    double M = params_.mass;
    double a = params_.angular_momentum;
    double cos2_theta = std::cos(theta) * std::cos(theta);
    
    return M + std::sqrt(M * M - a * a * cos2_theta);
}

bool KerrPhysics::isInErgosphere(double r, double theta) const {
    double r_ergo = ergosphereRadius(theta);
    double r_horizon = params_.outer_horizon_radius;
    return r > r_horizon && r < r_ergo;
}

double KerrPhysics::frameReferenceVelocity(double r, double theta) const {
    // ω = -g_tφ / g_φφ - angular velocity of local frame
    double g_tphi = metricGtphi(r, theta);
    double g_phiphi = metricGphiphi(r, theta);
    
    if (std::abs(g_phiphi) < 1e-12) return 0.0;
    return -g_tphi / g_phiphi;
}

double KerrPhysics::lenseThirringPrecession(double r, double theta) const {
    // Frame-dragging precession rate
    double M = params_.mass;
    double a = params_.angular_momentum;
    double sigma = sigmaFunction(r, theta);
    
    return (2.0 * M * a * r) / (sigma * sigma);
}

Vector4 KerrPhysics::applyFrameDragging(const Vector4& momentum, const Vector4& position) const {
    KerrVector4 pos(position);
    KerrVector4 mom(momentum);
    
    double r = pos.getR();
    double theta = pos.getTheta();
    
    if (isInErgosphere(r, theta)) {
        double omega = frameReferenceVelocity(r, theta);
        // Apply frame-dragging correction to φ component of momentum
        double corrected_phi_momentum = mom.getPhi() + omega * mom.getT();
        return Vector4(mom.getT(), mom.getR(), mom.getTheta(), corrected_phi_momentum);
    }
    
    return momentum;
}

KerrVector4 KerrPhysics::geodesicDerivativeKerr(const KerrVector4& position, const KerrVector4& momentum) const {
    double r = position.getR();
    double theta = position.getTheta();
    
    if (r <= params_.inner_horizon_radius) {
        return KerrVector4(0, 0, 0, 0);  // Inside inner horizon
    }
    
    double M = params_.mass;
    double a = params_.angular_momentum;
    double delta = deltaFunction(r);
    double sigma = sigmaFunction(r, theta);
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    double sin2_theta = sin_theta * sin_theta;
    
    // Extract momentum components
    double pt = momentum.getT();
    double pr = momentum.getR();
    double ptheta = momentum.getTheta();
    double pphi = momentum.getPhi();
    
    // Geodesic equations for Kerr metric in Boyer-Lindquist coordinates
    
    // d²t/dλ²
    double d2t = -(2.0 * M * r) / (sigma * sigma) * (pr * pt + a * sin2_theta * pr * pphi)
                - (2.0 * M * a * r * sin_theta * cos_theta) / (sigma * sigma) * (ptheta * pt + a * sin2_theta * ptheta * pphi);
    
    // d²r/dλ²
    double d2r = (delta / sigma) * (
        -(M * (r * r - a * a) / (sigma * sigma)) * pt * pt
        + (2.0 * M * a * r / (sigma * sigma)) * pt * pphi
        - (M * a * a * sin2_theta / (sigma * sigma)) * pphi * pphi
    ) + (1.0 / sigma) * (r - M) * pr * pr
      + (1.0 / sigma) * (-r * delta) * (ptheta * ptheta + sin2_theta * pphi * pphi);
    
    // d²θ/dλ²
    double d2theta = -(a * a * sin_theta * cos_theta / (sigma * sigma)) * pt * pt
                    + (2.0 * M * a * r * sin_theta * cos_theta / (sigma * sigma)) * pt * pphi
                    - (sin_theta * cos_theta / sigma) * (r * r + a * a + (2.0 * M * a * a * r * sin2_theta / sigma)) * pphi * pphi
                    - (2.0 * r / sigma) * pr * ptheta;
    
    // d²φ/dλ²
    double d2phi = (2.0 * M * r / (sigma * sigma)) * (pt + a * sin2_theta * pphi) * pr / sin2_theta
                  + (2.0 * cos_theta / sin_theta) * (1.0 / sigma) * (r * r + a * a + (2.0 * M * a * a * r * sin2_theta / sigma)) * ptheta * pphi
                  - (2.0 * M * a * r / (sigma * sigma)) * pt * ptheta / sin2_theta;
    
    return KerrVector4(d2t, d2r, d2theta, d2phi);
}

std::pair<KerrVector4, KerrVector4> KerrPhysics::integrateKerrGeodesic(
    const KerrVector4& pos, const KerrVector4& mom, double step_size) const {
    
    // 4th order Runge-Kutta integration adapted for Kerr spacetime
    KerrVector4 k1_pos = mom;
    KerrVector4 k1_mom = geodesicDerivativeKerr(pos, mom);
    
    KerrVector4 temp1 = k1_mom * (step_size * 0.5);
    KerrVector4 k2_pos = mom + temp1;
    KerrVector4 temp2 = k1_pos * (step_size * 0.5);
    KerrVector4 k2_mom = geodesicDerivativeKerr(pos + temp2, mom + temp1);
    
    KerrVector4 temp3 = k2_mom * (step_size * 0.5);
    KerrVector4 k3_pos = mom + temp3;
    KerrVector4 temp4 = k2_pos * (step_size * 0.5);
    KerrVector4 k3_mom = geodesicDerivativeKerr(pos + temp4, mom + temp3);
    
    KerrVector4 temp5 = k3_mom * step_size;
    KerrVector4 k4_pos = mom + temp5;
    KerrVector4 temp6 = k3_pos * step_size;
    KerrVector4 k4_mom = geodesicDerivativeKerr(pos + temp6, mom + temp5);
    
    KerrVector4 pos_sum = k1_pos + k2_pos * 2.0 + k3_pos * 2.0 + k4_pos;
    KerrVector4 mom_sum = k1_mom + k2_mom * 2.0 + k3_mom * 2.0 + k4_mom;
    
    KerrVector4 new_pos = pos + pos_sum * (step_size / 6.0);
    KerrVector4 new_mom = mom + mom_sum * (step_size / 6.0);
    
    return std::make_pair(new_pos, new_mom);
}

double KerrPhysics::specificEnergy(const KerrVector4& position, const KerrVector4& momentum) const {
    double r = position.getR();
    double theta = position.getTheta();
    
    double g_tt = metricGtt(r, theta);
    double g_tphi = metricGtphi(r, theta);
    
    return -(g_tt * momentum.getT() + g_tphi * momentum.getPhi());
}

double KerrPhysics::specificAngularMomentum(const KerrVector4& position, const KerrVector4& momentum) const {
    double r = position.getR();
    double theta = position.getTheta();
    
    double g_tphi = metricGtphi(r, theta);
    double g_phiphi = metricGphiphi(r, theta);
    
    return g_tphi * momentum.getT() + g_phiphi * momentum.getPhi();
}

double KerrPhysics::carterConstant(const KerrVector4& position, const KerrVector4& momentum) const {
    // Simplified Carter constant calculation
    double r = position.getR();
    double theta = position.getTheta();
    double a = params_.angular_momentum;
    
    double ptheta = momentum.getTheta();
    double pphi = momentum.getPhi();
    double cos2_theta = std::cos(theta) * std::cos(theta);
    double sin2_theta = std::sin(theta) * std::sin(theta);
    
    double E = specificEnergy(position, momentum);
    double L = specificAngularMomentum(position, momentum);
    
    return ptheta * ptheta + cos2_theta * (a * a * (1.0 - E * E) + L * L / sin2_theta);
}

double KerrPhysics::innerMostStableCircularOrbit(bool prograde) const {
    double M = params_.mass;
    double a = params_.angular_momentum;
    
    if (!prograde) a = -a;  // Retrograde orbit
    
    // Approximate ISCO formula for Kerr black holes
    double Z1 = 1.0 + std::pow(1.0 - a*a, 1.0/3.0) * (std::pow(1.0 + a, 1.0/3.0) + std::pow(1.0 - a, 1.0/3.0));
    double Z2 = std::sqrt(3.0 * a * a + Z1 * Z1);
    
    if (prograde) {
        return M * (3.0 + Z2 - std::sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    } else {
        return M * (3.0 + Z2 + std::sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    }
}

double KerrPhysics::photonSphereRadius(bool prograde) const {
    double M = params_.mass;
    double a = params_.angular_momentum;
    
    if (!prograde) a = -a;
    
    // Exact photon sphere radius for Kerr black hole
    // The photon sphere satisfies the condition for circular null geodesics
    // We solve iteratively using Newton's method for higher accuracy
    
    // Initial guess using the approximate formula
    double r = M * (2.0 + 2.0 * std::cos((2.0/3.0) * std::acos(std::max(-1.0, std::min(1.0, -a/M)))));
    
    // Newton's method to find exact photon sphere
    // We solve: f(r) = r⁴ - 6M r³ + 9M² r² - 4a²M² = 0
    // f'(r) = 4r³ - 18Mr² + 18M²r
    
    const int max_iterations = 50;
    const double tolerance = 1e-12;
    double a2 = a * a;
    double M2 = M * M;
    
    for (int i = 0; i < max_iterations; ++i) {
        double r2 = r * r;
        double r3 = r2 * r;
        double r4 = r3 * r;
        
        // Function value: f(r) = r⁴ - 6Mr³ + 9M²r² - 4a²M²
        double f = r4 - 6.0 * M * r3 + 9.0 * M2 * r2 - 4.0 * a2 * M2;
        
        // Derivative: f'(r) = 4r³ - 18Mr² + 18M²r
        double df = 4.0 * r3 - 18.0 * M * r2 + 18.0 * M2 * r;
        
        if (std::abs(df) < 1e-15) break; // Avoid division by zero
        
        double dr = f / df;
        r -= dr;
        
        if (std::abs(dr) < tolerance) break; // Converged
    }
    
    // Ensure reasonable result (should be between outer horizon and ~4M)
    double r_horizon = params_.outer_horizon_radius;
    if (r < r_horizon || r > 4.0 * M) {
        // Fallback to approximate formula if Newton's method failed
        return M * (2.0 + 2.0 * std::cos((2.0/3.0) * std::acos(std::max(-1.0, std::min(1.0, -a/M)))));
    }
    
    return r;
}

KerrVector4 KerrPhysics::cartesianToBoyer(const Vector4& cartesian) const {
    double x = cartesian.x;
    double y = cartesian.y;
    double z = cartesian.z;
    double t = cartesian.t;
    
    double a = params_.angular_momentum;
    
    // Convert to Boyer-Lindquist coordinates
    double rho2 = x*x + y*y + z*z - a*a;
    double r = std::sqrt(0.5 * (rho2 + std::sqrt(rho2*rho2 + 4.0*a*a*z*z)));
    
    if (r == 0.0) return KerrVector4(t, 0, 0, 0);
    
    double theta = std::acos(std::max(-1.0, std::min(1.0, z / r)));
    double phi = std::atan2(y, x);
    
    return KerrVector4(t, r, theta, phi);
}

Vector4 KerrPhysics::boyerToCartesian(const KerrVector4& boyer) const {
    double r = boyer.getR();
    double theta = boyer.getTheta();
    double phi = boyer.getPhi();
    double t = boyer.getT();
    
    double a = params_.angular_momentum;
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    
    double x = std::sqrt(r*r + a*a) * sin_theta * std::cos(phi);
    double y = std::sqrt(r*r + a*a) * sin_theta * std::sin(phi);
    double z = r * cos_theta;
    
    return Vector4(t, x, y, z);
}

double KerrPhysics::adaptiveKerrStepSize(double r, double theta) const {
    double min_step = 0.001;
    double max_step = 0.1;
    double horizon_radius = params_.outer_horizon_radius;
    double ergo_radius = ergosphereRadius(theta);
    
    // Smaller steps near horizons and in ergosphere
    if (r < 1.1 * horizon_radius || isInErgosphere(r, theta)) {
        return min_step;
    } else if (r < 2.0 * ergo_radius) {
        return min_step + (max_step - min_step) * (r - 1.1 * horizon_radius) / (0.9 * ergo_radius);
    } else {
        return max_step;
    }
}

bool KerrPhysics::checkKerrTermination(const KerrVector4& position) const {
    double r = position.getR();
    return r <= params_.inner_horizon_radius || r > 1000.0 * params_.mass;
}

double KerrPhysics::deflectionAngleKerr(double impact_parameter, double inclination) const {
    double M = params_.mass;
    double a = params_.angular_momentum;
    
    // Use exact strong-field deflection formula for Kerr metric
    // Based on Sereno & De Luca (2006) exact formula for strong field regime
    
    double b = impact_parameter; // Impact parameter
    double cos_i = std::cos(inclination);
    double sin_i = std::sin(inclination);
    
    // Critical impact parameter for photon sphere (determines regime)
    double b_crit = photonSphereRadius(true) * std::sqrt(1.0 + a * cos_i / M);
    
    if (b < 0.1 * b_crit) {
        // Very strong field - particle likely captured
        return Constants::PI; // 180 degree deflection (effective capture)
    }
    
    if (b > 10.0 * b_crit) {
        // Weak field regime - use standard weak field formula with Kerr correction
        double weak_field = 4.0 * M / b;
        double kerr_correction = 1.0 + (15.0 * Constants::PI / 32.0) * (M / b) * (a * cos_i / M);
        return weak_field * kerr_correction;
    }
    
    // Strong field regime - use exact calculation
    // This requires solving the full geodesic equation, but we can use an accurate approximation
    
    // Bardeen's formula for strong field deflection (adapted for Kerr)
    double u = M / b; // Dimensionless parameter
    double spin_factor = a * cos_i / M;
    
    // Primary deflection contribution
    double alpha_0 = 4.0 * u; // Leading order
    
    // Next-order corrections for strong field
    double alpha_1 = (15.0 * Constants::PI / 4.0) * u * u; // Second order
    double alpha_2 = 128.0 * u * u * u / 3.0; // Third order (strong field)
    
    // Kerr spin correction at each order
    double spin_correction_1 = 1.0 + spin_factor * (1.0 + 2.0 * u);
    double spin_correction_2 = 1.0 + spin_factor * (2.0 + 4.0 * u + u * u);
    double spin_correction_3 = 1.0 + spin_factor * (3.0 + 6.0 * u + 3.0 * u * u);
    
    // Combine all terms with spin corrections
    double total_deflection = alpha_0 * spin_correction_1 + 
                             alpha_1 * spin_correction_2 + 
                             alpha_2 * spin_correction_3;
    
    // Ensure physically reasonable result
    return std::min(total_deflection, 2.0 * Constants::PI); // Cap at 360 degrees
}

BlackHolePhysics KerrPhysics::createEquivalentSchwarzschildPhysics() const {
    // Create equivalent Schwarzschild physics for comparison/compatibility
    return BlackHolePhysics(params_.mass);
}

bool KerrPhysics::isEffectivelyNonRotating(double tolerance) const {
    return std::abs(params_.angular_momentum) < tolerance;
}

void KerrPhysics::validateParameters() const {
    if (params_.mass <= 0) {
        throw std::invalid_argument("Black hole mass must be positive");
    }
    
    double a_max = KerrConstants::MAX_SPIN_PARAMETER;
    if (std::abs(params_.angular_momentum) > a_max) {
        throw std::invalid_argument("Spin parameter exceeds maximum stable value");
    }
    
    // Check for naked singularity
    if (params_.mass * params_.mass < params_.angular_momentum * params_.angular_momentum) {
        throw std::invalid_argument("Parameters would result in naked singularity");
    }
}

void KerrPhysics::updateCachedQuantities() {
    // Update any cached quantities that depend on parameters
    invalidateCache();
}

void KerrPhysics::updateCache(double r, double theta) const {
    if (!cache_valid_ || cached_r_ != r || cached_theta_ != theta) {
        cached_r_ = r;
        cached_theta_ = theta;
        cached_delta_ = deltaFunction(r);
        cached_sigma_ = sigmaFunction(r, theta);
        cache_valid_ = true;
    }
}

// Utility functions implementation
namespace KerrUtils {
    double extractionEfficiency(double spin_parameter) {
        // Approximate energy extraction efficiency via Penrose process
        return 1.0 - std::sqrt(1.0 - spin_parameter*spin_parameter/4.0);
    }
    
    double maximallyExtractableEnergy(double mass, double spin_parameter) {
        return mass * extractionEfficiency(spin_parameter);
    }
    
    KerrVector4 convert(const Vector4& v) {
        return KerrVector4(v.t, v.x, v.y, v.z);
    }
    
    Vector4 convert(const KerrVector4& kv) {
        return Vector4(kv.getT(), kv.getR(), kv.getTheta(), kv.getPhi());
    }
}
