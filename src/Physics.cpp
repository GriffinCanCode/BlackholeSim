#include "Physics.h"
#include <cmath>
#include <algorithm>

BlackHolePhysics::BlackHolePhysics(double mass_solar_units) 
    : mass_(mass_solar_units) {
    // Calculate Schwarzschild radius in geometric units
    schwarzschild_radius_ = 2.0 * mass_solar_units;
}

double BlackHolePhysics::metricGtt(double r) const {
    return -(1.0 - schwarzschild_radius_ / r);
}

double BlackHolePhysics::metricGrr(double r) const {
    return 1.0 / (1.0 - schwarzschild_radius_ / r);
}

double BlackHolePhysics::metricGthetatheta(double r) const {
    return r * r;
}

double BlackHolePhysics::metricGphiphi(double r, double theta) const {
    return r * r * std::sin(theta) * std::sin(theta);
}

Vector4 BlackHolePhysics::geodesicDerivative(const Vector4& position, const Vector4& momentum) const {
    double r = position.x;  // Using x as radial coordinate in spherical
    double theta = position.y;  // Using y as theta coordinate
    double phi = position.z;    // Using z as phi coordinate
    
    if (r <= schwarzschild_radius_) {
        return Vector4(0, 0, 0, 0);  // Inside event horizon
    }
    
    double rs = schwarzschild_radius_;
    double r2 = r * r;
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

std::pair<Vector4, Vector4> BlackHolePhysics::integrateGeodesic(
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

Vector4 BlackHolePhysics::cartesianToSpherical(const Vector4& cartesian) const {
    double x = cartesian.x;
    double y = cartesian.y;
    double z = cartesian.z;
    
    double r = std::sqrt(x*x + y*y + z*z);
    if (r == 0.0) return Vector4(cartesian.t, 0, 0, 0);
    
    double theta = std::acos(std::max(-1.0, std::min(1.0, z / r)));
    double phi = std::atan2(y, x);
    
    return Vector4(cartesian.t, r, theta, phi);
}

Vector4 BlackHolePhysics::sphericalToCartesian(const Vector4& spherical) const {
    double r = spherical.x;
    double theta = spherical.y;
    double phi = spherical.z;
    
    double x = r * std::sin(theta) * std::cos(phi);
    double y = r * std::sin(theta) * std::sin(phi);
    double z = r * std::cos(theta);
    
    return Vector4(spherical.t, x, y, z);
}

double BlackHolePhysics::adaptiveStepSize(double r) const {
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

double BlackHolePhysics::deflectionAngle(double impact_parameter) const {
    // Approximate deflection angle for weak field limit
    // More accurate calculation would require full geodesic integration
    return 4.0 * schwarzschild_radius_ / impact_parameter;
}

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
