#include "RelativisticJets.h"
#include <cmath>

RelativisticJets::RelativisticJets(const BlackHolePhysics& physics)
    : physics_(physics), 
      jet_model_(Vector4(0, 0, 0, 1), 0.1),  // Default vertical jet
      spectral_resolution_(64), variability_enabled_(true), variability_timescale_(100.0) {
    
    // Set reasonable default jet parameters
    jet_model_.bulk_lorentz_factor = 10.0;     // Γ = 10 (highly relativistic)
    jet_model_.opening_angle = 0.1;            // ~6 degrees
    jet_model_.magnetic_field = 1e4;           // 10 kGauss
    jet_model_.particle_density = 1e12;        // 10^12 particles/m³
}

double RelativisticJets::calculateBeamingFactor(const Vector4& jet_velocity, const Vector4& observer_direction) const {
    // Relativistic beaming factor: δ = 1/(γ(1 - β·n̂))
    // This gives the factor by which emission is enhanced/diminished due to relativistic motion
    
    double gamma = lorentzFactor(jet_velocity);
    if (gamma <= 1.0) return 1.0;  // Non-relativistic case
    
    Vector4 normalized_velocity = VectorMath::normalize3D(jet_velocity);
    Vector4 normalized_observer = VectorMath::normalize3D(observer_direction);
    
    double beta_magnitude = VectorMath::magnitude3D(jet_velocity) / SpectralConstants::SPEED_OF_LIGHT;
    double beta_dot_n = VectorMath::dot3D(normalized_velocity, normalized_observer);
    
    // Beaming factor (Doppler factor)
    double beaming_factor = 1.0 / (gamma * (1.0 - beta_magnitude * beta_dot_n));
    
    // Clamp to reasonable values to avoid numerical issues
    return std::max(1e-6, std::min(1e6, beaming_factor));
}

double RelativisticJets::dopplerBoostingFactor(const Vector4& velocity, const Vector4& photon_direction) const {
    // Calculate the Doppler boosting factor for photon frequency/intensity
    // ν_observed = δ * ν_emitted, where δ is the Doppler factor
    
    double gamma = lorentzFactor(velocity);
    double beta = VectorMath::magnitude3D(velocity) / SpectralConstants::SPEED_OF_LIGHT;
    
    Vector4 norm_vel = VectorMath::normalize3D(velocity);
    Vector4 norm_photon = VectorMath::normalize3D(photon_direction);
    
    double cos_theta = VectorMath::dot3D(norm_vel, norm_photon);
    
    // Relativistic Doppler formula
    return gamma * (1.0 + beta * cos_theta);
}

Vector4 RelativisticJets::calculateJetVelocity(const Vector4& position, const RelativisticJetModel& model) const {
    // Calculate bulk jet velocity at given position
    // Assume jet accelerates from base to terminal Lorentz factor
    
    double horizon_radius = physics_.outerHorizonRadius();
    double r = VectorMath::magnitude3D(position);
    
    if (r < 2.0 * horizon_radius) return Vector4(0, 0, 0, 0);  // Too close to black hole
    
    // Height along jet axis
    Vector4 jet_direction = VectorMath::normalize3D(model.jet_axis);
    double height = VectorMath::dot3D(position, jet_direction);
    
    // Acceleration profile - jet accelerates from launch point
    double acceleration_length = 5.0 * horizon_radius;
    double velocity_fraction = std::min(1.0, height / acceleration_length);
    
    // Terminal velocity magnitude
    double beta_terminal = std::sqrt(1.0 - 1.0/(model.bulk_lorentz_factor * model.bulk_lorentz_factor));
    double velocity_magnitude = beta_terminal * velocity_fraction * SpectralConstants::SPEED_OF_LIGHT;
    
    // Base jet velocity (purely radial initially)
    Vector4 base_velocity = jet_direction * velocity_magnitude;
    
    // Add frame-dragging effects for rotating black holes
    if (physics_.isRotating()) {
        // Convert to spherical coordinates to get theta
        Vector4 spherical_pos = physics_.cartesianToSpherical(position);
        double theta = spherical_pos.y;
        
        // Check if position is in or near ergosphere
        if (physics_.isInErgosphere(r, theta)) {
            // Frame-dragging velocity component (azimuthal)
            double frame_drag_omega = 2.0 * physics_.getSpinParameter() * horizon_radius / 
                                     (r * r + physics_.getSpinParameter() * physics_.getSpinParameter() * std::cos(theta) * std::cos(theta));
            
            // Add rotational component to jet velocity
            // φ direction in spherical coordinates corresponds to rotation
            Vector4 phi_direction(0, 0, 0, 1);  // Simplified: assume φ-hat direction
            Vector4 frame_drag_velocity = phi_direction * (frame_drag_omega * r * std::sin(theta) * 0.3); // Scale factor
            
            base_velocity = base_velocity + frame_drag_velocity;
        }
        
        // Additional twist from frame-dragging even outside ergosphere
        if (r < 10.0 * horizon_radius) {
            double twist_strength = physics_.getSpinParameter() * horizon_radius * horizon_radius / (r * r);
            Vector4 twist_velocity = Vector4(0, -base_velocity.z * twist_strength, base_velocity.y * twist_strength, 0);
            base_velocity = base_velocity + twist_velocity;
        }
    }
    
    return base_velocity;
}

SpectralColor RelativisticJets::synchrotronEmission(const JetElement& element, double frequency) const {
    // Calculate synchrotron emission from relativistic electrons in magnetic field
    
    if (element.magnetic_field <= 0 || element.density <= 0) {
        return SpectralColor(0, 0, 0, 0);
    }
    
    // Characteristic synchrotron frequency
    double B = element.magnetic_field;
    double gamma_min = jet_model_.min_electron_energy;
    double freq_sync = synchrotronFrequency(gamma_min, B);
    
    // Power-law electron energy distribution
    double p = jet_model_.power_law_index;
    
    // Synchrotron emissivity - simplified but physically motivated
    double emissivity = 0.0;
    
    if (frequency > freq_sync) {
        // Above characteristic frequency - power law behavior
        emissivity = element.density * pow(frequency / freq_sync, -(p-1)/2.0);
    } else {
        // Below characteristic frequency - self-absorbed regime
        emissivity = element.density * pow(frequency / freq_sync, 5.0/2.0);
    }
    
    // Convert to spectral color - approximate mapping
    double wavelength = SpectralConstants::SPEED_OF_LIGHT / frequency * 1e9;  // Convert to nm
    SpectralColor spectral_color;
    
    if (wavelength >= 400 && wavelength <= 700) {
        // Visible light
        // Simple temperature to RGB conversion (placeholder)
        if (element.temperature > 10000) {
            spectral_color = SpectralColor(0.8, 0.8, 1.0);  // Blue-white
        } else if (element.temperature > 6000) {
            spectral_color = SpectralColor(1.0, 1.0, 1.0);  // White
        } else {
            spectral_color = SpectralColor(1.0, 0.6, 0.2);  // Orange-red
        }
    } else {
        // Non-visible - map to intensity with red shift for radio, blue for X-ray
        if (wavelength > 700) {
            spectral_color = SpectralColor(1.0, 0.5, 0.2);  // Red for radio
        } else {
            spectral_color = SpectralColor(0.2, 0.5, 1.0);  // Blue for high energy
        }
    }
    
    return spectral_color * emissivity * 1e-20;  // Scale factor for reasonable intensities
}

SpectralColor RelativisticJets::inverseComptonEmission(const JetElement& element, double photon_energy) const {
    // Inverse Compton scattering of soft photons by relativistic electrons
    
    if (element.density <= 0) return SpectralColor(0, 0, 0, 0);
    
    // Thomson scattering cross-section
    double sigma_T = THOMSON_CROSS_SECTION;
    
    // Energy of scattered photon (simplified Klein-Nishina formula)
    double gamma_e = jet_model_.min_electron_energy;
    double scattered_energy = photon_energy * gamma_e * gamma_e;
    
    // Inverse Compton emissivity
    double ic_emissivity = element.density * sigma_T * scattered_energy / (photon_energy * photon_energy);
    
    // Map energy to color
    double wavelength = SpectralConstants::PLANCK_CONSTANT * SpectralConstants::SPEED_OF_LIGHT / 
                       (scattered_energy * 1.602176634e-19) * 1e9;  // Convert to nm
    
    // Simple temperature to RGB conversion for harder spectrum (placeholder)
    double temp = element.temperature * 2.0;
    SpectralColor spectral_color;
    if (temp > 15000) {
        spectral_color = SpectralColor(0.6, 0.6, 1.0);  // Blue
    } else if (temp > 8000) {
        spectral_color = SpectralColor(0.9, 0.9, 1.0);  // Blue-white
    } else {
        spectral_color = SpectralColor(1.0, 0.8, 0.4);  // Yellow-orange
    }
    
    return spectral_color * ic_emissivity * 1e-15;
}

double RelativisticJets::synchrotronFrequency(double electron_energy, double magnetic_field) const {
    // Characteristic synchrotron frequency: ν = (3/2) * (e*B*γ²)/(2π*me*c)
    
    double gamma = electron_energy;
    double freq = (3.0/2.0) * (ELECTRON_CHARGE * magnetic_field * gamma * gamma) / 
                  (2.0 * Constants::PI * ELECTRON_MASS * SpectralConstants::SPEED_OF_LIGHT);
    
    return freq;
}

std::vector<JetElement> RelativisticJets::generateJetElements(const RelativisticJetModel& model, int resolution) const {
    std::vector<JetElement> elements;
    elements.reserve(resolution);
    
    double horizon_radius = physics_.outerHorizonRadius();
    
    // Generate elements along jet axis and perpendicular directions
    for (int i = 0; i < resolution; ++i) {
        double height_fraction = double(i) / (resolution - 1);
        double height = 2.0 * horizon_radius + height_fraction * model.jet_height * horizon_radius;
        
        // Position along jet axis
        Vector4 position = model.jet_axis * height;
        
        JetElement element(position);
        
        // Calculate bulk velocity
        element.velocity = calculateJetVelocity(position, model);
        
        // Density profile - decreases with height and radius
        element.density = jetDensityProfile(position, model);
        
        // Magnetic field - may be ordered or turbulent, enhanced by frame-dragging near rotating BH
        double field_strength = model.magnetic_field * std::pow(horizon_radius/height, 1.5);  // Decreases with distance
        
        // Frame-dragging can enhance magnetic field strength near ergosphere
        if (physics_.isRotating()) {
            Vector4 spherical_pos = physics_.cartesianToSpherical(position);
            double r = spherical_pos.x;
            double theta = spherical_pos.y;
            
            if (physics_.isInErgosphere(r, theta)) {
                // Field amplification due to frame-dragging shear
                double amplification = 1.0 + 0.5 * physics_.getSpinParameter() * horizon_radius / r;
                field_strength *= amplification;
            }
        }
        
        element.magnetic_field = field_strength;
        
        // Temperature from equipartition
        element.temperature = 1e9 * std::pow(element.magnetic_field / model.magnetic_field, 2.0/3.0);
        
        // Synchrotron power
        element.synchrotron_power = calculateSynchrotronPower(model.min_electron_energy, element.magnetic_field);
        
        // Beaming factor toward observer (simplified - assuming observer at infinity in z direction)
        Vector4 observer_direction(0, 0, 0, 1);
        element.beaming_factor = calculateBeamingFactor(element.velocity, observer_direction);
        
        elements.push_back(element);
    }
    
    return elements;
}

bool RelativisticJets::isPointInJet(const Vector4& position, const RelativisticJetModel& model) const {
    // Check if a point is within the jet cone
    
    double height = VectorMath::dot3D(position, model.jet_axis);
    if (height < 0 || height > model.jet_height * physics_.outerHorizonRadius()) {
        return false;
    }
    
    // Distance from jet axis  
    Vector4 axis_component = model.jet_axis * height;
    Vector4 perpendicular = Vector4(0, position.x - axis_component.x, position.y - axis_component.y, position.z - axis_component.z);
    double perp_distance = VectorMath::magnitude3D(perpendicular);
    
    // Jet radius at this height (conical expansion)
    double local_radius = model.base_radius * physics_.outerHorizonRadius() * 
                         (1.0 + height * std::tan(model.opening_angle) / (model.base_radius * physics_.outerHorizonRadius()));
    
    return perp_distance <= local_radius;
}

double RelativisticJets::jetDensityProfile(const Vector4& position, const RelativisticJetModel& model) const {
    if (!isPointInJet(position, model)) return 0.0;
    
    double height = VectorMath::dot3D(position, model.jet_axis);
    double horizon_radius = physics_.outerHorizonRadius();
    
    // Density decreases with height and radius
    Vector4 axis_component = model.jet_axis * height;
    Vector4 perpendicular = Vector4(0, position.x - axis_component.x, position.y - axis_component.y, position.z - axis_component.z);
    double r_perp = VectorMath::magnitude3D(perpendicular);
    double local_radius = model.base_radius * horizon_radius * (1.0 + height * std::tan(model.opening_angle) / (model.base_radius * horizon_radius));
    
    // Gaussian profile across jet width
    double radial_profile = std::exp(-(r_perp * r_perp) / (0.5 * local_radius * local_radius));
    
    // Power-law decrease along jet
    double height_profile = std::pow(2.0*horizon_radius/std::max(height, 2.0*horizon_radius), 2.0);
    
    return model.particle_density * radial_profile * height_profile;
}

double RelativisticJets::beamingConeAngle(double lorentz_factor) const {
    // Half-opening angle of relativistic beaming cone: θ ≈ 1/γ
    if (lorentz_factor <= 1.0) return Constants::PI;  // Non-relativistic - isotropic
    
    return 1.0 / lorentz_factor;
}

SpectralColor RelativisticJets::calculateBeamedEmission(const JetElement& element, const Vector4& observer_position) const {
    // Calculate emission as seen by observer, including relativistic beaming
    
    Vector4 to_observer = Vector4(0, observer_position.x - element.position.x, observer_position.y - element.position.y, observer_position.z - element.position.z);
    Vector4 observer_direction = VectorMath::normalize3D(to_observer);
    
    // Beaming factor
    double beaming = calculateBeamingFactor(element.velocity, observer_direction);
    
    // Angle between jet velocity and observer direction
    Vector4 norm_velocity = VectorMath::normalize3D(element.velocity);
    double cos_viewing_angle = VectorMath::dot3D(norm_velocity, observer_direction);
    double viewing_angle = std::acos(std::max(-1.0, std::min(1.0, cos_viewing_angle)));
    
    // Check if within beaming cone
    double beaming_cone = beamingConeAngle(lorentzFactor(element.velocity));
    
    SpectralColor intrinsic_emission = element.emission;
    
    if (viewing_angle > 2.0 * beaming_cone) {
        // Outside beaming cone - heavily suppressed
        return intrinsic_emission * (beaming * 1e-6);
    } else if (viewing_angle < beaming_cone) {
        // Within main beaming cone - strong enhancement
        // Frequency boosting: ν_obs = δ * ν_emit
        // Intensity boosting: I_obs = δ^3 * I_emit (for moving source)
        return intrinsic_emission * std::pow(beaming, 3.0);
    } else {
        // Transition region
        double suppression = std::exp(-(viewing_angle - beaming_cone) / beaming_cone);
        return intrinsic_emission * beaming * suppression;
    }
}

double RelativisticJets::lightAberrationAngle(const Vector4& velocity, const Vector4& emission_direction) const {
    // Calculate aberration angle due to relativistic motion
    
    double gamma = lorentzFactor(velocity);
    double beta = VectorMath::magnitude3D(velocity) / SpectralConstants::SPEED_OF_LIGHT;
    
    Vector4 norm_vel = VectorMath::normalize3D(velocity);
    Vector4 norm_emission = VectorMath::normalize3D(emission_direction);
    
    double cos_theta_emit = VectorMath::dot3D(norm_vel, norm_emission);
    
    // Relativistic aberration formula
    double cos_theta_obs = (cos_theta_emit + beta) / (1.0 + beta * cos_theta_emit);
    
    return std::acos(std::max(-1.0, std::min(1.0, cos_theta_obs)));
}

SpectralColor RelativisticJets::sampleJetEmission(const Vector4& ray_position, const Vector4& ray_direction,
                                                 const Vector4& observer_position, const RelativisticJetModel& model) const {
    // Sample jet emission along a ray path
    
    if (!model.active) return SpectralColor(0, 0, 0, 0);
    
    SpectralColor total_emission(0, 0, 0, 0);
    
    // Generate jet elements for sampling
    auto elements = generateJetElements(model, 50);
    
    for (const auto& element : elements) {
        // Check if ray passes close to this jet element
        Vector4 to_element = Vector4(0, element.position.x - ray_position.x, element.position.y - ray_position.y, element.position.z - ray_position.z);
        double distance_to_ray = VectorMath::magnitude3D(VectorMath::cross3D(to_element, ray_direction)) / 
                                VectorMath::magnitude3D(ray_direction);
        
        double element_size = model.base_radius * physics_.outerHorizonRadius();
        if (distance_to_ray > element_size) continue;
        
        // Calculate emission with relativistic beaming
        SpectralColor beamed_emission = calculateBeamedEmission(element, observer_position);
        
        // Add synchrotron emission at multiple frequencies
        double base_frequency = 1e15;  // ~1 PHz (near-infrared)
        for (int freq_idx = 0; freq_idx < 3; ++freq_idx) {
            double frequency = base_frequency * std::pow(10.0, freq_idx - 1);
            SpectralColor sync_emission = synchrotronEmission(element, frequency);
            beamed_emission = beamed_emission + sync_emission;
        }
        
        // Weight by distance and element density
        double weight = element.density / (1.0 + distance_to_ray * distance_to_ray);
        total_emission = total_emission + beamed_emission * weight * 1e-10;
    }
    
    return total_emission;
}

void RelativisticJets::updateJetElements(std::vector<JetElement>& elements, double time_step) const {
    for (auto& element : elements) {
        // Apply variability if enabled
        if (variability_enabled_) {
            double variability = calculateVariabilityFactor(element.position, time_step);
            element.emission = element.emission * variability;
        }
        
        // Apply cooling losses
        applyCooling(element, time_step);
        
        // Update position (jet expansion)
        Vector4 velocity = calculateJetVelocity(element.position, jet_model_);
        element.position = element.position + velocity * (time_step / SpectralConstants::SPEED_OF_LIGHT);
    }
}

// Private helper functions implementation

double RelativisticJets::lorentzFactor(const Vector4& velocity) const {
    double v_mag = VectorMath::magnitude3D(velocity);
    double beta = v_mag / SpectralConstants::SPEED_OF_LIGHT;
    beta = std::min(beta, 0.9999);  // Avoid division by zero
    
    return 1.0 / std::sqrt(1.0 - beta * beta);
}

Vector4 RelativisticJets::fourVelocity(const Vector4& three_velocity) const {
    double gamma = lorentzFactor(three_velocity);
    return Vector4(gamma, gamma * three_velocity.x / SpectralConstants::SPEED_OF_LIGHT,
                   gamma * three_velocity.y / SpectralConstants::SPEED_OF_LIGHT,
                   gamma * three_velocity.z / SpectralConstants::SPEED_OF_LIGHT);
}

Matrix4 RelativisticJets::lorentzBoostMatrix(const Vector4& velocity) const {
    Matrix4 boost;
    
    double gamma = lorentzFactor(velocity);
    Vector4 beta = velocity * (1.0 / SpectralConstants::SPEED_OF_LIGHT);
    double beta_mag = VectorMath::magnitude3D(beta);
    
    if (beta_mag < 1e-10) return boost;  // Identity for small velocities
    
    Vector4 n = VectorMath::normalize3D(beta);
    
    // Lorentz boost matrix along direction n
    boost.elements[0][0] = gamma;
    boost.elements[0][1] = -gamma * beta.x;
    boost.elements[0][2] = -gamma * beta.y;
    boost.elements[0][3] = -gamma * beta.z;
    
    boost.elements[1][0] = -gamma * beta.x;
    boost.elements[2][0] = -gamma * beta.y;
    boost.elements[3][0] = -gamma * beta.z;
    
    for (int i = 1; i < 4; ++i) {
        for (int j = 1; j < 4; ++j) {
            double delta_ij = (i == j) ? 1.0 : 0.0;
            double n_i = (i == 1) ? n.x : (i == 2) ? n.y : n.z;
            double n_j = (j == 1) ? n.x : (j == 2) ? n.y : n.z;
            
            boost.elements[i][j] = delta_ij + (gamma - 1.0) * n_i * n_j;
        }
    }
    
    return boost;
}

double RelativisticJets::calculateSynchrotronPower(double electron_energy, double magnetic_field) const {
    // Synchrotron power: P = (2/3) * r_e * c * β² * γ² * B²
    double gamma = electron_energy;
    double beta = std::sqrt(1.0 - 1.0/(gamma*gamma));
    
    return (2.0/3.0) * CLASSICAL_ELECTRON_RADIUS * SpectralConstants::SPEED_OF_LIGHT * 
           beta * beta * gamma * gamma * magnetic_field * magnetic_field;
}

double RelativisticJets::synchrotronCoolingTime(double electron_energy, double magnetic_field) const {
    double power = calculateSynchrotronPower(electron_energy, magnetic_field);
    double gamma = electron_energy;
    double rest_energy = ELECTRON_MASS * SpectralConstants::SPEED_OF_LIGHT * SpectralConstants::SPEED_OF_LIGHT;
    
    return gamma * rest_energy / power;
}

void RelativisticJets::applyCooling(JetElement& element, double time_step) const {
    // Simplified cooling - reduce temperature and modify emission
    double cooling_time = synchrotronCoolingTime(jet_model_.min_electron_energy, element.magnetic_field);
    double cooling_factor = std::exp(-time_step / cooling_time);
    
    element.temperature *= cooling_factor;
    element.emission = element.emission * cooling_factor;
}

double RelativisticJets::calculateVariabilityFactor(const Vector4& position, double time) const {
    if (!variability_enabled_) return 1.0;
    
    // Simple turbulent variability model
    double x = VectorMath::magnitude3D(position) / physics_.outerHorizonRadius();
    double phase = time / variability_timescale_ + x * 0.1;
    
    return 1.0 + 0.3 * std::sin(2.0 * Constants::PI * phase) * std::cos(Constants::PI * phase * 0.7);
}

// Namespace implementations
namespace JetPhysics {
    double synchrotronSelfAbsorptionFrequency(double magnetic_field, double density) {
        // Simplified self-absorption frequency
        return 1e9 * std::pow(magnetic_field, 2.5) * std::pow(density, 0.5);
    }
    
    double plasmaCutoffFrequency(double density) {
        // Plasma frequency: ω_p = √(e²n/(ε₀m_e))
        const double epsilon_0 = 8.8541878128e-12;  // F/m
        return std::sqrt(RelativisticJets::ELECTRON_CHARGE * RelativisticJets::ELECTRON_CHARGE * density / 
                        (epsilon_0 * RelativisticJets::ELECTRON_MASS)) / (2.0 * Constants::PI);
    }
    
    double characteristicSynchrotronFrequency(double electron_energy, double magnetic_field) {
        double gamma = electron_energy;
        return (3.0/2.0) * (RelativisticJets::ELECTRON_CHARGE * magnetic_field * gamma * gamma) / 
               (2.0 * Constants::PI * RelativisticJets::ELECTRON_MASS * SpectralConstants::SPEED_OF_LIGHT);
    }
}
