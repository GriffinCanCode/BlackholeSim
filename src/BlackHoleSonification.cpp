#include "BlackHoleSonification.h"
#include <cmath>
#include <algorithm>
#include <random>
#include <sstream>
#include <iomanip>
#include <iostream>

// Standard library for advanced math
using namespace std;

BlackHoleSonification::BlackHoleSonification(const BlackHolePhysics& physics,
                                           const SimpleLightSystem& light_sim,
                                           const RelativisticJets& jets)
    : physics_(physics), light_sim_(light_sim), jets_(jets),
      melodious_calc_(MusicalScale::MAJOR_PENTATONIC),
      volume_scale_(1.0), min_freq_(SonificationConstants::MIN_AUDIBLE_FREQ),
      max_freq_(SonificationConstants::MAX_AUDIBLE_FREQ),
      relativistic_effects_(true), spatial_audio_(true),
      dynamic_range_db_(SonificationConstants::DYNAMIC_RANGE_DB),
      current_scale_(MusicalScale::MAJOR_PENTATONIC),
      mapping_method_(FrequencyMappingMethod::QUANTIZED),
      use_melodious_frequencies_(true),
      high_freq_rolloff_(8000.0),           // Default 8kHz rolloff
      high_freq_reduction_db_(-6.0),       // 6dB reduction
      compression_threshold_(0.7),         // Moderate compression threshold
      compression_ratio_(3.0),             // 3:1 compression ratio
      compression_attack_(0.01),           // 10ms attack
      compression_release_(0.1),           // 100ms release
      enable_eq_(true),                    // Enable EQ by default
      enable_compression_(true),           // Enable compression by default
      enable_deessing_(true),              // Enable de-essing by default
      cache_valid_(false) {
    
    // Initialize melodious calculator with appropriate frequency range
    melodious_calc_.setFrequencyRange(min_freq_, max_freq_);
    
    updatePhysicsCache();
}

// ====== CORE SONIFICATION FUNCTIONS ======

BlackHoleAudioSpectrum BlackHoleSonification::generateAudioSpectrum(const SonificationContext& context) const {
    BlackHoleAudioSpectrum spectrum;
    
    // Generate each physics-based audio component
    spectrum.gravitational_wave = generateGravitationalWaveAudio(context);
    spectrum.electromagnetic = generateElectromagneticAudio(context);
    spectrum.accretion_disk = generateAccretionDiskAudio(context.observer_position);
    spectrum.relativistic_jets = generateJetAudio(context.observer_position, context.observer_position);
    spectrum.hawking_radiation = generateHawkingRadiationAudio();
    
    // Apply relativistic effects if enabled
    if (relativistic_effects_) {
        applyGravitationalRedshiftToSpectrum(spectrum, context);
        applyDopplerShiftToSpectrum(spectrum, context);
        applyTimeDilationToSpectrum(spectrum, context);
    }
    
    // Apply special region effects
    if (context.approaching_horizon) {
        applyEventHorizonEffects(spectrum, context.schwarzschild_ratio);
    }
    
    if (context.inside_ergosphere) {
        applyErgosphereEffects(spectrum, context.observer_position);
    }
    
    // Calculate composite properties
    calculateCompositeSpectrum(spectrum);
    
    // Apply spatial audio effects
    if (spatial_audio_) {
        applySpatialEffects(spectrum, context.observer_position);
    }
    
    // Normalize and compress dynamic range
    normalizeAudioSpectrum(spectrum);
    applyDynamicRangeCompression(spectrum);
    
    // Apply advanced audio processing to reduce harshness
    if (enable_eq_) {
        applyFrequencyEQ(spectrum);
    }
    
    if (enable_compression_) {
        applyMultibandCompression(spectrum);
    }
    
    if (enable_deessing_) {
        applyDeEssing(spectrum);
    }
    
    // Apply high frequency rolloff last to smooth everything
    applyHighFrequencyRolloff(spectrum);
    
    return spectrum;
}

SonificationContext BlackHoleSonification::updateContext(const Vector4& observer_pos,
                                                        const Vector4& observer_vel,
                                                        double time) const {
    SonificationContext context;
    
    // Basic position and time
    context.observer_position = observer_pos;
    context.observer_velocity = observer_vel;
    context.coordinate_time = time;
    
    // Calculate position-dependent parameters
    double r = VectorMath::magnitude3D(observer_pos);
    double rs = cached_schwarzschild_radius_;
    
    context.schwarzschild_ratio = r / rs;
    context.escape_velocity = calculateEscapeVelocity(r);
    context.tidal_acceleration = calculateTidalAcceleration(observer_pos);
    
    // Proper time calculation using metric
    double g_tt = physics_.metricGtt(r);
    context.proper_time = time * sqrt(-g_tt);
    
    // Region detection
    context.approaching_horizon = (r < 10.0 * rs);
    context.inside_ergosphere = physics_.isRotating() && 
                               physics_.isInErgosphere(r, acos(observer_pos.z / r));
    
    // Check if in accretion disk (simplified: near equatorial plane within disk radius)
    double disk_outer = 20.0 * rs;  // Simple fallback: 20 Schwarzschild radii
    double height = abs(observer_pos.z);
    context.in_disk_region = (r > rs && r < disk_outer && height < 0.5 * r);
    
    // Check if in jet region (simplified: near polar axis)
    Vector4 jet_axis = jets_.getJetModel().jet_axis;
    double alignment = abs(VectorMath::dot3D(VectorMath::normalize3D(observer_pos), jet_axis));
    context.in_jet_region = (alignment > 0.8 && r < jets_.getJetModel().jet_height);
    
    return context;
}

std::vector<AudioSample> BlackHoleSonification::generateAudioSamples(
    const BlackHoleAudioSpectrum& spectrum,
    double duration_seconds,
    int sample_rate) const {
    
    std::vector<AudioSample> samples;
    int num_samples = static_cast<int>(duration_seconds * sample_rate);
    samples.reserve(num_samples);
    
    double dt = 1.0 / sample_rate;
    
    // Phase accumulators for each component to avoid discontinuities
    double gw_phase = spectrum.gravitational_wave.phase;
    double em_phase = spectrum.electromagnetic.phase;
    double disk_phase = spectrum.accretion_disk.phase;
    double jet_phase = spectrum.relativistic_jets.phase;
    double hawking_phase = spectrum.hawking_radiation.phase;
    
    // Calculate proper amplitude scaling to prevent clipping
    double total_component_amplitude = spectrum.gravitational_wave.amplitude + 
                                     spectrum.electromagnetic.amplitude +
                                     spectrum.accretion_disk.amplitude +
                                     spectrum.relativistic_jets.amplitude +
                                     spectrum.hawking_radiation.amplitude;
    
    double amplitude_scale = (total_component_amplitude > 1.0) ? (0.8 / total_component_amplitude) : 0.8;
    
    for (int i = 0; i < num_samples; ++i) {
        double t = i * dt;
        AudioSample sample;
        
        // Generate each component with proper phase continuity
        double gw_component = 0.0, em_component = 0.0, disk_component = 0.0, 
               jet_component = 0.0, hawking_component = 0.0;
        
        if (spectrum.gravitational_wave.amplitude > 0.001) {
            gw_component = spectrum.gravitational_wave.amplitude * amplitude_scale *
                          generateSineWave(spectrum.gravitational_wave.frequency, gw_phase, t);
        }
        
        if (spectrum.electromagnetic.amplitude > 0.001) {
            // Use sine wave instead of harsh sawtooth for electromagnetic
            em_component = spectrum.electromagnetic.amplitude * amplitude_scale *
                          generateSineWave(spectrum.electromagnetic.frequency, em_phase, t);
        }
        
        if (spectrum.accretion_disk.amplitude > 0.001) {
            disk_component = spectrum.accretion_disk.amplitude * amplitude_scale *
                           generateSineWave(spectrum.accretion_disk.frequency, disk_phase, t);
        }
        
        if (spectrum.relativistic_jets.amplitude > 0.001) {
            jet_component = spectrum.relativistic_jets.amplitude * amplitude_scale *
                           generateSineWave(spectrum.relativistic_jets.frequency, jet_phase, t);
        }
        
        if (spectrum.hawking_radiation.amplitude > 0.001) {
            // Very quiet thermal noise for Hawking radiation
            hawking_component = spectrum.hawking_radiation.amplitude * amplitude_scale * 0.1 *
                               generateNoiseWave(1.0, t);
        }
        
        // Mix components with proper weighting
        double composite_amplitude = gw_component + em_component + disk_component + 
                                   jet_component + hawking_component;
        
        // Apply volume scaling with soft limiting to prevent harsh clipping
        double final_amplitude = composite_amplitude * volume_scale_;
        final_amplitude = tanh(final_amplitude);  // Soft limiting
        
        sample.amplitude = clamp(final_amplitude, -1.0, 1.0);
        sample.frequency = spectrum.dominant_frequency;
        sample.phase = 2.0 * Constants::PI * spectrum.dominant_frequency * t;
        sample.duration = dt;
        
        samples.push_back(sample);
    }
    
    return samples;
}

// ====== PHYSICS-TO-FREQUENCY CONVERSIONS ======
// All formulas based on established physics

double BlackHoleSonification::orbitalFrequencyToAudio(double orbital_radius) const {
    if (!cache_valid_) updatePhysicsCache();
    
    // Kepler's third law: f_orbital = (1/2π) * √(GM/r³)
    double GM = Constants::G * cached_mass_ * Constants::SOLAR_MASS;
    double orbital_freq = (1.0 / (2.0 * Constants::PI)) * sqrt(GM / pow(orbital_radius, 3));
    
    if (use_melodious_frequencies_) {
        // Use melodious calculator to map physics frequency to musical scale
        double physics_min = 1e-6;  // Hz - very low orbital frequencies
        double physics_max = 1e3;   // Hz - high orbital frequencies
        return melodious_calc_.convertPhysicsToAudio(orbital_freq, physics_min, physics_max, current_scale_);
    } else {
        // Original scaling method
        double scaled_freq = SonificationConstants::ORBITAL_FREQ_SCALE * 
                            log(orbital_freq / (1e-6)) / log(1e6);
        return clampFrequency(scaled_freq);
    }
}

double BlackHoleSonification::electromagneticFrequencyToAudio(double wavelength_meters) const {
    // EM frequency: f = c/λ
    double em_freq = SpectralConstants::SPEED_OF_LIGHT / wavelength_meters;
    
    if (use_melodious_frequencies_) {
        // Use melodious calculator to map EM frequency to musical scale
        double physics_min = 1e6;   // Hz - radio waves
        double physics_max = 1e18;  // Hz - X-rays
        return melodious_calc_.convertPhysicsToAudio(em_freq, physics_min, physics_max, current_scale_);
    } else {
        // Original logarithmic scaling from EM spectrum to audio
        // X-rays (~10^18 Hz) -> high audio, Radio (~10^6 Hz) -> low audio
        double log_em_freq = log10(em_freq);
        double log_audio_min = log10(min_freq_);
        double log_audio_max = log10(max_freq_);
        
        // Map from electromagnetic range (log10(10^6) to log10(10^18)) to audio range
        double scaled_freq = pow(10, log_audio_min + 
                                (log_em_freq - 6.0) / (18.0 - 6.0) * 
                                (log_audio_max - log_audio_min));
        
        return clampFrequency(scaled_freq);
    }
}

double BlackHoleSonification::temperatureToAudioFrequency(double temperature_kelvin) const {
    // Thermal frequency: f = k_B*T/h
    double thermal_freq = SpectralConstants::BOLTZMANN_CONSTANT * temperature_kelvin / 
                         SpectralConstants::PLANCK_CONSTANT;
    
    if (use_melodious_frequencies_) {
        // Use melodious calculator for temperature-based frequencies
        double physics_min = 1e10;  // Hz - very cold temperatures
        double physics_max = 1e14;  // Hz - very hot temperatures
        return melodious_calc_.convertPhysicsToAudio(thermal_freq, physics_min, physics_max, current_scale_);
    } else {
        // Original scaling using temperature-frequency scaling constant
        double scaled_freq = SonificationConstants::TEMP_FREQ_SCALE * 
                            sqrt(thermal_freq / 1e12);  // Geometric scaling
        
        return clampFrequency(scaled_freq);
    }
}

double BlackHoleSonification::synchrotronFrequencyToAudio(double electron_energy, double magnetic_field) const {
    // Synchrotron frequency: f_sync = (eB/2πm_e) * γ²
    constexpr double ELECTRON_CHARGE = 1.602176634e-19;  // C
    constexpr double ELECTRON_MASS = 9.1093837015e-31;   // kg
    
    double cyclotron_freq = (ELECTRON_CHARGE * magnetic_field) / (2.0 * Constants::PI * ELECTRON_MASS);
    double synch_freq = cyclotron_freq * electron_energy * electron_energy;  // γ² factor
    
    if (use_melodious_frequencies_) {
        // Use melodious calculator for synchrotron frequencies
        double physics_min = 1e9;   // Hz - low-energy electrons
        double physics_max = 1e15;  // Hz - high-energy electrons
        return melodious_calc_.convertPhysicsToAudio(synch_freq, physics_min, physics_max, current_scale_);
    } else {
        // Original scaling to audio range
        double scaled_freq = SonificationConstants::SYNCH_FREQ_SCALE * 
                            log(synch_freq / 1e9) / log(1e6);
        
        return clampFrequency(scaled_freq);
    }
}

double BlackHoleSonification::gravitationalWaveFrequency(double orbital_radius) const {
    // Gravitational waves: f_GW = 2 * f_orbital (quadrupole emission)
    double orbital_freq_hz = calculateOrbitalFrequency(orbital_radius);
    double gw_freq = 2.0 * orbital_freq_hz;
    
    if (use_melodious_frequencies_) {
        // Use melodious calculator for gravitational wave frequencies
        double physics_min = 1e-6;  // Hz - very low GW frequencies
        double physics_max = 1e3;   // Hz - high GW frequencies (LIGO range)
        return melodious_calc_.convertPhysicsToAudio(gw_freq, physics_min, physics_max, current_scale_);
    } else {
        // Original scaling to audible range
        double scaled_freq = SonificationConstants::ORBITAL_FREQ_SCALE * 
                            pow(gw_freq / 1e3, 0.3);  // Power-law scaling
        
        return clampFrequency(scaled_freq);
    }
}

// ====== RELATIVISTIC EFFECTS ON SOUND ======

double BlackHoleSonification::applyGravitationalRedshift(double base_frequency,
                                                        double emission_radius,
                                                        double observer_radius) const {
    // Gravitational redshift: f_obs = f_emit * √(-g_tt(r_obs) / -g_tt(r_emit))
    double g_tt_emission = physics_.metricGtt(emission_radius);
    double g_tt_observer = physics_.metricGtt(observer_radius);
    
    if (g_tt_emission >= 0 || g_tt_observer >= 0) {
        return SonificationConstants::REDSHIFT_FREQ_LIMIT;  // Near/inside photon sphere
    }
    
    double redshift_factor = sqrt(-g_tt_observer / -g_tt_emission);
    double shifted_freq = base_frequency * redshift_factor;
    
    return clampFrequency(shifted_freq);
}

double BlackHoleSonification::applyDopplerShift(double base_frequency,
                                              const Vector4& source_velocity,
                                              const Vector4& observer_direction) const {
    // Relativistic Doppler: f_obs = f_emit * γ(1 + β·n̂)
    double gamma = lorentzFactor(source_velocity);
    double beta = VectorMath::magnitude3D(source_velocity) / SpectralConstants::SPEED_OF_LIGHT;
    
    Vector4 norm_vel = VectorMath::normalize3D(source_velocity);
    Vector4 norm_obs = VectorMath::normalize3D(observer_direction);
    double cos_angle = VectorMath::dot3D(norm_vel, norm_obs);
    
    double doppler_factor = gamma * (1.0 + beta * cos_angle);
    double shifted_freq = base_frequency * doppler_factor;
    
    // Limit extreme shifts
    shifted_freq = std::max(shifted_freq, SonificationConstants::REDSHIFT_FREQ_LIMIT);
    shifted_freq = std::min(shifted_freq, base_frequency * SonificationConstants::BLUESHIFT_FREQ_LIMIT);
    
    return clampFrequency(shifted_freq);
}

double BlackHoleSonification::calculateFrameDraggingPan(const Vector4& position) const {
    if (!physics_.isRotating()) return 0.0;
    
    double r = VectorMath::magnitude3D(position);
    Vector4 spherical = physics_.cartesianToSpherical(position);
    double theta = spherical.y;
    
    // Lense-Thirring precession: Ω = 2Ma/(r³ + a²r + 2Ma²)
    double a = physics_.getSpinParameter();
    double M = cached_mass_;
    double rs = cached_schwarzschild_radius_;
    
    double omega_LT = (2.0 * M * a * rs) / (r*r*r + a*a*r + 2.0*M*a*a*rs);
    
    // Convert to stereo pan (-1 to 1)
    double max_omega = 0.1;  // Normalize to reasonable range
    double pan = tanh(omega_LT / max_omega) * sin(theta);
    
    return clamp(pan, -1.0, 1.0);
}

// ====== SPECIALIZED PHYSICS AUDIO ======

AudioSample BlackHoleSonification::generateAccretionDiskAudio(const Vector4& position) const {
    double r = VectorMath::magnitude3D(position);
    AudioSample sample;
    
    // Shakura-Sunyaev temperature profile: T(r) ∝ r^(-3/4)
    double rs = cached_schwarzschild_radius_;
    double r_isco = 3.0 * rs;  // Innermost stable circular orbit
    
    if (r < r_isco) {
        // Inside ISCO - chaotic, high-energy
        sample.frequency = temperatureToAudioFrequency(1e8);  // 100 million K
        sample.amplitude = 0.8;
    } else {
        // Standard disk temperature profile
        double central_temp = 3e7;  // 30 million Kelvin fallback
        double temperature = central_temp * pow(r_isco / r, 0.75);
        sample.frequency = temperatureToAudioFrequency(temperature);
        
        // Orbital velocity for amplitude: v ∝ r^(-1/2)
        double orbital_freq = orbitalFrequencyToAudio(r);
        sample.amplitude = 0.3 * sqrt(r_isco / r);
    }
    
    sample.phase = 0.0;
    sample.duration = 0.1;
    
    return sample;
}

AudioSample BlackHoleSonification::generateJetAudio(const Vector4& position, const Vector4& observer) const {
    AudioSample sample;
    
    const RelativisticJetModel& jet_model = jets_.getJetModel();
    
    if (!jets_.isPointInJet(position, jet_model)) {
        sample.amplitude = 0.0;
        sample.frequency = 440.0;
        return sample;
    }
    
    // Jet bulk motion creates Doppler-shifted synchrotron emission
    Vector4 jet_velocity = jets_.calculateJetVelocity(position, jet_model);
    double jet_speed = VectorMath::magnitude3D(jet_velocity);
    
    // Synchrotron emission from power-law electrons
    double electron_energy = jet_model.min_electron_energy * 100.0;  // Characteristic energy
    double sync_freq = synchrotronFrequencyToAudio(electron_energy, jet_model.magnetic_field);
    
    // Apply relativistic beaming
    Vector4 observer_dir(observer.x - position.x, observer.y - position.y, 
                        observer.z - position.z, observer.t - position.t);
    double beaming_factor = jets_.calculateBeamingFactor(jet_velocity, observer_dir);
    
    sample.frequency = sync_freq * sqrt(beaming_factor);  // Frequency boost
    sample.amplitude = 0.4 * log(1.0 + beaming_factor) / log(10.0);  // Amplitude boost
    sample.phase = 0.0;
    sample.duration = 0.05;  // Short bursts characteristic of jet variability
    
    return sample;
}

AudioSample BlackHoleSonification::generateHawkingRadiationAudio() const {
    AudioSample sample;
    
    // Hawking temperature: T_H = ℏc³/(8πGMk_B)
    double hbar = SpectralConstants::PLANCK_CONSTANT / (2.0 * Constants::PI);
    double c3 = pow(SpectralConstants::SPEED_OF_LIGHT, 3);
    double GM = Constants::G * cached_mass_ * Constants::SOLAR_MASS;
    double kB = SpectralConstants::BOLTZMANN_CONSTANT;
    
    double T_hawking = (hbar * c3) / (8.0 * Constants::PI * GM * kB);
    
    // For stellar mass black holes, T_H ~ 10^(-7) K - incredibly cold!
    // Scale dramatically for audible representation
    double scaled_temp = T_hawking * 1e15;  // Massive scaling factor
    
    sample.frequency = temperatureToAudioFrequency(scaled_temp);
    sample.amplitude = 0.01;  // Very quiet, like the real thing
    sample.phase = 0.0;
    sample.duration = 1.0;    // Long duration like thermal noise
    
    return sample;
}

AudioSample BlackHoleSonification::generateTidalForceAudio(const Vector4& position) const {
    AudioSample sample;
    
    double tidal_accel = calculateTidalAcceleration(position);
    
    // Convert tidal acceleration to rhythmic frequency
    // Higher tidal forces = faster rhythm
    sample.frequency = 50.0 + 200.0 * tanh(tidal_accel / 1e10);  // Scale to 50-250 Hz range
    sample.amplitude = 0.2 * sqrt(tidal_accel / 1e8);
    sample.phase = 0.0;
    sample.duration = 0.2;
    
    return sample;
}

// ====== DYNAMIC AUDIO EFFECTS ======

double BlackHoleSonification::calculateTimeDilationFactor(const Vector4& position) const {
    double r = VectorMath::magnitude3D(position);
    double g_tt = physics_.metricGtt(r);
    
    // Time dilation factor: dt_proper = dt_coordinate * √(-g_tt)
    return sqrt(-g_tt);
}

void BlackHoleSonification::applyEventHorizonEffects(BlackHoleAudioSpectrum& spectrum,
                                                    double schwarzschild_ratio) const {
    if (schwarzschild_ratio > 10.0) return;  // Far from horizon
    
    // Dramatic redshift as r → r_s
    double redshift_factor = sqrt(schwarzschild_ratio - 1.0);
    redshift_factor = std::max(redshift_factor, 0.01);  // Prevent zero
    
    // Apply to all frequencies
    spectrum.gravitational_wave.frequency *= redshift_factor;
    spectrum.electromagnetic.frequency *= redshift_factor;
    spectrum.accretion_disk.frequency *= redshift_factor;
    spectrum.relativistic_jets.frequency *= redshift_factor;
    
    // Increase amplitude of low-frequency components (redshift effect)
    spectrum.gravitational_wave.amplitude *= (1.5 / redshift_factor);
    
    // Add dramatic time dilation effects
    if (schwarzschild_ratio < 2.0) {
        // Near horizon - extreme time dilation
        double time_stretch = 10.0 / schwarzschild_ratio;
        spectrum.total_amplitude *= 0.1;  // Quieter due to time dilation
        
        // All frequencies shift way down
        spectrum.dominant_frequency /= time_stretch;
    }
}

void BlackHoleSonification::applyErgosphereEffects(BlackHoleAudioSpectrum& spectrum,
                                                  const Vector4& position) const {
    // Frame-dragging creates swirling stereo effects
    spectrum.stereo_pan = calculateFrameDraggingPan(position);
    
    // Add frequency modulation from frame-dragging
    double frame_drag_freq = 0.5;  // Hz, slow swirling
    double mod_depth = 0.1;
    
    spectrum.gravitational_wave.frequency *= (1.0 + mod_depth * sin(2.0 * Constants::PI * frame_drag_freq * position.t));
    spectrum.electromagnetic.frequency *= (1.0 + mod_depth * sin(2.0 * Constants::PI * frame_drag_freq * position.t + Constants::PI/3));
}

// ====== HELPER FUNCTIONS ======

double BlackHoleSonification::calculateOrbitalFrequency(double radius) const {
    // Kepler's third law for circular orbits
    double GM = Constants::G * cached_mass_ * Constants::SOLAR_MASS;
    return (1.0 / (2.0 * Constants::PI)) * sqrt(GM / pow(radius, 3));
}

double BlackHoleSonification::calculateEscapeVelocity(double radius) const {
    double GM = Constants::G * cached_mass_ * Constants::SOLAR_MASS;
    return sqrt(2.0 * GM / radius);
}

double BlackHoleSonification::calculateTidalAcceleration(const Vector4& position) const {
    double r = VectorMath::magnitude3D(position);
    double GM = Constants::G * cached_mass_ * Constants::SOLAR_MASS;
    
    // Tidal acceleration: a_tidal = 2GM/r³ * Δr (assume Δr ~ 1m for human scale)
    return (2.0 * GM / pow(r, 3)) * 1.0;  // 1 meter separation
}

double BlackHoleSonification::lorentzFactor(const Vector4& velocity) const {
    double v_mag = VectorMath::magnitude3D(velocity);
    double beta = v_mag / SpectralConstants::SPEED_OF_LIGHT;
    return 1.0 / sqrt(1.0 - beta * beta);
}

// ====== WAVEFORM GENERATORS ======

double BlackHoleSonification::generateSineWave(double frequency, double phase, double time) const {
    return sin(2.0 * Constants::PI * frequency * time + phase);
}

double BlackHoleSonification::generateSquareWave(double frequency, double phase, double time) const {
    double sine_val = sin(2.0 * Constants::PI * frequency * time + phase);
    return (sine_val >= 0.0) ? 1.0 : -1.0;
}

double BlackHoleSonification::generateSawtoothWave(double frequency, double phase, double time) const {
    double period = 1.0 / frequency;
    double t_mod = fmod(time + phase / (2.0 * Constants::PI * frequency), period);
    return 2.0 * (t_mod / period) - 1.0;
}

double BlackHoleSonification::generateNoiseWave(double amplitude, double time) const {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(-1.0, 1.0);
    return amplitude * dis(gen);
}

double BlackHoleSonification::accretionDiskWaveform(double time, double radius) const {
    // Combine orbital motion with thermal fluctuations - smoother version
    double orbital_freq = orbitalFrequencyToAudio(radius);
    double orbital_phase = 2.0 * Constants::PI * orbital_freq * time;
    
    // Add harmonics for richer sound but keep it musical
    double fundamental = 0.6 * sin(orbital_phase);
    double second_harmonic = 0.3 * sin(2.0 * orbital_phase);
    double third_harmonic = 0.15 * sin(3.0 * orbital_phase);
    
    // Very quiet thermal fluctuation (filtered noise)
    static double prev_noise = 0.0;
    double new_noise = generateNoiseWave(1.0, time);
    double filtered_noise = 0.9 * prev_noise + 0.1 * new_noise;  // Simple low-pass filter
    prev_noise = filtered_noise;
    
    return fundamental + second_harmonic + third_harmonic + 0.05 * filtered_noise;
}

double BlackHoleSonification::jetEmissionWaveform(double time, const Vector4& position) const {
    // Smoother jet emission with less harsh bursts
    double burst_freq = 3.0;  // Slower 3 Hz burst rate for less harsh sound
    double burst_phase = 2.0 * Constants::PI * burst_freq * time;
    
    // Smoother envelope using cosine-based function instead of sharp exponential
    double t_in_cycle = fmod(time, 1.0/burst_freq);
    double cycle_progress = t_in_cycle * burst_freq;
    double envelope = (cycle_progress < 0.5) ? 
                      0.5 * (1.0 + cos(2.0 * Constants::PI * cycle_progress)) :
                      0.1 * exp(-10.0 * (cycle_progress - 0.5));  // Gentle decay
    
    // Use sine wave with subtle modulation instead of harsh noise
    double carrier = sin(burst_phase);
    double modulation = 1.0 + 0.1 * sin(13.7 * burst_phase);  // Subtle frequency modulation
    
    return envelope * carrier * modulation;
}

double BlackHoleSonification::gravitationalWaveform(double time, double orbital_freq) const {
    // Gravitational wave chirp - more realistic and musical
    double chirp_rate = 0.05;  // Even slower, less aggressive chirp
    double instantaneous_freq = orbital_freq * (1.0 + chirp_rate * time);
    
    // Two-polarization gravitational wave (+ and × modes)
    double plus_mode = sin(2.0 * Constants::PI * instantaneous_freq * time);
    double cross_mode = sin(2.0 * Constants::PI * instantaneous_freq * time + Constants::PI/2);
    
    // Blend the polarizations for a richer sound
    double blend_factor = 0.5 + 0.3 * sin(0.1 * time);  // Slowly varying blend
    return 0.4 * (blend_factor * plus_mode + (1.0 - blend_factor) * cross_mode);
}

double BlackHoleSonification::hawkingRadiationWaveform(double time) const {
    // Thermal noise characteristic of black body radiation - heavily filtered
    static double prev_noise1 = 0.0, prev_noise2 = 0.0;
    double new_noise = generateNoiseWave(1.0, time);
    
    // Two-stage low-pass filter to make thermal noise very smooth
    double filtered1 = 0.95 * prev_noise1 + 0.05 * new_noise;
    double filtered2 = 0.95 * prev_noise2 + 0.05 * filtered1;
    
    prev_noise1 = filtered1;
    prev_noise2 = filtered2;
    
    return filtered2 * 0.05;  // Even quieter than before
}

// ====== UTILITY FUNCTIONS ======

double BlackHoleSonification::clampFrequency(double frequency) const {
    // Ensure frequency is within audible range and below Nyquist frequency for 44.1kHz
    const double NYQUIST_FREQ = 22050.0;  // Half of 44.1kHz sample rate
    const double SAFE_MAX_FREQ = std::min(max_freq_, NYQUIST_FREQ * 0.8);  // 80% of Nyquist for safety
    
    // Also ensure minimum frequency isn't too low (below human hearing)
    const double SAFE_MIN_FREQ = std::max(min_freq_, 40.0);  // Minimum 40Hz
    
    return std::max(SAFE_MIN_FREQ, std::min(SAFE_MAX_FREQ, frequency));
}

double BlackHoleSonification::logarithmicScale(double input_freq, double min_in, double max_in,
                                             double min_out, double max_out) const {
    double log_in = log(input_freq / min_in) / log(max_in / min_in);
    return min_out * pow(max_out / min_out, log_in);
}

void BlackHoleSonification::updatePhysicsCache() const {
    cached_schwarzschild_radius_ = physics_.schwarzschildRadius();
    cached_mass_ = physics_.getMass();
    cached_spin_parameter_ = physics_.getSpinParameter();
    cache_valid_ = true;
}

void BlackHoleSonification::normalizeAudioSpectrum(BlackHoleAudioSpectrum& spectrum) const {
    // Find peak amplitude
    double max_amp = std::max({
        spectrum.gravitational_wave.amplitude,
        spectrum.electromagnetic.amplitude,
        spectrum.accretion_disk.amplitude,
        spectrum.relativistic_jets.amplitude,
        spectrum.hawking_radiation.amplitude
    });
    
    if (max_amp > 1.0) {
        double scale = 1.0 / max_amp;
        spectrum.gravitational_wave.amplitude *= scale;
        spectrum.electromagnetic.amplitude *= scale;
        spectrum.accretion_disk.amplitude *= scale;
        spectrum.relativistic_jets.amplitude *= scale;
        spectrum.hawking_radiation.amplitude *= scale;
    }
    
    // Calculate total amplitude
    spectrum.total_amplitude = spectrum.gravitational_wave.amplitude +
                              spectrum.electromagnetic.amplitude +
                              spectrum.accretion_disk.amplitude +
                              spectrum.relativistic_jets.amplitude +
                              spectrum.hawking_radiation.amplitude;
    
    // Find dominant frequency (weighted by amplitude)
    double freq_sum = 0.0, weight_sum = 0.0;
    
    freq_sum += spectrum.gravitational_wave.frequency * spectrum.gravitational_wave.amplitude;
    weight_sum += spectrum.gravitational_wave.amplitude;
    
    freq_sum += spectrum.electromagnetic.frequency * spectrum.electromagnetic.amplitude;
    weight_sum += spectrum.electromagnetic.amplitude;
    
    freq_sum += spectrum.accretion_disk.frequency * spectrum.accretion_disk.amplitude;
    weight_sum += spectrum.accretion_disk.amplitude;
    
    if (weight_sum > 0.0) {
        spectrum.dominant_frequency = freq_sum / weight_sum;
    } else {
        spectrum.dominant_frequency = SonificationConstants::CONCERT_A;
    }
}

void BlackHoleSonification::calculateCompositeSpectrum(BlackHoleAudioSpectrum& spectrum) const {
    // This is called by generateAudioSpectrum, implement the missing functions it needs
    // For now, implement basic versions
}

void BlackHoleSonification::applyGravitationalRedshiftToSpectrum(BlackHoleAudioSpectrum& spectrum, 
                                                               const SonificationContext& context) const {
    double r_obs = VectorMath::magnitude3D(context.observer_position);
    double r_emit = 10.0 * cached_schwarzschild_radius_;  // Assume emission from 10 r_s
    
    double redshift_factor = sqrt(-physics_.metricGtt(r_obs) / -physics_.metricGtt(r_emit));
    
    spectrum.gravitational_wave.frequency *= redshift_factor;
    spectrum.electromagnetic.frequency *= redshift_factor;
    spectrum.accretion_disk.frequency *= redshift_factor;
}

void BlackHoleSonification::applyDopplerShiftToSpectrum(BlackHoleAudioSpectrum& spectrum,
                                                       const SonificationContext& context) const {
    // Apply Doppler shift based on observer velocity
    double gamma = lorentzFactor(context.observer_velocity);
    double beta = VectorMath::magnitude3D(context.observer_velocity) / SpectralConstants::SPEED_OF_LIGHT;
    
    double doppler_factor = gamma * (1.0 + beta * 0.5);  // Simplified radial component
    
    spectrum.electromagnetic.frequency *= doppler_factor;
    spectrum.relativistic_jets.frequency *= doppler_factor;
}

void BlackHoleSonification::applyTimeDilationToSpectrum(BlackHoleAudioSpectrum& spectrum,
                                                       const SonificationContext& context) const {
    double time_dilation = calculateTimeDilationFactor(context.observer_position);
    
    // Time dilation affects all temporal aspects
    spectrum.gravitational_wave.frequency *= time_dilation;
    spectrum.electromagnetic.frequency *= time_dilation;
}

AudioSample BlackHoleSonification::generateGravitationalWaveAudio(const SonificationContext& context) const {
    AudioSample sample;
    double r = VectorMath::magnitude3D(context.observer_position);
    
    sample.frequency = gravitationalWaveFrequency(r);
    sample.amplitude = 0.3 * exp(-r / (10.0 * cached_schwarzschild_radius_));  // Decay with distance
    sample.phase = 0.0;
    sample.duration = 0.1;
    
    return sample;
}

AudioSample BlackHoleSonification::generateElectromagneticAudio(const SonificationContext& context) const {
    AudioSample sample;
    
    // Use accretion disk temperature to generate EM audio
    double r = VectorMath::magnitude3D(context.observer_position);
    double central_temp = 3e7;  // 30 million Kelvin fallback
    double temp = central_temp * pow(3.0 / (r / cached_schwarzschild_radius_), 0.75);
    
    sample.frequency = temperatureToAudioFrequency(temp);
    sample.amplitude = 0.4;
    sample.phase = 0.0;
    sample.duration = 0.1;
    
    return sample;
}

void BlackHoleSonification::applyDynamicRangeCompression(BlackHoleAudioSpectrum& spectrum) const {
    // Apply logarithmic compression for realistic dynamic range
    double compress_ratio = 0.6;
    
    spectrum.gravitational_wave.amplitude = pow(spectrum.gravitational_wave.amplitude, compress_ratio);
    spectrum.electromagnetic.amplitude = pow(spectrum.electromagnetic.amplitude, compress_ratio);
    spectrum.accretion_disk.amplitude = pow(spectrum.accretion_disk.amplitude, compress_ratio);
    spectrum.relativistic_jets.amplitude = pow(spectrum.relativistic_jets.amplitude, compress_ratio);
    spectrum.hawking_radiation.amplitude = pow(spectrum.hawking_radiation.amplitude, compress_ratio);
}

void BlackHoleSonification::applySpatialEffects(BlackHoleAudioSpectrum& spectrum, const Vector4& position) const {
    // Frame-dragging stereo effect
    spectrum.stereo_pan = calculateFrameDraggingPan(position);
}

std::string BlackHoleSonification::getPhysicsBreakdown(const SonificationContext& context) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);
    
    oss << "=== BLACKHOLE SONIFICATION PHYSICS ===" << std::endl;
    oss << "Position: r/r_s = " << context.schwarzschild_ratio << std::endl;
    oss << "Escape velocity: " << context.escape_velocity / 1000.0 << " km/s" << std::endl;
    oss << "Tidal acceleration: " << context.tidal_acceleration << " m/s²" << std::endl;
    oss << "Time dilation factor: " << calculateTimeDilationFactor(context.observer_position) << std::endl;
    
    if (context.inside_ergosphere) oss << "INSIDE ERGOSPHERE - Frame dragging active" << std::endl;
    if (context.approaching_horizon) oss << "APPROACHING HORIZON - Extreme redshift" << std::endl;
    
    return oss.str();
}

// ====== UTILITY NAMESPACE IMPLEMENTATIONS ======

namespace SonificationUtils {
    std::string frequencyToMusicalNote(double frequency) {
        const std::vector<std::string> notes = {
            "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"
        };
        
        double A4 = SonificationConstants::CONCERT_A;
        double semitones = 12.0 * log2(frequency / A4);
        int note_index = (int(round(semitones)) + 9) % 12;  // A is index 9
        int octave = 4 + (int(round(semitones)) + 9) / 12;
        
        return notes[note_index] + std::to_string(octave);
    }
    
    double calculateBeatFrequency(double freq1, double freq2) {
        return abs(freq1 - freq2);
    }
}

// ====== AUDIO BUFFER CLASS ======

BlackHoleAudioBuffer::BlackHoleAudioBuffer(int buffer_size, int sample_rate)
    : buffer_size_(buffer_size), sample_rate_(sample_rate),
      read_position_(0), write_position_(0), current_time_(0.0), buffer_ready_(false), is_playing_(false)
#ifdef __APPLE__
    , audio_initialized_(false)
#endif
{
    left_buffer_.resize(buffer_size);
    right_buffer_.resize(buffer_size);
    
#ifdef __APPLE__
    if (initializeCoreAudio()) {
        startPlayback();
    }
#endif
}

void BlackHoleAudioBuffer::fillBuffer(const BlackHoleSonification& sonifier,
                                     const SonificationContext& context) {
    BlackHoleAudioSpectrum spectrum = sonifier.generateAudioSpectrum(context);
    
    double dt = 1.0 / sample_rate_;
    double buffer_duration = buffer_size_ * dt;
    
    // Generate proper audio samples using the sonifier's full system
    std::vector<AudioSample> audio_samples = sonifier.generateAudioSamples(spectrum, buffer_duration, sample_rate_);
    
    for (int i = 0; i < buffer_size_ && i < audio_samples.size(); ++i) {
        double mono_sample = audio_samples[i].amplitude;
        
        // Apply stereo panning with smooth transitions
        double pan = spectrum.stereo_pan;
        double left_gain = sqrt((1.0 - pan) / 2.0);
        double right_gain = sqrt((1.0 + pan) / 2.0);
        
        left_buffer_[i] = mono_sample * left_gain;
        right_buffer_[i] = mono_sample * right_gain;
    }
    
    // Fill any remaining buffer space with gentle fade to silence (avoid clicks)
    for (int i = audio_samples.size(); i < buffer_size_; ++i) {
        double fade = exp(-5.0 * (i - audio_samples.size()) / buffer_size_);
        left_buffer_[i] = left_buffer_[std::max(0, i-1)] * fade * 0.9;
        right_buffer_[i] = right_buffer_[std::max(0, i-1)] * fade * 0.9;
    }
    
    buffer_ready_ = true;
    current_time_ += buffer_size_ * dt;
}

void BlackHoleAudioBuffer::getNextSamples(std::vector<double>& left_channel,
                                         std::vector<double>& right_channel,
                                         int num_samples) {
    left_channel.clear();
    right_channel.clear();
    
    for (int i = 0; i < num_samples && read_position_ < buffer_size_; ++i) {
        left_channel.push_back(left_buffer_[read_position_]);
        right_channel.push_back(right_buffer_[read_position_]);
        read_position_++;
    }
    
    if (read_position_ >= buffer_size_) {
        read_position_ = 0;
        buffer_ready_ = false;  // Need to refill
    }
}

BlackHoleAudioBuffer::~BlackHoleAudioBuffer() {
#ifdef __APPLE__
    stopPlayback();
    cleanupCoreAudio();
#endif
}

void BlackHoleAudioBuffer::startPlayback() {
#ifdef __APPLE__
    if (audio_initialized_ && !is_playing_) {
        OSStatus status = AudioUnitInitialize(audio_unit_);
        if (status == noErr) {
            status = AudioOutputUnitStart(audio_unit_);
            if (status == noErr) {
                is_playing_ = true;
                std::cout << "🎵 Audio playback started!" << std::endl;
            } else {
                std::cerr << "🎵 Failed to start audio playback: " << status << std::endl;
            }
        } else {
            std::cerr << "🎵 Failed to initialize audio unit: " << status << std::endl;
        }
    }
#endif
}

void BlackHoleAudioBuffer::stopPlayback() {
#ifdef __APPLE__
    if (audio_initialized_ && is_playing_) {
        AudioOutputUnitStop(audio_unit_);
        AudioUnitUninitialize(audio_unit_);
        is_playing_ = false;
        std::cout << "🎵 Audio playback stopped!" << std::endl;
    }
#endif
}

#ifdef __APPLE__
bool BlackHoleAudioBuffer::initializeCoreAudio() {
    OSStatus status;
    
    // Create audio component description
    AudioComponentDescription desc;
    desc.componentType = kAudioUnitType_Output;
    desc.componentSubType = kAudioUnitSubType_DefaultOutput;
    desc.componentManufacturer = kAudioUnitManufacturer_Apple;
    desc.componentFlags = 0;
    desc.componentFlagsMask = 0;
    
    // Find the default audio component
    AudioComponent component = AudioComponentFindNext(NULL, &desc);
    if (!component) {
        std::cerr << "🎵 Could not find default audio component!" << std::endl;
        return false;
    }
    
    // Create audio unit instance
    status = AudioComponentInstanceNew(component, &audio_unit_);
    if (status != noErr) {
        std::cerr << "🎵 Could not create audio unit instance: " << status << std::endl;
        return false;
    }
    
    // Set up the audio format (stereo 44.1kHz)
    AudioStreamBasicDescription audioFormat;
    audioFormat.mSampleRate = sample_rate_;
    audioFormat.mFormatID = kAudioFormatLinearPCM;
    audioFormat.mFormatFlags = kAudioFormatFlagIsFloat | kAudioFormatFlagIsPacked | kAudioFormatFlagIsNonInterleaved;
    audioFormat.mBytesPerPacket = sizeof(Float32);
    audioFormat.mFramesPerPacket = 1;
    audioFormat.mBytesPerFrame = sizeof(Float32);
    audioFormat.mChannelsPerFrame = 2; // Stereo
    audioFormat.mBitsPerChannel = 32;
    audioFormat.mReserved = 0;
    
    // Set the audio format on the audio unit
    status = AudioUnitSetProperty(audio_unit_,
                                 kAudioUnitProperty_StreamFormat,
                                 kAudioUnitScope_Input,
                                 0,
                                 &audioFormat,
                                 sizeof(audioFormat));
    if (status != noErr) {
        std::cerr << "🎵 Could not set audio format: " << status << std::endl;
        AudioComponentInstanceDispose(audio_unit_);
        return false;
    }
    
    // Set up the render callback
    AURenderCallbackStruct callback;
    callback.inputProc = audioCallback;
    callback.inputProcRefCon = this;
    
    status = AudioUnitSetProperty(audio_unit_,
                                 kAudioUnitProperty_SetRenderCallback,
                                 kAudioUnitScope_Input,
                                 0,
                                 &callback,
                                 sizeof(callback));
    if (status != noErr) {
        std::cerr << "🎵 Could not set render callback: " << status << std::endl;
        AudioComponentInstanceDispose(audio_unit_);
        return false;
    }
    
    audio_initialized_ = true;
    return true;
}

void BlackHoleAudioBuffer::cleanupCoreAudio() {
    if (audio_initialized_) {
        AudioComponentInstanceDispose(audio_unit_);
        audio_initialized_ = false;
    }
}

OSStatus BlackHoleAudioBuffer::audioCallback(void* inRefCon,
                                            AudioUnitRenderActionFlags* ioActionFlags,
                                            const AudioTimeStamp* inTimeStamp,
                                            UInt32 inBusNumber,
                                            UInt32 inNumberFrames,
                                            AudioBufferList* ioData) {
    BlackHoleAudioBuffer* buffer = static_cast<BlackHoleAudioBuffer*>(inRefCon);
    return buffer->renderCallback(ioActionFlags, inTimeStamp, inBusNumber, inNumberFrames, ioData);
}

OSStatus BlackHoleAudioBuffer::renderCallback(AudioUnitRenderActionFlags* ioActionFlags,
                                             const AudioTimeStamp* inTimeStamp,
                                             UInt32 inBusNumber,
                                             UInt32 inNumberFrames,
                                             AudioBufferList* ioData) {
    // Get audio data pointers
    Float32* leftChannel = (Float32*)ioData->mBuffers[0].mData;
    Float32* rightChannel = (Float32*)ioData->mBuffers[1].mData;
    
    // Fill the audio buffers
    for (UInt32 frame = 0; frame < inNumberFrames; ++frame) {
        if (read_position_ < buffer_size_ && buffer_ready_) {
            // Use the generated audio samples
            leftChannel[frame] = static_cast<Float32>(left_buffer_[read_position_]);
            rightChannel[frame] = static_cast<Float32>(right_buffer_[read_position_]);
            read_position_++;
            
            if (read_position_ >= buffer_size_) {
                read_position_ = 0;
                buffer_ready_ = false;  // Need to refill
            }
        } else {
            // Silence if no buffer ready
            leftChannel[frame] = 0.0f;
            rightChannel[frame] = 0.0f;
        }
    }
    
    return noErr;
}
#endif

// ====== ADVANCED AUDIO PROCESSING IMPLEMENTATIONS ======

void BlackHoleSonification::applyFrequencyEQ(BlackHoleAudioSpectrum& spectrum) const {
    // Apply EQ to reduce harsh high frequencies and enhance pleasant mid frequencies
    
    // Reduce harsh frequencies above 4kHz
    if (spectrum.gravitational_wave.frequency > 4000.0) {
        spectrum.gravitational_wave.amplitude *= 0.7;
    }
    if (spectrum.electromagnetic.frequency > 4000.0) {
        spectrum.electromagnetic.amplitude *= 0.7;
    }
    if (spectrum.accretion_disk.frequency > 4000.0) {
        spectrum.accretion_disk.amplitude *= 0.7;
    }
    if (spectrum.relativistic_jets.frequency > 4000.0) {
        spectrum.relativistic_jets.amplitude *= 0.7;
    }
    
    // Enhance pleasant mid-range frequencies (200Hz - 2kHz)
    if (spectrum.gravitational_wave.frequency >= 200.0 && spectrum.gravitational_wave.frequency <= 2000.0) {
        spectrum.gravitational_wave.amplitude *= 1.2;
    }
    if (spectrum.electromagnetic.frequency >= 200.0 && spectrum.electromagnetic.frequency <= 2000.0) {
        spectrum.electromagnetic.amplitude *= 1.2;
    }
    if (spectrum.accretion_disk.frequency >= 200.0 && spectrum.accretion_disk.frequency <= 2000.0) {
        spectrum.accretion_disk.amplitude *= 1.2;
    }
    
    // Gentle low-end boost for warmth (80Hz - 200Hz)
    if (spectrum.gravitational_wave.frequency >= 80.0 && spectrum.gravitational_wave.frequency <= 200.0) {
        spectrum.gravitational_wave.amplitude *= 1.1;
    }
}

void BlackHoleSonification::applyMultibandCompression(BlackHoleAudioSpectrum& spectrum) const {
    // Apply different compression to different frequency bands
    
    // Low band compression (20Hz - 250Hz) - gentle compression for warmth
    double low_threshold = 0.8;
    double low_ratio = 2.0;
    
    if (spectrum.gravitational_wave.frequency <= 250.0) {
        if (spectrum.gravitational_wave.amplitude > low_threshold) {
            double excess = spectrum.gravitational_wave.amplitude - low_threshold;
            spectrum.gravitational_wave.amplitude = low_threshold + excess / low_ratio;
        }
    }
    
    // Mid band compression (250Hz - 4kHz) - moderate compression for consistency
    double mid_threshold = 0.7;
    double mid_ratio = 3.0;
    
    auto applyMidCompression = [&](AudioSample& sample) {
        if (sample.frequency >= 250.0 && sample.frequency <= 4000.0) {
            if (sample.amplitude > mid_threshold) {
                double excess = sample.amplitude - mid_threshold;
                sample.amplitude = mid_threshold + excess / mid_ratio;
            }
        }
    };
    
    applyMidCompression(spectrum.electromagnetic);
    applyMidCompression(spectrum.accretion_disk);
    applyMidCompression(spectrum.relativistic_jets);
    
    // High band compression (4kHz+) - aggressive compression to reduce harshness
    double high_threshold = 0.5;
    double high_ratio = 4.0;
    
    auto applyHighCompression = [&](AudioSample& sample) {
        if (sample.frequency > 4000.0) {
            if (sample.amplitude > high_threshold) {
                double excess = sample.amplitude - high_threshold;
                sample.amplitude = high_threshold + excess / high_ratio;
            }
        }
    };
    
    applyHighCompression(spectrum.gravitational_wave);
    applyHighCompression(spectrum.electromagnetic);
    applyHighCompression(spectrum.accretion_disk);
    applyHighCompression(spectrum.relativistic_jets);
}

void BlackHoleSonification::applyHighFrequencyRolloff(BlackHoleAudioSpectrum& spectrum) const {
    // Apply smooth high frequency rolloff to reduce harshness
    
    auto applyRolloff = [&](AudioSample& sample) {
        if (sample.frequency > high_freq_rolloff_) {
            // Calculate rolloff amount based on how far above the rolloff frequency
            double freq_ratio = sample.frequency / high_freq_rolloff_;
            double rolloff_factor = 1.0 / (1.0 + (freq_ratio - 1.0) * 2.0); // Gentle rolloff curve
            
            // Apply additional dB reduction
            double db_reduction = std::min(0.0, high_freq_reduction_db_);
            double linear_reduction = pow(10.0, db_reduction / 20.0);
            
            sample.amplitude *= rolloff_factor * linear_reduction;
        }
    };
    
    applyRolloff(spectrum.gravitational_wave);
    applyRolloff(spectrum.electromagnetic);
    applyRolloff(spectrum.accretion_disk);
    applyRolloff(spectrum.relativistic_jets);
    applyRolloff(spectrum.hawking_radiation);
}

void BlackHoleSonification::applyDeEssing(BlackHoleAudioSpectrum& spectrum) const {
    // Reduce harsh sibilant frequencies (typically 5kHz - 8kHz)
    const double sibilant_freq_low = 5000.0;
    const double sibilant_freq_high = 8000.0;
    const double deessing_reduction = 0.6; // Reduce by 40%
    
    auto applyDeEss = [&](AudioSample& sample) {
        if (sample.frequency >= sibilant_freq_low && sample.frequency <= sibilant_freq_high) {
            // Apply more aggressive reduction to harsh sibilant range
            sample.amplitude *= deessing_reduction;
        }
    };
    
    applyDeEss(spectrum.electromagnetic);
    applyDeEss(spectrum.relativistic_jets);
    applyDeEss(spectrum.hawking_radiation);
}

double BlackHoleSonification::lowPassFilter(double input, double cutoff_freq, double sample_rate) const {
    // Simple first-order low-pass filter
    static double prev_output = 0.0;
    
    double rc = 1.0 / (cutoff_freq * 2.0 * Constants::PI);
    double dt = 1.0 / sample_rate;
    double alpha = dt / (rc + dt);
    
    double output = prev_output + alpha * (input - prev_output);
    prev_output = output;
    
    return output;
}

double BlackHoleSonification::highShelfFilter(double input, double frequency, double gain_db) const {
    // High shelf filter for frequency-dependent gain adjustment
    double gain_linear = pow(10.0, gain_db / 20.0);
    
    // Simple implementation - in practice this would use biquad filters
    return input * gain_linear;
}

double BlackHoleSonification::compressor(double input, double threshold, double ratio, 
                                       double attack, double release) const {
    // Simple dynamic range compressor
    static double envelope = 0.0;
    
    double input_level = abs(input);
    
    // Update envelope
    if (input_level > envelope) {
        envelope += (input_level - envelope) * attack;
    } else {
        envelope += (input_level - envelope) * release;
    }
    
    // Apply compression if above threshold
    if (envelope > threshold) {
        double excess = envelope - threshold;
        double compressed_excess = excess / ratio;
        double gain_reduction = (threshold + compressed_excess) / envelope;
        return input * gain_reduction;
    }
    
    return input;
}
