#pragma once
#include "Physics.h"
#include "KerrPhysics.h"
#include "LightSystem.h"
#include "RelativisticJets.h"
#include "MelodiousCalculator.h"
#include <vector>
#include <string>
#ifdef __APPLE__
#include <AudioUnit/AudioUnit.h>
#include <CoreAudio/CoreAudio.h>
#endif

// Forward declarations
class BlackHolePhysics;
class SimpleLightSystem;
class RelativisticJets;

/**
 * SCIENTIFIC SONIFICATION CONSTANTS
 * Based on real astronomical phenomena and established physics
 */
namespace SonificationConstants {
    // Human hearing range (Hz)
    constexpr double MIN_AUDIBLE_FREQ = 20.0;
    constexpr double MAX_AUDIBLE_FREQ = 20000.0;
    
    // Musical constants for harmonic scaling
    constexpr double CONCERT_A = 440.0;           // A4 (Hz)
    constexpr double TWELFTH_ROOT_OF_2 = 1.059463094359; // Semitone ratio
    
    // Physical scaling factors - derived from real blackhole physics
    constexpr double ORBITAL_FREQ_SCALE = 1e12;   // Scale orbital frequencies to audio range
    constexpr double EM_FREQ_SCALE = 1e12;        // Scale electromagnetic frequencies 
    constexpr double TEMP_FREQ_SCALE = 2.08e10;   // k_B*T/h scaling (Boltzmann/Planck)
    constexpr double SYNCH_FREQ_SCALE = 1e9;      // Synchrotron frequency scaling
    
    // Dynamic range constants for realistic audio
    constexpr double DYNAMIC_RANGE_DB = 60.0;     // 60dB dynamic range
    constexpr double REFERENCE_AMPLITUDE = 0.1;    // Reference sound amplitude
    
    // Doppler shift parameters from real physics
    constexpr double REDSHIFT_FREQ_LIMIT = 0.01;  // Minimum frequency after redshift
    constexpr double BLUESHIFT_FREQ_LIMIT = 100.0; // Maximum frequency boost
}

/**
 * AUDIO SAMPLE REPRESENTATION
 * Real-time audio generation with physics-based parameters
 */
struct AudioSample {
    double amplitude;       // Sound amplitude (0.0 to 1.0)
    double frequency;       // Instantaneous frequency (Hz)
    double phase;          // Phase offset (radians)
    double duration;       // Sample duration (seconds)
    
    AudioSample(double amp = 0.0, double freq = 440.0, double ph = 0.0, double dur = 0.0)
        : amplitude(amp), frequency(freq), phase(ph), duration(dur) {}
};

/**
 * MULTI-LAYERED AUDIO SPECTRUM
 * Represents the complex soundscape of blackhole phenomena
 */
struct BlackHoleAudioSpectrum {
    // Primary frequency layers based on physics
    AudioSample gravitational_wave;     // From orbital motion (20-2000 Hz)
    AudioSample electromagnetic;        // From light emissions (200-8000 Hz)
    AudioSample accretion_disk;        // From temperature variations (50-1000 Hz)  
    AudioSample relativistic_jets;     // From jet activity (100-4000 Hz)
    AudioSample hawking_radiation;     // Theoretical (very high freq, scaled down)
    
    // Environmental effects
    AudioSample doppler_shift;         // Frequency modulation from motion
    AudioSample gravitational_redshift; // Frequency shifting from gravity
    AudioSample frame_dragging;        // Stereo panning from rotation
    
    // Composite output
    double total_amplitude;
    double dominant_frequency;
    double stereo_pan;                 // -1.0 (left) to 1.0 (right)
    
    BlackHoleAudioSpectrum() : total_amplitude(0.0), dominant_frequency(440.0), stereo_pan(0.0) {}
};

/**
 * POSITION-DEPENDENT SONIFICATION PARAMETERS
 * Changes dynamically as observer moves through spacetime
 */
struct SonificationContext {
    Vector4 observer_position;         // Current position in spacetime
    Vector4 observer_velocity;         // Observer's 4-velocity
    double proper_time;               // Observer's proper time
    double coordinate_time;           // Coordinate time
    
    // Distance-dependent scaling
    double schwarzschild_ratio;       // r/r_s (distance in Schwarzschild radii)
    double escape_velocity;           // Local escape velocity
    double tidal_acceleration;        // Tidal force magnitude
    
    // Physics state
    bool inside_ergosphere;           // Frame-dragging region
    bool approaching_horizon;         // Within 10 r_s of event horizon
    bool in_jet_region;              // Within relativistic jet
    bool in_disk_region;             // Within accretion disk
    
    SonificationContext() 
        : proper_time(0.0), coordinate_time(0.0), schwarzschild_ratio(100.0),
          escape_velocity(0.0), tidal_acceleration(0.0),
          inside_ergosphere(false), approaching_horizon(false),
          in_jet_region(false), in_disk_region(false) {}
};

/**
 * MAIN BLACKHOLE SONIFICATION ENGINE
 * Scientifically accurate audio generation from physics simulations
 */
class BlackHoleSonification {
public:
    BlackHoleSonification(const BlackHolePhysics& physics,
                         const SimpleLightSystem& light_sim,
                         const RelativisticJets& jets);
    
    // ====== CORE SONIFICATION FUNCTIONS ======
    
    /**
     * Generate real-time audio spectrum based on current physics state
     * Uses 100% real physics calculations
     */
    BlackHoleAudioSpectrum generateAudioSpectrum(const SonificationContext& context) const;
    
    /**
     * Update sonification context based on observer position
     * Calculates all physics-dependent parameters
     */
    SonificationContext updateContext(const Vector4& observer_pos, 
                                     const Vector4& observer_vel,
                                     double time) const;
    
    /**
     * Generate audio samples for real-time playback
     * Sample rate: 44.1 kHz standard
     */
    std::vector<AudioSample> generateAudioSamples(const BlackHoleAudioSpectrum& spectrum,
                                                 double duration_seconds,
                                                 int sample_rate = 44100) const;
    
    // ====== PHYSICS-TO-FREQUENCY CONVERSIONS ======
    // All based on established physical formulas
    
    /**
     * ORBITAL FREQUENCY TO AUDIO
     * Formula: f_orbital = (1/2π) * √(GM/r³)  [Kepler's third law]
     * Scaled to audible range using logarithmic mapping
     */
    double orbitalFrequencyToAudio(double orbital_radius) const;
    
    /**
     * ELECTROMAGNETIC FREQUENCY TO AUDIO  
     * Formula: f_em = c/λ, scaled from ~10^14 Hz to audible
     * Uses Wien's displacement law: λ_peak = b/T
     */
    double electromagneticFrequencyToAudio(double wavelength_meters) const;
    
    /**
     * TEMPERATURE TO FREQUENCY
     * Formula: f = k_B*T/h  [Thermal energy frequency]
     * Boltzmann constant / Planck constant scaling
     */
    double temperatureToAudioFrequency(double temperature_kelvin) const;
    
    /**
     * SYNCHROTRON FREQUENCY TO AUDIO
     * Formula: f_sync = (eB/2πm_e) * γ²  [Cyclotron frequency with relativistic γ]
     */
    double synchrotronFrequencyToAudio(double electron_energy, double magnetic_field) const;
    
    /**
     * GRAVITATIONAL WAVE FREQUENCY
     * Formula: f_GW = 2*f_orbital  [Quadrupole emission]
     * For binary inspiral: f(t) = (1/8π) * (5GM/c³)^(3/8) * (c³(t_c-t)/5GM)^(-3/8)
     */
    double gravitationalWaveFrequency(double orbital_radius) const;
    
    // ====== RELATIVISTIC EFFECTS ON SOUND ======
    
    /**
     * GRAVITATIONAL REDSHIFT ON FREQUENCY
     * Formula: f_observed = f_emitted * √(-g_tt(r_obs) / -g_tt(r_emit))
     */
    double applyGravitationalRedshift(double base_frequency, 
                                     double emission_radius,
                                     double observer_radius) const;
    
    /**
     * DOPPLER SHIFT ON FREQUENCY  
     * Formula: f_obs = f_emit * γ(1 + β·n̂)  [Relativistic Doppler]
     */
    double applyDopplerShift(double base_frequency,
                            const Vector4& source_velocity,
                            const Vector4& observer_direction) const;
    
    /**
     * FRAME-DRAGGING EFFECT ON STEREO POSITION
     * Formula: Ω = 2Ma/(r³ + a²r + 2Ma²)  [Lense-Thirring precession]
     */
    double calculateFrameDraggingPan(const Vector4& position) const;
    
    // ====== SPECIALIZED PHYSICS AUDIO ======
    
    /**
     * ACCRETION DISK SONIFICATION
     * Temperature profile: T(r) ∝ r^(-3/4)  [Shakura-Sunyaev disk]
     * Orbital velocity: v(r) = √(GM/r)
     */
    AudioSample generateAccretionDiskAudio(const Vector4& position) const;
    
    /**
     * RELATIVISTIC JET SONIFICATION
     * Synchrotron emission from power-law electron distribution
     * Beaming factor: δ = 1/[γ(1-β·n̂)]
     */
    AudioSample generateJetAudio(const Vector4& position, const Vector4& observer) const;
    
    /**
     * HAWKING RADIATION SONIFICATION
     * Temperature: T_H = ℏc³/(8πGMk_B)  [Hawking temperature]
     * Heavily scaled down from ~10^(-7) K for stellar mass black holes
     */
    AudioSample generateHawkingRadiationAudio() const;
    
    /**
     * TIDAL FORCE SONIFICATION
     * Tidal acceleration: a_tidal = 2GM/r³ * Δr  [Differential gravity]
     * Creates rhythmic patterns representing stretching forces
     */
    AudioSample generateTidalForceAudio(const Vector4& position) const;
    
    // ====== DYNAMIC AUDIO EFFECTS ======
    
    /**
     * TIME DILATION EFFECTS ON AUDIO PLAYBACK
     * Formula: dt_proper = dt_coordinate * √(-g_tt)
     * Audio slows down near event horizon due to time dilation
     */
    double calculateTimeDilationFactor(const Vector4& position) const;
    
    /**
     * APPROACH TO EVENT HORIZON EFFECTS
     * Infinite redshift as r → r_s creates dramatic audio effects
     */
    void applyEventHorizonEffects(BlackHoleAudioSpectrum& spectrum, 
                                 double schwarzschild_ratio) const;
    
    /**
     * PLUNGE THROUGH ERGOSPHERE
     * Frame-dragging creates swirling audio effects
     */
    void applyErgosphereEffects(BlackHoleAudioSpectrum& spectrum,
                               const Vector4& position) const;
    
    // ====== CONFIGURATION ======
    
    void setVolumeScale(double scale) { volume_scale_ = scale; }
    void setFrequencyRange(double min_hz, double max_hz) { 
        min_freq_ = min_hz; 
        max_freq_ = max_hz; 
        melodious_calc_.setFrequencyRange(min_hz, max_hz);
    }
    void enableRelativisticEffects(bool enabled) { relativistic_effects_ = enabled; }
    void setSpatialAudio(bool enabled) { spatial_audio_ = enabled; }
    void setDynamicRange(double db) { dynamic_range_db_ = db; }
    
    // ====== MUSICAL CONFIGURATION ======
    
    /**
     * Set the musical scale for sonification
     * Options: CHROMATIC, MAJOR_PENTATONIC, MINOR_PENTATONIC, etc.
     */
    void setMusicalScale(MusicalScale scale) { 
        current_scale_ = scale; 
        melodious_calc_.setDefaultScale(scale);
    }
    
    /**
     * Set frequency mapping method
     * Options: LINEAR, LOGARITHMIC, QUANTIZED, INTERPOLATED, HARMONIC_WEIGHTED
     */
    void setFrequencyMappingMethod(FrequencyMappingMethod method) { mapping_method_ = method; }
    
    /**
     * Enable/disable melodious frequency conversion
     * When enabled, all frequencies are mapped to musical scales
     */
    void enableMelodiousFrequencies(bool enabled) { use_melodious_frequencies_ = enabled; }
    
    /**
     * Set consonance weight for harmonic frequency selection
     * Higher values prefer more consonant (pleasant) intervals
     */
    void setConsonanceWeight(double weight) { melodious_calc_.setConsonanceWeight(weight); }
    
    // ====== AUDIO PROCESSING CONFIGURATION ======
    
    /**
     * Configure high frequency rolloff to reduce harshness
     * @param rolloff_freq Frequency above which to reduce volume (Hz)
     * @param reduction_db Amount of reduction in decibels
     */
    void setHighFrequencyRolloff(double rolloff_freq, double reduction_db) {
        high_freq_rolloff_ = rolloff_freq;
        high_freq_reduction_db_ = reduction_db;
    }
    
    /**
     * Configure dynamic range compression to smooth audio
     * @param threshold Compression threshold (0.0 to 1.0)
     * @param ratio Compression ratio (1.0 = no compression, 4.0 = 4:1 compression)
     * @param attack_ms Attack time in milliseconds
     * @param release_ms Release time in milliseconds
     */
    void setCompression(double threshold, double ratio, double attack_ms, double release_ms) {
        compression_threshold_ = threshold;
        compression_ratio_ = ratio;
        compression_attack_ = attack_ms / 1000.0; // Convert to seconds
        compression_release_ = release_ms / 1000.0;
    }
    
    /**
     * Enable/disable audio processing features
     */
    void enableEQ(bool enabled) { enable_eq_ = enabled; }
    void enableCompression(bool enabled) { enable_compression_ = enabled; }
    void enableDeEssing(bool enabled) { enable_deessing_ = enabled; }
    
    /**
     * Set audio processing to "gentle" preset for pleasant listening
     */
    void setGentleAudioProcessing() {
        setHighFrequencyRolloff(8000.0, -6.0); // Roll off above 8kHz by 6dB
        setCompression(0.7, 3.0, 10.0, 100.0); // Gentle 3:1 compression
        enableEQ(true);
        enableCompression(true);
        enableDeEssing(true);
    }
    
    // Getters
    double getVolumeScale() const { return volume_scale_; }
    bool areRelativisticEffectsEnabled() const { return relativistic_effects_; }
    bool isSpatialAudioEnabled() const { return spatial_audio_; }
    MusicalScale getCurrentMusicalScale() const { return current_scale_; }
    FrequencyMappingMethod getFrequencyMappingMethod() const { return mapping_method_; }
    bool areMelodiousFrequenciesEnabled() const { return use_melodious_frequencies_; }
    
    // ====== DIAGNOSTIC FUNCTIONS ======
    
    /**
     * Get detailed physics information for current audio generation
     */
    std::string getPhysicsBreakdown(const SonificationContext& context) const;
    
    /**
     * Calculate total acoustic power based on physics
     */
    double calculateAcousticPower(const BlackHoleAudioSpectrum& spectrum) const;

private:
    // Physics engine references
    const BlackHolePhysics& physics_;
    const SimpleLightSystem& light_sim_;
    const RelativisticJets& jets_;
    
    // Melodious frequency calculator
    MelodiousCalculator melodious_calc_;
    
    // Configuration parameters
    double volume_scale_;
    double min_freq_;
    double max_freq_;
    bool relativistic_effects_;
    bool spatial_audio_;
    double dynamic_range_db_;
    
    // Musical configuration
    MusicalScale current_scale_;
    FrequencyMappingMethod mapping_method_;
    bool use_melodious_frequencies_;
    
    // Audio processing configuration
    mutable double high_freq_rolloff_; // Hz - frequency above which to reduce volume
    mutable double high_freq_reduction_db_; // dB - how much to reduce high frequencies
    mutable double compression_threshold_; // Threshold for compression
    mutable double compression_ratio_; // Compression ratio
    mutable double compression_attack_; // Attack time for compressor
    mutable double compression_release_; // Release time for compressor
    mutable bool enable_eq_; // Enable frequency EQ
    mutable bool enable_compression_; // Enable dynamic compression
    mutable bool enable_deessing_; // Enable de-essing for harsh frequencies
    
    // Cached physics parameters for performance
    mutable double cached_schwarzschild_radius_;
    mutable double cached_mass_;
    mutable double cached_spin_parameter_;
    mutable bool cache_valid_;
    
    // ====== INTERNAL HELPER FUNCTIONS ======
    
    // Frequency scaling utilities
    double logarithmicScale(double input_freq, double min_in, double max_in,
                           double min_out, double max_out) const;
    double musicalScale(double input_freq, double reference_freq) const;
    
    // Audio synthesis helpers
    double generateSineWave(double frequency, double phase, double time) const;
    double generateSquareWave(double frequency, double phase, double time) const;
    double generateSawtoothWave(double frequency, double phase, double time) const;
    double generateNoiseWave(double amplitude, double time) const;
    
    // Physics calculation helpers
    double calculateOrbitalFrequency(double radius) const;
    double calculateEscapeVelocity(double radius) const;
    double calculateTidalAcceleration(const Vector4& position) const;
    
    // Relativistic calculation helpers
    double lorentzFactor(const Vector4& velocity) const;
    double fourVelocityMagnitude(const Vector4& velocity) const;
    
    // Audio processing
    void normalizeAudioSpectrum(BlackHoleAudioSpectrum& spectrum) const;
    void applyDynamicRangeCompression(BlackHoleAudioSpectrum& spectrum) const;
    void applySpatialEffects(BlackHoleAudioSpectrum& spectrum, const Vector4& position) const;
    
    // Advanced audio processing for harsh frequency reduction
    void applyFrequencyEQ(BlackHoleAudioSpectrum& spectrum) const;
    void applyMultibandCompression(BlackHoleAudioSpectrum& spectrum) const;
    void applyHighFrequencyRolloff(BlackHoleAudioSpectrum& spectrum) const;
    void applyDeEssing(BlackHoleAudioSpectrum& spectrum) const;
    double lowPassFilter(double input, double cutoff_freq, double sample_rate) const;
    double highShelfFilter(double input, double frequency, double gain_db) const;
    double compressor(double input, double threshold, double ratio, double attack, double release) const;
    
    // Cache management
    void updatePhysicsCache() const;
    void invalidateCache() const { cache_valid_ = false; }
    
    // Mathematical utilities
    double clampFrequency(double frequency) const;
    double decibelsToLinear(double db) const;
    double linearToDecibels(double linear) const;
    template<typename T>
    T clamp(T value, T min_val, T max_val) const {
        return std::max(min_val, std::min(max_val, value));
    }
    
    // Missing function declarations
    void calculateCompositeSpectrum(BlackHoleAudioSpectrum& spectrum) const;
    void applyGravitationalRedshiftToSpectrum(BlackHoleAudioSpectrum& spectrum, const SonificationContext& context) const;
    void applyDopplerShiftToSpectrum(BlackHoleAudioSpectrum& spectrum, const SonificationContext& context) const;
    void applyTimeDilationToSpectrum(BlackHoleAudioSpectrum& spectrum, const SonificationContext& context) const;
    AudioSample generateGravitationalWaveAudio(const SonificationContext& context) const;
    AudioSample generateElectromagneticAudio(const SonificationContext& context) const;
    
    // Waveform generators for different physics phenomena
    double accretionDiskWaveform(double time, double radius) const;
    double jetEmissionWaveform(double time, const Vector4& position) const;
    double gravitationalWaveform(double time, double orbital_freq) const;
    double hawkingRadiationWaveform(double time) const;
};

// ====== UTILITY FUNCTIONS FOR SONIFICATION ======

namespace SonificationUtils {
    /**
     * Convert electromagnetic spectrum to musical notes
     * Based on frequency ratios and harmonic series
     */
    std::string frequencyToMusicalNote(double frequency);
    
    /**
     * Calculate beat frequency between two close frequencies
     * f_beat = |f1 - f2|
     */
    double calculateBeatFrequency(double freq1, double freq2);
    
    /**
     * Generate chord progressions based on orbital harmonics
     * Uses Kepler's harmonic law: T² ∝ r³
     */
    std::vector<double> generateOrbitalHarmonics(double fundamental_freq, int num_harmonics);
    
    /**
     * MIDI note conversion for musical analysis
     */
    int frequencyToMIDINote(double frequency);
    double MIDINoteToFrequency(int note);
    
    /**
     * Psychoacoustic scaling for human perception
     * Mel scale: m = 2595 * log₁₀(1 + f/700)
     */
    double frequencyToMelScale(double frequency);
    double melScaleToFrequency(double mel);
    
    /**
     * Bark scale for critical bands
     * z = 13 * arctan(0.00076*f) + 3.5 * arctan((f/7500)²)
     */
    double frequencyToBarkScale(double frequency);
}

/**
 * REAL-TIME AUDIO BUFFER MANAGEMENT
 * For continuous playback of blackhole sonification
 */
class BlackHoleAudioBuffer {
public:
    BlackHoleAudioBuffer(int buffer_size = 4096, int sample_rate = 44100);
    ~BlackHoleAudioBuffer();
    
    // Buffer management
    void fillBuffer(const BlackHoleSonification& sonifier,
                   const SonificationContext& context);
    void getNextSamples(std::vector<double>& left_channel,
                       std::vector<double>& right_channel,
                       int num_samples);
    
    // Playback control
    void startPlayback();
    void stopPlayback();
    bool isPlaying() const { return is_playing_; }
    
    // Status
    bool isBufferReady() const { return buffer_ready_; }
    double getCurrentPlaybackTime() const { return current_time_; }
    
private:
    std::vector<double> left_buffer_;
    std::vector<double> right_buffer_;
    int buffer_size_;
    int sample_rate_;
    int read_position_;
    int write_position_;
    double current_time_;
    bool buffer_ready_;
    bool is_playing_;
    
#ifdef __APPLE__
    // CoreAudio members
    AudioUnit audio_unit_;
    bool audio_initialized_;
    
    // CoreAudio methods
    bool initializeCoreAudio();
    void cleanupCoreAudio();
    static OSStatus audioCallback(void* inRefCon,
                                AudioUnitRenderActionFlags* ioActionFlags,
                                const AudioTimeStamp* inTimeStamp,
                                UInt32 inBusNumber,
                                UInt32 inNumberFrames,
                                AudioBufferList* ioData);
    OSStatus renderCallback(AudioUnitRenderActionFlags* ioActionFlags,
                          const AudioTimeStamp* inTimeStamp,
                          UInt32 inBusNumber,
                          UInt32 inNumberFrames,
                          AudioBufferList* ioData);
#endif
    
    void updateBufferPosition();
};
