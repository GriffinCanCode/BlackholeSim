#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <map>

/**
 * MELODIOUS FREQUENCY CALCULATOR
 * Converts arbitrary physical frequencies to musically pleasant scales
 * Based on equal temperament, chromatic scales, and consonant intervals
 */

namespace MusicalConstants {
    // Standard musical reference frequencies
    constexpr double CONCERT_A4 = 440.0;              // A4 reference pitch (Hz)
    constexpr double SEMITONE_RATIO = 1.059463094359; // 12th root of 2
    constexpr double OCTAVE_RATIO = 2.0;               // Frequency doubles per octave
    
    // Frequency ranges for different octaves
    constexpr double C0_FREQ = 16.35;                  // C0 (lowest piano note region)
    constexpr double C8_FREQ = 4186.01;               // C8 (highest piano note region)
    
    // Pleasant frequency ranges for sonification
    constexpr double MIN_PLEASANT_FREQ = 55.0;        // A1 (low but audible)
    constexpr double MAX_PLEASANT_FREQ = 1760.0;      // A6 (high but not piercing)
}

/**
 * MUSICAL SCALE TYPES
 * Different scales for different moods and consonance levels
 */
enum class MusicalScale {
    CHROMATIC,          // All 12 semitones (most flexible, can be dissonant)
    MAJOR_PENTATONIC,   // 5 notes (C, D, E, G, A) - very pleasant, no dissonance
    MINOR_PENTATONIC,   // 5 notes (C, Eb, F, G, Bb) - bluesy, pleasant
    MAJOR_SCALE,        // 7 notes (Do, Re, Mi, Fa, Sol, La, Ti) - classical
    MINOR_SCALE,        // 7 notes (natural minor) - melancholic but pleasant
    DORIAN_MODE,        // 7 notes (jazz/modal) - sophisticated
    WHOLE_TONE,         // 6 notes (dreamy, ethereal) - impressionist
    BLUES_SCALE,        // 6 notes (minor pentatonic + blue note) - expressive
    HARMONIC_SERIES     // Natural harmonics (most consonant) - pure physics
};

/**
 * FREQUENCY MAPPING METHODS
 * Different ways to map input data to musical scales
 */
enum class FrequencyMappingMethod {
    LINEAR,             // Direct linear mapping to scale
    LOGARITHMIC,        // Log mapping (better for wide dynamic ranges)
    EXPONENTIAL,        // Exponential mapping (emphasizes high values)
    QUANTIZED,          // Snap to nearest scale note (most musical)
    INTERPOLATED,       // Smooth interpolation between scale notes
    HARMONIC_WEIGHTED   // Weight by harmonic consonance
};

/**
 * MUSICAL NOTE REPRESENTATION
 * Contains both frequency and musical information
 */
struct MusicalNote {
    double frequency;           // Frequency in Hz
    int midi_note_number;      // MIDI note number (C4 = 60)
    std::string note_name;     // Human readable name (e.g., "A4", "C#5")
    int octave;                // Octave number
    int semitone_offset;       // Semitones from A4
    double consonance_score;   // How consonant this note is (0.0 to 1.0)
    
    MusicalNote(double freq = 440.0) 
        : frequency(freq), midi_note_number(69), note_name("A4"), 
          octave(4), semitone_offset(0), consonance_score(1.0) {}
};

/**
 * SCALE DEFINITION
 * Defines the intervals and characteristics of a musical scale
 */
struct ScaleDefinition {
    std::vector<int> semitone_intervals;    // Semitone intervals from root
    std::string name;                       // Scale name
    double consonance_rating;               // Overall consonance (0.0 to 1.0)
    std::string description;                // Musical description
    
    ScaleDefinition() : consonance_rating(0.5) {}
};

/**
 * MAIN MELODIOUS CALCULATOR CLASS
 * Converts physics-based frequencies to musically pleasant frequencies
 */
class MelodiousCalculator {
public:
    MelodiousCalculator(MusicalScale default_scale = MusicalScale::MAJOR_PENTATONIC);
    
    // ====== CORE MAPPING FUNCTIONS ======
    
    /**
     * Convert arbitrary frequency to nearest musical note
     * Returns the closest note in the selected scale
     */
    MusicalNote mapFrequencyToScale(double input_freq, 
                                   MusicalScale scale = MusicalScale::MAJOR_PENTATONIC,
                                   FrequencyMappingMethod method = FrequencyMappingMethod::QUANTIZED) const;
    
    /**
     * Map a range of frequencies to a musical scale
     * Useful for mapping physics data ranges to pleasant frequency ranges
     */
    std::vector<MusicalNote> mapFrequencyRange(double min_freq, double max_freq,
                                              int num_points,
                                              MusicalScale scale = MusicalScale::MAJOR_PENTATONIC,
                                              FrequencyMappingMethod method = FrequencyMappingMethod::LINEAR) const;
    
    /**
     * Map normalized data value (0.0 to 1.0) to musical frequency
     * Most useful for sonification of physics parameters
     */
    double mapNormalizedToMusicalFreq(double normalized_value,
                                     MusicalScale scale = MusicalScale::MAJOR_PENTATONIC,
                                     double base_octave = 4.0) const;
    
    /**
     * Convert physics frequency to pleasant audio frequency
     * Uses logarithmic mapping with musical quantization
     */
    double convertPhysicsToAudio(double physics_freq,
                                double physics_min, double physics_max,
                                MusicalScale scale = MusicalScale::MAJOR_PENTATONIC) const;
    
    // ====== SCALE GENERATION FUNCTIONS ======
    
    /**
     * Generate all frequencies for a given scale in a frequency range
     * Returns all available notes in the scale within the range
     */
    std::vector<double> generateScaleFrequencies(MusicalScale scale,
                                                double min_freq = MusicalConstants::MIN_PLEASANT_FREQ,
                                                double max_freq = MusicalConstants::MAX_PLEASANT_FREQ) const;
    
    /**
     * Get scale definition for a musical scale
     * Returns the semitone intervals and characteristics
     */
    ScaleDefinition getScaleDefinition(MusicalScale scale) const;
    
    /**
     * Calculate harmonic series based on fundamental frequency
     * Returns natural harmonic frequencies (most consonant)
     */
    std::vector<double> generateHarmonicSeries(double fundamental, int num_harmonics = 8) const;
    
    // ====== MUSICAL ANALYSIS FUNCTIONS ======
    
    /**
     * Calculate consonance score between two frequencies
     * Higher score = more consonant (pleasant sounding together)
     */
    double calculateConsonanceScore(double freq1, double freq2) const;
    
    /**
     * Find the most consonant frequency near the input frequency
     * Uses harmonic ratios and interval theory
     */
    double findMostConsonantFrequency(double input_freq, 
                                     const std::vector<double>& reference_frequencies) const;
    
    /**
     * Calculate the "pleasantness" score of a frequency
     * Based on psychoacoustic principles and musical theory
     */
    double calculatePleasantnessScore(double frequency) const;
    
    // ====== UTILITY FUNCTIONS ======
    
    /**
     * Convert frequency to MIDI note number
     * A4 (440 Hz) = MIDI note 69
     */
    int frequencyToMIDI(double frequency) const;
    
    /**
     * Convert MIDI note number to frequency
     * MIDI note 69 = A4 (440 Hz)
     */
    double MIDIToFrequency(int midi_note) const;
    
    /**
     * Get musical note name from frequency
     * Returns string like "A4", "C#5", etc.
     */
    std::string frequencyToNoteName(double frequency) const;
    
    /**
     * Calculate frequency ratio for musical interval
     * Returns ratio for perfect consonant intervals
     */
    double getIntervalRatio(const std::string& interval_name) const;
    
    // ====== CONFIGURATION ======
    
    void setDefaultScale(MusicalScale scale) { default_scale_ = scale; }
    void setFrequencyRange(double min_freq, double max_freq) {
        min_frequency_ = min_freq;
        max_frequency_ = max_freq;
    }
    void setConsonanceWeight(double weight) { consonance_weight_ = weight; }
    
    // Getters
    MusicalScale getDefaultScale() const { return default_scale_; }
    std::pair<double, double> getFrequencyRange() const { 
        return {min_frequency_, max_frequency_}; 
    }
    
    // ====== DIAGNOSTIC FUNCTIONS ======
    
    /**
     * Get detailed analysis of a frequency mapping
     */
    std::string getFrequencyAnalysis(double input_freq, MusicalScale scale) const;
    
    /**
     * Get information about available scales
     */
    std::vector<std::string> getAvailableScales() const;

private:
    MusicalScale default_scale_;
    double min_frequency_;
    double max_frequency_;
    double consonance_weight_;
    
    // Cached scale definitions for performance
    mutable std::map<MusicalScale, ScaleDefinition> scale_cache_;
    
    // ====== INTERNAL HELPER FUNCTIONS ======
    
    // Scale generation helpers
    ScaleDefinition buildScaleDefinition(MusicalScale scale) const;
    std::vector<int> getScaleIntervals(MusicalScale scale) const;
    
    // Frequency calculation helpers
    double calculateSemitoneFrequency(int semitones_from_a4) const;
    int frequencyToSemitones(double frequency) const;
    double quantizeToScale(double frequency, const std::vector<double>& scale_frequencies) const;
    
    // Musical theory helpers
    double calculateFrequencyRatio(double freq1, double freq2) const;
    bool isConsonantInterval(double ratio) const;
    double getConsonanceWeight(double ratio) const;
    
    // Mapping algorithm implementations
    double linearMapping(double input, double input_min, double input_max,
                        double output_min, double output_max) const;
    double logarithmicMapping(double input, double input_min, double input_max,
                             double output_min, double output_max) const;
    double exponentialMapping(double input, double input_min, double input_max,
                             double output_min, double output_max) const;
    
    // Utility functions
    template<typename T>
    T clamp(T value, T min_val, T max_val) const {
        return std::max(min_val, std::min(max_val, value));
    }
    
    double findClosestValue(double target, const std::vector<double>& values) const;
    int findClosestIndex(double target, const std::vector<double>& values) const;
};

// ====== UTILITY NAMESPACE ======

namespace MelodiousUtils {
    /**
     * Quick conversion functions for common use cases
     */
    
    // Convert Hz to musical note name
    std::string hzToNoteName(double frequency);
    
    // Convert note name to Hz (e.g., "A4" -> 440.0)
    double noteNameToHz(const std::string& note_name);
    
    // Create pleasant chord from fundamental frequency
    std::vector<double> createConsonantChord(double fundamental, 
                                           const std::string& chord_type = "major");
    
    // Generate pleasant frequency sequence for time-varying data
    std::vector<double> createMusicalSequence(const std::vector<double>& data,
                                             MusicalScale scale = MusicalScale::MAJOR_PENTATONIC);
    
    // Calculate beat frequency for two notes (pleasant for modulation)
    double calculateBeatFrequency(double freq1, double freq2);
    
    // Check if frequency is in pleasant hearing range
    bool isInPleasantRange(double frequency);
    
    // Get recommended frequency range for different physics phenomena
    std::pair<double, double> getRecommendedRange(const std::string& phenomenon_type);
}
