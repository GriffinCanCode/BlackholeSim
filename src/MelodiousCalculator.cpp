#include "MelodiousCalculator.h"
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <iostream>

using namespace std;

MelodiousCalculator::MelodiousCalculator(MusicalScale default_scale)
    : default_scale_(default_scale)
    , min_frequency_(MusicalConstants::MIN_PLEASANT_FREQ)
    , max_frequency_(MusicalConstants::MAX_PLEASANT_FREQ)
    , consonance_weight_(1.0) {
}

// ====== CORE MAPPING FUNCTIONS ======

MusicalNote MelodiousCalculator::mapFrequencyToScale(double input_freq, 
                                                    MusicalScale scale,
                                                    FrequencyMappingMethod method) const {
    MusicalNote result;
    
    // Clamp input to reasonable range
    input_freq = clamp(input_freq, MusicalConstants::C0_FREQ, MusicalConstants::C8_FREQ);
    
    // Get available scale frequencies in our range
    vector<double> scale_frequencies = generateScaleFrequencies(scale, min_frequency_, max_frequency_);
    
    double mapped_freq;
    
    switch (method) {
        case FrequencyMappingMethod::QUANTIZED: {
            // Find closest scale note
            mapped_freq = quantizeToScale(input_freq, scale_frequencies);
            break;
        }
        case FrequencyMappingMethod::INTERPOLATED: {
            // Smooth interpolation between scale notes
            auto it = lower_bound(scale_frequencies.begin(), scale_frequencies.end(), input_freq);
            if (it == scale_frequencies.begin()) {
                mapped_freq = scale_frequencies[0];
            } else if (it == scale_frequencies.end()) {
                mapped_freq = scale_frequencies.back();
            } else {
                double freq_high = *it;
                double freq_low = *(it - 1);
                double ratio = (input_freq - freq_low) / (freq_high - freq_low);
                // Use musical interpolation (geometric mean)
                mapped_freq = freq_low * pow(freq_high / freq_low, ratio);
            }
            break;
        }
        case FrequencyMappingMethod::HARMONIC_WEIGHTED: {
            // Weight by harmonic consonance
            mapped_freq = findMostConsonantFrequency(input_freq, scale_frequencies);
            break;
        }
        default: {
            // Default to quantized
            mapped_freq = quantizeToScale(input_freq, scale_frequencies);
            break;
        }
    }
    
    // Fill in musical note information
    result.frequency = mapped_freq;
    result.midi_note_number = frequencyToMIDI(mapped_freq);
    result.note_name = frequencyToNoteName(mapped_freq);
    
    // Calculate semitone offset from A4
    result.semitone_offset = frequencyToSemitones(mapped_freq);
    result.octave = (result.midi_note_number - 12) / 12;
    
    // Calculate consonance score
    result.consonance_score = calculatePleasantnessScore(mapped_freq);
    
    return result;
}

double MelodiousCalculator::mapNormalizedToMusicalFreq(double normalized_value,
                                                      MusicalScale scale,
                                                      double base_octave) const {
    // Clamp input
    normalized_value = clamp(normalized_value, 0.0, 1.0);
    
    // Get scale frequencies around the base octave
    double octave_min = MusicalConstants::CONCERT_A4 * pow(2.0, base_octave - 4.0 - 1.0); // One octave below
    double octave_max = MusicalConstants::CONCERT_A4 * pow(2.0, base_octave - 4.0 + 1.0); // One octave above
    
    vector<double> scale_frequencies = generateScaleFrequencies(scale, octave_min, octave_max);
    
    // Map normalized value to scale index
    int scale_index = static_cast<int>(normalized_value * (scale_frequencies.size() - 1));
    scale_index = clamp(scale_index, 0, static_cast<int>(scale_frequencies.size() - 1));
    
    return scale_frequencies[scale_index];
}

double MelodiousCalculator::convertPhysicsToAudio(double physics_freq,
                                                 double physics_min, double physics_max,
                                                 MusicalScale scale) const {
    // Logarithmic mapping from physics range to musical range
    double log_physics = log10(clamp(physics_freq, physics_min, physics_max));
    double log_min = log10(physics_min);
    double log_max = log10(physics_max);
    
    // Normalize to 0-1
    double normalized = (log_physics - log_min) / (log_max - log_min);
    
    // Map to musical frequency using scale
    return mapNormalizedToMusicalFreq(normalized, scale, 4.0);
}

// ====== SCALE GENERATION FUNCTIONS ======

vector<double> MelodiousCalculator::generateScaleFrequencies(MusicalScale scale,
                                                           double min_freq,
                                                           double max_freq) const {
    vector<double> frequencies;
    ScaleDefinition scale_def = getScaleDefinition(scale);
    
    // Find starting octave
    double base_freq = MusicalConstants::CONCERT_A4; // A4 = 440 Hz
    
    // Generate frequencies for multiple octaves
    for (int octave = -2; octave <= 3; ++octave) {
        double octave_base = base_freq * pow(2.0, octave);
        
        // Generate all scale notes in this octave
        for (int interval : scale_def.semitone_intervals) {
            double freq = octave_base * pow(MusicalConstants::SEMITONE_RATIO, interval);
            if (freq >= min_freq && freq <= max_freq) {
                frequencies.push_back(freq);
            }
        }
    }
    
    // Sort and remove duplicates
    sort(frequencies.begin(), frequencies.end());
    frequencies.erase(unique(frequencies.begin(), frequencies.end()), frequencies.end());
    
    return frequencies;
}

ScaleDefinition MelodiousCalculator::getScaleDefinition(MusicalScale scale) const {
    // Check cache first
    if (scale_cache_.find(scale) != scale_cache_.end()) {
        return scale_cache_[scale];
    }
    
    ScaleDefinition def = buildScaleDefinition(scale);
    scale_cache_[scale] = def;
    return def;
}

vector<double> MelodiousCalculator::generateHarmonicSeries(double fundamental, int num_harmonics) const {
    vector<double> harmonics;
    harmonics.reserve(num_harmonics);
    
    for (int i = 1; i <= num_harmonics; ++i) {
        double harmonic_freq = fundamental * i;
        if (harmonic_freq <= max_frequency_) {
            harmonics.push_back(harmonic_freq);
        }
    }
    
    return harmonics;
}

// ====== MUSICAL ANALYSIS FUNCTIONS ======

double MelodiousCalculator::calculateConsonanceScore(double freq1, double freq2) const {
    double ratio = max(freq1, freq2) / min(freq1, freq2);
    
    // Consonant ratios and their scores
    vector<pair<double, double>> consonant_ratios = {
        {1.0, 1.0},      // Unison
        {2.0, 0.95},     // Octave
        {1.5, 0.9},      // Perfect fifth
        {1.333, 0.85},   // Perfect fourth
        {1.25, 0.8},     // Major third
        {1.2, 0.75},     // Minor third
        {1.125, 0.7},    // Major second
        {1.067, 0.65}    // Minor second
    };
    
    double best_score = 0.0;
    for (const auto& consonant_pair : consonant_ratios) {
        double target_ratio = consonant_pair.first;
        double max_score = consonant_pair.second;
        
        // Check both the ratio and its reciprocal
        double error1 = abs(ratio - target_ratio) / target_ratio;
        double error2 = abs(ratio - (1.0 / target_ratio)) / (1.0 / target_ratio);
        double min_error = min(error1, error2);
        
        if (min_error < 0.05) { // Within 5% tolerance
            double score = max_score * (1.0 - min_error / 0.05);
            best_score = max(best_score, score);
        }
    }
    
    return best_score;
}

double MelodiousCalculator::findMostConsonantFrequency(double input_freq,
                                                      const vector<double>& reference_frequencies) const {
    if (reference_frequencies.empty()) return input_freq;
    
    double best_freq = reference_frequencies[0];
    double best_score = 0.0;
    
    for (double ref_freq : reference_frequencies) {
        double consonance = calculateConsonanceScore(input_freq, ref_freq);
        double distance_penalty = abs(log(ref_freq / input_freq)); // Logarithmic distance
        double total_score = consonance - 0.1 * distance_penalty; // Weight consonance over proximity
        
        if (total_score > best_score) {
            best_score = total_score;
            best_freq = ref_freq;
        }
    }
    
    return best_freq;
}

double MelodiousCalculator::calculatePleasantnessScore(double frequency) const {
    // Psychoacoustic pleasantness based on:
    // 1. Frequency range (middle frequencies are most pleasant)
    // 2. Harmonic content
    // 3. Avoiding harsh frequencies
    
    double range_score = 1.0;
    
    // Optimal range around 200-800 Hz (human voice range)
    if (frequency >= 200.0 && frequency <= 800.0) {
        range_score = 1.0;
    } else if (frequency >= 100.0 && frequency <= 1600.0) {
        range_score = 0.8;
    } else if (frequency >= 50.0 && frequency <= 3200.0) {
        range_score = 0.6;
    } else {
        range_score = 0.3;
    }
    
    // Penalize frequencies that are too low or too high
    if (frequency < 40.0) range_score *= 0.1;
    if (frequency > 4000.0) range_score *= 0.2;
    
    // Bonus for frequencies close to harmonic series of A4
    double harmonic_bonus = 0.0;
    for (int harmonic = 1; harmonic <= 16; ++harmonic) {
        double harmonic_freq = MusicalConstants::CONCERT_A4 * harmonic;
        if (abs(frequency - harmonic_freq) / harmonic_freq < 0.02) {
            harmonic_bonus = 0.2;
            break;
        }
    }
    
    return clamp(range_score + harmonic_bonus, 0.0, 1.0);
}

// ====== UTILITY FUNCTIONS ======

int MelodiousCalculator::frequencyToMIDI(double frequency) const {
    if (frequency <= 0) return 0;
    
    // MIDI formula: note = 69 + 12 * log2(frequency / 440)
    double midi_float = 69.0 + 12.0 * log2(frequency / MusicalConstants::CONCERT_A4);
    return static_cast<int>(round(midi_float));
}

double MelodiousCalculator::MIDIToFrequency(int midi_note) const {
    // Frequency = 440 * 2^((MIDI - 69) / 12)
    return MusicalConstants::CONCERT_A4 * pow(2.0, (midi_note - 69) / 12.0);
}

string MelodiousCalculator::frequencyToNoteName(double frequency) const {
    const vector<string> note_names = {
        "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"
    };
    
    int midi_note = frequencyToMIDI(frequency);
    int octave = (midi_note - 12) / 12;
    int note_index = (midi_note - 12) % 12;
    
    return note_names[note_index] + to_string(octave);
}

double MelodiousCalculator::getIntervalRatio(const string& interval_name) const {
    map<string, double> intervals = {
        {"unison", 1.0},
        {"minor_second", 16.0/15.0},
        {"major_second", 9.0/8.0},
        {"minor_third", 6.0/5.0},
        {"major_third", 5.0/4.0},
        {"perfect_fourth", 4.0/3.0},
        {"tritone", sqrt(2.0)},
        {"perfect_fifth", 3.0/2.0},
        {"minor_sixth", 8.0/5.0},
        {"major_sixth", 5.0/3.0},
        {"minor_seventh", 16.0/9.0},
        {"major_seventh", 15.0/8.0},
        {"octave", 2.0}
    };
    
    auto it = intervals.find(interval_name);
    return (it != intervals.end()) ? it->second : 1.0;
}

// ====== INTERNAL HELPER FUNCTIONS ======

ScaleDefinition MelodiousCalculator::buildScaleDefinition(MusicalScale scale) const {
    ScaleDefinition def;
    
    switch (scale) {
        case MusicalScale::CHROMATIC:
            def.semitone_intervals = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
            def.name = "Chromatic";
            def.consonance_rating = 0.5;
            def.description = "All twelve semitones - maximum flexibility";
            break;
            
        case MusicalScale::MAJOR_PENTATONIC:
            def.semitone_intervals = {0, 2, 4, 7, 9}; // C, D, E, G, A
            def.name = "Major Pentatonic";
            def.consonance_rating = 0.95;
            def.description = "Very pleasant, no dissonant intervals";
            break;
            
        case MusicalScale::MINOR_PENTATONIC:
            def.semitone_intervals = {0, 3, 5, 7, 10}; // C, Eb, F, G, Bb
            def.name = "Minor Pentatonic";
            def.consonance_rating = 0.9;
            def.description = "Bluesy and expressive";
            break;
            
        case MusicalScale::MAJOR_SCALE:
            def.semitone_intervals = {0, 2, 4, 5, 7, 9, 11}; // C, D, E, F, G, A, B
            def.name = "Major Scale";
            def.consonance_rating = 0.8;
            def.description = "Classical major scale - bright and happy";
            break;
            
        case MusicalScale::MINOR_SCALE:
            def.semitone_intervals = {0, 2, 3, 5, 7, 8, 10}; // C, D, Eb, F, G, Ab, Bb
            def.name = "Natural Minor";
            def.consonance_rating = 0.8;
            def.description = "Classical minor scale - melancholic";
            break;
            
        case MusicalScale::DORIAN_MODE:
            def.semitone_intervals = {0, 2, 3, 5, 7, 9, 10}; // C, D, Eb, F, G, A, Bb
            def.name = "Dorian Mode";
            def.consonance_rating = 0.85;
            def.description = "Jazz/modal sound - sophisticated";
            break;
            
        case MusicalScale::WHOLE_TONE:
            def.semitone_intervals = {0, 2, 4, 6, 8, 10}; // C, D, E, F#, G#, A#
            def.name = "Whole Tone";
            def.consonance_rating = 0.7;
            def.description = "Dreamy, ethereal - impressionist";
            break;
            
        case MusicalScale::BLUES_SCALE:
            def.semitone_intervals = {0, 3, 5, 6, 7, 10}; // C, Eb, F, F#, G, Bb
            def.name = "Blues Scale";
            def.consonance_rating = 0.75;
            def.description = "Expressive with blue notes";
            break;
            
        case MusicalScale::HARMONIC_SERIES:
            // Natural harmonics - most consonant
            def.semitone_intervals = {0, 12, 19, 24, 28, 31, 34}; // Approximate semitones for harmonics
            def.name = "Harmonic Series";
            def.consonance_rating = 1.0;
            def.description = "Natural harmonics - pure consonance";
            break;
            
        default:
            def = buildScaleDefinition(MusicalScale::MAJOR_PENTATONIC);
            break;
    }
    
    return def;
}

double MelodiousCalculator::calculateSemitoneFrequency(int semitones_from_a4) const {
    return MusicalConstants::CONCERT_A4 * pow(MusicalConstants::SEMITONE_RATIO, semitones_from_a4);
}

int MelodiousCalculator::frequencyToSemitones(double frequency) const {
    if (frequency <= 0) return 0;
    return static_cast<int>(round(12.0 * log2(frequency / MusicalConstants::CONCERT_A4)));
}

double MelodiousCalculator::quantizeToScale(double frequency, const vector<double>& scale_frequencies) const {
    if (scale_frequencies.empty()) return frequency;
    
    return findClosestValue(frequency, scale_frequencies);
}

double MelodiousCalculator::findClosestValue(double target, const vector<double>& values) const {
    if (values.empty()) return target;
    
    auto it = lower_bound(values.begin(), values.end(), target);
    
    if (it == values.begin()) return *it;
    if (it == values.end()) return values.back();
    
    // Check which is closer: the lower or higher value
    double higher = *it;
    double lower = *(it - 1);
    
    // Use logarithmic distance for musical perception
    double dist_to_lower = abs(log(target / lower));
    double dist_to_higher = abs(log(higher / target));
    
    return (dist_to_lower < dist_to_higher) ? lower : higher;
}

string MelodiousCalculator::getFrequencyAnalysis(double input_freq, MusicalScale scale) const {
    MusicalNote note = mapFrequencyToScale(input_freq, scale);
    ScaleDefinition scale_def = getScaleDefinition(scale);
    
    ostringstream oss;
    oss << fixed << setprecision(2);
    oss << "=== FREQUENCY ANALYSIS ===" << endl;
    oss << "Input frequency: " << input_freq << " Hz" << endl;
    oss << "Mapped to: " << note.note_name << " (" << note.frequency << " Hz)" << endl;
    oss << "Scale: " << scale_def.name << endl;
    oss << "Consonance rating: " << scale_def.consonance_rating << endl;
    oss << "Pleasantness score: " << note.consonance_score << endl;
    oss << "MIDI note: " << note.midi_note_number << endl;
    
    return oss.str();
}

vector<string> MelodiousCalculator::getAvailableScales() const {
    return {
        "Chromatic (all notes)",
        "Major Pentatonic (very pleasant)",
        "Minor Pentatonic (bluesy)",
        "Major Scale (classical bright)",
        "Minor Scale (classical sad)",
        "Dorian Mode (jazzy)",
        "Whole Tone (dreamy)",
        "Blues Scale (expressive)",
        "Harmonic Series (pure consonance)"
    };
}

// ====== UTILITY NAMESPACE IMPLEMENTATIONS ======

namespace MelodiousUtils {
    string hzToNoteName(double frequency) {
        MelodiousCalculator calc;
        return calc.frequencyToNoteName(frequency);
    }
    
    double noteNameToHz(const string& note_name) {
        // Simple parser for note names like "A4", "C#5", etc.
        // This is a simplified implementation
        if (note_name == "A4") return 440.0;
        if (note_name == "C4") return 261.63;
        if (note_name == "C5") return 523.25;
        // Add more as needed...
        return 440.0; // Default to A4
    }
    
    vector<double> createConsonantChord(double fundamental, const string& chord_type) {
        vector<double> chord;
        chord.push_back(fundamental);
        
        if (chord_type == "major") {
            chord.push_back(fundamental * 1.25);  // Major third
            chord.push_back(fundamental * 1.5);   // Perfect fifth
        } else if (chord_type == "minor") {
            chord.push_back(fundamental * 1.2);   // Minor third
            chord.push_back(fundamental * 1.5);   // Perfect fifth
        } else if (chord_type == "perfect") {
            chord.push_back(fundamental * 1.333); // Perfect fourth
            chord.push_back(fundamental * 1.5);   // Perfect fifth
            chord.push_back(fundamental * 2.0);   // Octave
        }
        
        return chord;
    }
    
    vector<double> createMusicalSequence(const vector<double>& data, MusicalScale scale) {
        MelodiousCalculator calc(scale);
        vector<double> musical_sequence;
        
        if (data.empty()) return musical_sequence;
        
        // Find data range
        auto minmax = minmax_element(data.begin(), data.end());
        double data_min = *minmax.first;
        double data_max = *minmax.second;
        
        for (double value : data) {
            double normalized = (data_max != data_min) ? 
                               (value - data_min) / (data_max - data_min) : 0.5;
            double musical_freq = calc.mapNormalizedToMusicalFreq(normalized, scale);
            musical_sequence.push_back(musical_freq);
        }
        
        return musical_sequence;
    }
    
    double calculateBeatFrequency(double freq1, double freq2) {
        return abs(freq1 - freq2);
    }
    
    bool isInPleasantRange(double frequency) {
        return frequency >= MusicalConstants::MIN_PLEASANT_FREQ && 
               frequency <= MusicalConstants::MAX_PLEASANT_FREQ;
    }
    
    pair<double, double> getRecommendedRange(const string& phenomenon_type) {
        if (phenomenon_type == "gravitational_waves") {
            return {55.0, 220.0};   // Low bass range
        } else if (phenomenon_type == "electromagnetic") {
            return {220.0, 880.0};  // Mid range
        } else if (phenomenon_type == "accretion_disk") {
            return {110.0, 440.0};  // Low-mid range
        } else if (phenomenon_type == "jets") {
            return {440.0, 1760.0}; // High range
        } else {
            return {MusicalConstants::MIN_PLEASANT_FREQ, MusicalConstants::MAX_PLEASANT_FREQ};
        }
    }
}
