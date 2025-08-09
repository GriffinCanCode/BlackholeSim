#pragma once
#include "Physics.h"
#include "LightSystem.h"
#include <array>
#include <vector>


// Forward declarations
class BlackHolePhysics;


// Simple 4x4 matrix for Lorentz transformations
struct Matrix4 {
    std::array<std::array<double, 4>, 4> elements;
    
    Matrix4() {
        // Initialize as identity matrix
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                elements[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }
    
    Vector4 multiply(const Vector4& v) const {
        return Vector4(
            elements[0][0]*v.t + elements[0][1]*v.x + elements[0][2]*v.y + elements[0][3]*v.z,
            elements[1][0]*v.t + elements[1][1]*v.x + elements[1][2]*v.y + elements[1][3]*v.z,
            elements[2][0]*v.t + elements[2][1]*v.x + elements[2][2]*v.y + elements[2][3]*v.z,
            elements[3][0]*v.t + elements[3][1]*v.x + elements[3][2]*v.y + elements[3][3]*v.z
        );
    }
};

// Relativistic jet model parameters
struct RelativisticJetModel {
    Vector4 jet_axis;           // Unit vector along jet axis (usually polar)
    double base_radius;         // Jet radius at launch point (in Schwarzschild radii)
    double jet_height;          // Maximum jet extension
    double bulk_lorentz_factor; // Γ = 1/√(1-β²) for bulk jet motion
    double opening_angle;       // Half-opening angle of jet cone (radians)
    double magnetic_field;      // Characteristic magnetic field strength (Tesla)
    double particle_density;    // Number density of relativistic particles (m⁻³)
    double power_law_index;     // Spectral index for particle energy distribution
    double min_electron_energy; // Minimum electron Lorentz factor
    double max_electron_energy; // Maximum electron Lorentz factor
    bool active;               // Whether jet is currently emitting

    RelativisticJetModel(const Vector4& axis = Vector4(0, 0, 0, 1), double radius = 0.1)
        : jet_axis(VectorMath::normalize3D(axis)), base_radius(radius), 
          jet_height(20.0), bulk_lorentz_factor(10.0), opening_angle(0.1),
          magnetic_field(1e4), particle_density(1e12), power_law_index(2.5),
          min_electron_energy(1.0), max_electron_energy(1e6), active(true) {}
};

// Individual jet emission element for ray tracing
struct JetElement {
    Vector4 position;           // Position in space
    Vector4 velocity;           // 4-velocity of bulk motion
    double density;             // Local particle density
    double magnetic_field;      // Local B-field strength
    double temperature;         // Local electron temperature
    double synchrotron_power;   // Total synchrotron power
    double beaming_factor;      // Relativistic beaming enhancement
    SpectralColor emission;     // Spectral emission profile

    JetElement(const Vector4& pos = Vector4(), const Vector4& vel = Vector4())
        : position(pos), velocity(vel), density(0), magnetic_field(0),
          temperature(0), synchrotron_power(0), beaming_factor(1.0) {}
};

class RelativisticJets {
public:
    RelativisticJets(const BlackHolePhysics& physics);
    
    // Core relativistic beaming calculations
    double calculateBeamingFactor(const Vector4& jet_velocity, const Vector4& observer_direction) const;
    double dopplerBoostingFactor(const Vector4& velocity, const Vector4& photon_direction) const;
    Vector4 calculateJetVelocity(const Vector4& position, const RelativisticJetModel& model) const;
    
    // Jet emission physics
    SpectralColor synchrotronEmission(const JetElement& element, double frequency) const;
    SpectralColor inverseComptonEmission(const JetElement& element, double photon_energy) const;
    double synchrotronFrequency(double electron_energy, double magnetic_field) const;
    
    // Jet geometry and structure
    std::vector<JetElement> generateJetElements(const RelativisticJetModel& model, int resolution = 100) const;
    bool isPointInJet(const Vector4& position, const RelativisticJetModel& model) const;
    double jetDensityProfile(const Vector4& position, const RelativisticJetModel& model) const;
    
    // Relativistic beaming cone calculations  
    double beamingConeAngle(double lorentz_factor) const;
    SpectralColor calculateBeamedEmission(const JetElement& element, const Vector4& observer_position) const;
    double lightAberrationAngle(const Vector4& velocity, const Vector4& emission_direction) const;
    
    // Integration with existing systems
    SpectralColor sampleJetEmission(const Vector4& ray_position, const Vector4& ray_direction,
                                   const Vector4& observer_position, const RelativisticJetModel& model) const;
    void updateJetElements(std::vector<JetElement>& elements, double time_step) const;
    
    // Multi-wavelength emission
    SpectralColor calculateSpectralEnergyDistribution(const JetElement& element, 
                                                     const std::vector<double>& frequencies) const;
    double fluxDensity(const JetElement& element, double frequency, double distance) const;
    
    // Advanced relativistic effects
    SpectralColor applyRelativisticCorrections(const SpectralColor& intrinsic_emission,
                                              const Vector4& bulk_velocity, 
                                              const Vector4& observer_direction) const;
    Vector4 transformToLabFrame(const Vector4& jet_frame_vector, const Vector4& bulk_velocity) const;
    Vector4 transformToJetFrame(const Vector4& lab_frame_vector, const Vector4& bulk_velocity) const;
    
    // Configuration and settings
    void setJetModel(const RelativisticJetModel& model) { jet_model_ = model; }
    const RelativisticJetModel& getJetModel() const { return jet_model_; }
    void setSpectralResolution(int resolution) { spectral_resolution_ = resolution; }
    void enableVariability(bool enabled) { variability_enabled_ = enabled; }
    void setTimescale(double timescale) { variability_timescale_ = timescale; }
    
    // Diagnostic and visualization
    std::vector<Vector4> getJetStreamlines(const RelativisticJetModel& model, int num_lines = 10) const;
    double totalJetLuminosity(const RelativisticJetModel& model) const;
    SpectralColor integratedJetSpectrum(const RelativisticJetModel& model, 
                                       const Vector4& observer_position) const;

private:
    const BlackHolePhysics& physics_;

    RelativisticJetModel jet_model_;
    int spectral_resolution_;
    bool variability_enabled_;
    double variability_timescale_;
    
public:
    // Physical constants for jet physics
    static constexpr double ELECTRON_CHARGE = 1.602176634e-19;     // C
    static constexpr double ELECTRON_MASS = 9.1093837015e-31;      // kg
    static constexpr double THOMSON_CROSS_SECTION = 6.6524587e-29; // m²
    static constexpr double CLASSICAL_ELECTRON_RADIUS = 2.8179403262e-15; // m

private:
    
    // Internal calculation helpers
    double calculateSynchrotronPower(double electron_energy, double magnetic_field) const;
    double electronEnergyDistribution(double energy, const RelativisticJetModel& model) const;
    SpectralColor integrateEmissionOverElectrons(const JetElement& element, double frequency) const;
    
    // Relativistic transformation utilities
    double lorentzFactor(const Vector4& velocity) const;
    Vector4 fourVelocity(const Vector4& three_velocity) const;
    Matrix4 lorentzBoostMatrix(const Vector4& velocity) const;
    
    // Jet dynamics
    Vector4 calculateJetAcceleration(const Vector4& position, const RelativisticJetModel& model) const;
    double jetPressure(const Vector4& position, const RelativisticJetModel& model) const;
    Vector4 magneticFieldVector(const Vector4& position, const RelativisticJetModel& model) const;
    
    // Cooling and energy loss
    double synchrotronCoolingTime(double electron_energy, double magnetic_field) const;
    double inverseComptonCoolingTime(double electron_energy, double photon_density) const;
    void applyCooling(JetElement& element, double time_step) const;
    
    // Variability and turbulence
    double calculateVariabilityFactor(const Vector4& position, double time) const;
    SpectralColor addTurbulentFluctuations(const SpectralColor& base_emission, 
                                          const Vector4& position, double time) const;
};

// Utility functions for jet physics
namespace JetPhysics {
    // Critical frequencies
    double synchrotronSelfAbsorptionFrequency(double magnetic_field, double density);
    double plasmaCutoffFrequency(double density);
    double characteristicSynchrotronFrequency(double electron_energy, double magnetic_field);
    
    // Emission coefficients
    double synchrotronEmissivity(double frequency, double electron_density, 
                                double magnetic_field, double power_law_index);
    double synchrotronAbsorptionCoefficient(double frequency, double electron_density,
                                          double magnetic_field, double power_law_index);
    
    // Relativistic plasma physics
    double plasmaDispersionRelation(double frequency, double density, double magnetic_field);
    Vector4 calculateJetPolarization(const Vector4& magnetic_field, const Vector4& propagation_direction);
}
