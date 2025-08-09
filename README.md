# Black Hole Simulator

A physically accurate 3D black hole visualization and sonification system using C++ and OpenGL, implementing real general relativity physics with scientifically accurate audio generation.

## Features

### Physics Accuracy
- **Dual Physics Engines**: Factory pattern supporting both Schwarzschild (non-rotating) and Kerr (rotating) black holes
- **Schwarzschild Metric**: Accurate spacetime curvature around non-rotating black holes
- **Kerr Metric**: Full rotating black hole physics with frame dragging and ergosphere effects
- **Ergosphere Physics**: Regions where spacetime itself is dragged, enabling energy extraction
- **Dual Event Horizons**: Outer and inner horizons for rotating black holes with spin-dependent geometry
- **Geodesic Ray Tracing**: Light paths calculated using Einstein's field equations
- **4th-Order Runge-Kutta Integration**: Numerical solution of geodesic differential equations
- **Adaptive Step Sizes**: Smaller integration steps near the event horizon for accuracy

### Visual Effects
- **Gravitational Lensing**: Realistic light bending around the black hole
- **Event Horizon**: True black region where light cannot escape (dual horizons for Kerr black holes)
- **Photon Sphere**: Critical orbital radius at 1.5 × Schwarzschild radius
- **Ergosphere Visualization**: Frame-dragging region where spacetime is dragged around rotating black holes
- **Relativistic Jets**: Highly collimated plasma streams with synchrotron emission and relativistic beaming
- **Accretion Disk**: Hot plasma with temperature gradient and Doppler shift effects
- **Background Stars**: Procedural star field showing gravitational lensing distortion

### Music & Sonification
- **Scientific Sonification**: Real-time audio generation from physics calculations
- **Gravitational Wave Audio**: Convert orbital frequencies using Kepler's laws (`f = (1/2π)√(GM/r³)`)
- **Electromagnetic Sonification**: Temperature-based frequencies from accretion disk (`f = k_B T/h`)
- **Relativistic Effects**: Gravitational redshift, Doppler shifts, and frame-dragging audio
- **Musical Scales**: Physics frequencies mapped to pleasant musical scales (pentatonic, chromatic, harmonic series)
- **Multi-Physics Audio**: Simultaneous sonification of jets, tidal forces, and Hawking radiation
- **Spatial Audio**: Stereo positioning based on frame-dragging and observer position
- **Real-time Processing**: 44.1kHz audio with professional-grade EQ, compression, and de-essing

### Technical Implementation
- **Factory Pattern Architecture**: Modular physics engines with runtime switching between Schwarzschild and Kerr
- **Real-time Ray Marching**: GPU-accelerated fragment shader ray tracing
- **Camera System**: Precision orbital camera with emergency escape controls
- **Performance System**: Adaptive quality scaling, turbo mode, and automatic FPS optimization
- **Performance Logging**: Real-time monitoring and analysis of rendering performance
- **Cross-platform**: macOS native with OpenGL 3.3 Core Profile
- **Modular Audio Engine**: Separate sonification system with EQ, compression, and spatial audio

## Physics Background

### Schwarzschild Metric (Non-rotating Black Holes)
**Event Horizon**: rs = 2GM/c² = 2M (in geometric units)
**Photon Sphere**: rph = 1.5rs = 3M
**ISCO**: risco = 3rs = 6M

**Metric Components** (exact implementation):
```cpp
g_tt = -(1 - rs/r)
g_rr = 1/(1 - rs/r) 
g_θθ = r²
g_φφ = r²sin²θ
```

**Geodesic Equations** (actual implementation):
```cpp
d²t/dλ² = -2rs/(r(r-rs)) × dr/dλ × dt/dλ
d²r/dλ² = -rs(r-rs)/(2r³)(dt/dλ)² + rs/(2r(r-rs))(dr/dλ)² + (r-rs)[(dθ/dλ)² + sin²θ(dφ/dλ)²]  
d²θ/dλ² = -2/r × dr/dλ × dθ/dλ + sinθcosθ(dφ/dλ)²
d²φ/dλ² = -2/r × dr/dλ × dφ/dλ - 2cosθ/sinθ × dθ/dλ × dφ/dλ
```

### Kerr Metric (Rotating Black Holes)
**Boyer-Lindquist Coordinates** with spin parameter a:
**Horizons**: r± = M ± √(M² - a²)  
**Ergosphere**: r_ergo = M + √(M² - a²cos²θ)
**Critical Functions**:
- Δ = r² - 2Mr + a²
- Σ = r² + a²cos²θ  
- A = (r² + a²)² - a²Δsin²θ

**Metric Components**:
```cpp
g_tt = -(1 - 2Mr/Σ)
g_rr = Σ/Δ
g_θθ = Σ
g_φφ = (A sin²θ)/Σ
g_tφ = -2Mar sin²θ/Σ  // Frame-dragging term
```

**Frame-Dragging Angular Velocity**:
```cpp
ω = -g_tφ/g_φφ = 2Mar/(r² + a²)² - a²Δsin²θ)
```

### 4th-Order Runge-Kutta Integration
**Implementation** (both Schwarzschild and Kerr):
```cpp
k1 = h × f(x_n, y_n)
k2 = h × f(x_n + h/2, y_n + k1/2)  
k3 = h × f(x_n + h/2, y_n + k2/2)
k4 = h × f(x_n + h, y_n + k3)
y_{n+1} = y_n + (k1 + 2k2 + 2k3 + k4)/6
```

### Sonification Physics (Actual Implementations)

**Orbital Frequency to Audio**:
```cpp
f_orbital = (1/(2π)) × √(GM/r³)    [Kepler's Law]
audio_freq = scale × log(f_orbital/1e-6) / log(1e6)
```

**Electromagnetic to Audio**:
```cpp
f_em = c/λ
log_audio = log_min + (log_em - 6)/(18-6) × (log_max - log_min)
audio_freq = 10^log_audio
```

**Temperature to Audio** (Planck's relation):
```cpp
f_thermal = k_B × T / h
audio_freq = scale × √(f_thermal/1e12)
```

**Synchrotron Emission** (relativistic jets):
```cpp  
f_sync = (eB/2πm_e) × γ²    [Cyclotron × Lorentz²]
audio_freq = scale × log(f_sync/1e9) / log(1e6)
```

### Relativistic Jets Physics

**Beaming Factor** (exact implementation):
```cpp
δ = 1/(γ(1 - β·n̂))    [Relativistic Doppler]
γ = 1/√(1 - β²)       [Lorentz factor]
β = v/c               [Normalized velocity]
```

**Synchrotron Power**:
```cpp
P_sync ∝ B²γ²β²sin²α   [Larmor formula]
Temperature: T = 1e9 × (B/B_0)^(2/3)  [Equipartition]
```

**Beaming Cone Angle**:
```cpp
θ_cone ≈ 1/γ          [Natural beaming cone]
```

### Musical Scale Mapping
**Physics-to-Audio Conversion**:
```cpp
log_normalized = (log₁₀(f_physics) - log₁₀(f_min)) / (log₁₀(f_max) - log₁₀(f_min))
f_musical = A₄ × 2^((scale_note - 69)/12)  [Equal temperament]
```

### GPU Shader Calculations (Real-time Ray Tracing)

**Einstein Deflection Angle** (weak-field approximation):
```glsl
α ≈ 4GM/(c²b) = 2rs/b    [Impact parameter b]
deflectionAngle = 2.0 * schwarzschildRadius / max(impactParameter, 0.1 * rs)
```

**Strong-Field Deflection** (Darwin's exact formula):
```glsl
δφ ≈ 4Rs/b + 15πRs²/(4b²) + 128Rs³/(3b³) + ...  [Higher-order terms]
```

**Impact Parameter Calculation**:
```glsl
impactParameter = distanceToBlackHole × sin(acos(dot(rayDir, toBlackHole)))
```

**Coordinate Transformations** (Boyer-Lindquist ↔ Cartesian):
```cpp
r = √(0.5 × (ρ² + √(ρ²⁺ + 4a²z²)))  where ρ² = x² + y² + z² - a²
x = √(r² + a²) × sinθ × cosφ
y = √(r² + a²) × sinθ × sinφ  
z = r × cosθ
```

### Adaptive Integration Methods

**Schwarzschild Step Size**:
```cpp
if (r < 2rs):     h = 0.001         [Event horizon region]
if (2rs < r < 5rs): h = linear interpolation  [Transition zone] 
if (r > 5rs):     h = 0.1           [Far field]
```

**Kerr Step Size** (includes ergosphere):
```cpp
if (r < 1.1r₊ || inErgosphere): h = 0.001    [Critical regions]
if (r < 2r_ergo):                h = interpolated  [Near ergosphere]
if (r > 2r_ergo):                h = 0.1      [Safe region]
```

**Performance-Adaptive Quality**:
```cpp
if (frameTime > 100ms): performance_scale = 0.2  [Emergency mode]
if (frameTime > 50ms):  performance_scale = 0.4  [Severe degradation] 
zoom_speed = base_speed × (1 - 0.95 × journey_progress²)  [Dramatic slowdown]
```

**Strong-Field Deflection** (exact Darwin series):
```cpp
u = M/b                          [Dimensionless parameter]
α₀ = 4u                         [Weak field (Newton)]
α₁ = (15π/4)u²                  [Post-Newtonian correction]
α₂ = (128/3)u³                  [Strong field correction] 
α₃ = (3465π/64)u⁴               [Higher-order relativistic]
α_total = α₀ + α₁ + α₂ + α₃      [Complete deflection]
```

## System Requirements

- macOS 10.15+ (optimized for Apple Silicon M4 Max)
- OpenGL 3.3 Core Profile support
- CoreAudio framework (built-in on macOS)
- Audio output device (headphones recommended for spatial audio)
- CMake 3.20+
- C++17 compatible compiler (Clang recommended)

## Building

### Quick Build
```bash
./build.sh
```

### Manual Build
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(sysctl -n hw.ncpu)
```

## Running

```bash
cd build
./BlackHoleSimulator
```

## Controls

### Navigation Controls
- **SPACEBAR**: Toggle auto-movement through the black hole (enable/disable autopilot)
- **Mouse Wheel**: Zoom in/out (allows extremely close approach to event horizon)
- **W/S**: Zoom in/out (hold Ctrl for fine control)
- **Q/E**: Precision zoom control (ultra-fine steps)
- **BACKSPACE**: Emergency reverse to safe distance (instant escape)

### Visual Controls
- **1-3**: Change render modes
- **L**: Toggle light ray visualization
- **J**: Toggle relativistic jets on/off
- **U/I**: Increase/decrease jet Lorentz factor
- **O/P**: Increase/decrease jet opening angle
- **V**: Show dramatic visual effects status
- **F**: Test jump to close distance (1.5 Rs)
- **ESC**: Exit

### Performance Controls
- **Q**: Toggle optimization mode (OPTIMIZED ↔ STANDARD)
- **R**: Reset performance settings to defaults
- **T**: TURBO mode (ultra performance for smooth movement)
- **Y**: QUALITY mode (all features enabled)
- **E**: Toggle control panel

### Audio Controls
- **S**: Toggle black hole sonification on/off
- **+/-**: Adjust audio volume
- **M**: Toggle spatial audio effects
- **N**: Toggle relativistic audio effects

## Physics Parameters

### Schwarzschild (Non-rotating) Black Hole
Default configuration with mass equal to our Sun:
- **Mass**: 1 Solar Mass (1.989 × 10³⁰ kg)
- **Spin Parameter**: a = 0.0 (no rotation)
- **Event Horizon**: ~2,954 meters (single horizon)
- **Photon Sphere**: ~4,431 meters (exact calculation using Newton's method)
- **ISCO** (Innermost Stable Circular Orbit): ~8,862 meters

### Kerr (Rotating) Black Hole
Configurable spin parameter up to near-extremal rotation:
- **Spin Parameter Range**: -0.998 ≤ a ≤ 0.998 (dimensionless)
- **Outer Event Horizon**: r₊ = M + √(M² - a²) (depends on spin)
- **Inner Event Horizon**: r₋ = M - √(M² - a²) (Cauchy horizon)
- **Ergosphere**: r_ergo = M + √(M² - a²cos²θ) (frame-dragging region)
- **ISCO**: Spin-dependent, ranges from ~6M (a=0) to ~1.2M (a→1)
- **Frame Dragging**: Spacetime rotation effects with angular velocity Ω

### Audio Parameters
- **Frequency Range**: 20 Hz - 20 kHz (typical hearing range)
- **Sample Rate**: 44.1 kHz professional audio standard
- **Default Musical Scale**: Major Pentatonic (most consonant)
- **Dynamic Range**: 60 dB with professional audio compression
- **Spatial Audio**: Stereo positioning based on frame-dragging effects

## Performance

Optimized for real-time rendering on M4 Max MacBook:
- **Adaptive Ray Marching**: Variable step sizes based on gravitational field strength
- **GPU Fragment Shaders**: Parallel ray computation
- **Efficient Integration**: 4th-order Runge-Kutta with optimized step selection
- **Auto-FPS Optimization**: Dynamic quality scaling to maintain target framerate
- **Turbo Mode**: Ultra-performance mode for smooth extreme proximity navigation
- **Emergency Controls**: Instant escape mechanisms for safe exploration

## Software Architecture

### Factory Pattern Physics Engine
- **Modular Design**: Separate implementations for Schwarzschild and Kerr physics
- **Runtime Switching**: Dynamic selection between physics engines based on spin parameter
- **Interface Abstraction**: Clean separation between physics calculations and rendering
- **Extensible Framework**: Easy addition of new black hole types (e.g., charged Reissner-Nordström)

### Component System
- **BlackHolePhysics**: Core physics calculations and metric computations  
- **RayTracer**: GPU-accelerated ray marching with geodesic integration
- **BlackHoleSonification**: Real-time audio generation from physics parameters
- **RelativisticJets**: Synchrotron emission and relativistic beaming calculations
- **LightSystem**: Photon sphere and gravitational lensing effects
- **PerformanceLogger**: Real-time monitoring and adaptive optimization
- **MelodiousCalculator**: Musical scale mapping for pleasant sonification

### Audio Pipeline
- **Multi-Physics Sources**: Orbital dynamics, thermal emission, jets, and gravitational waves
- **Professional Processing**: EQ, compression, de-essing, and spatial audio
- **Musical Mapping**: Physics frequencies converted to consonant musical scales
- **Real-time Generation**: 44.1 kHz sample rate with low-latency CoreAudio integration

## Scientific Accuracy

This simulation implements exact solutions with validated numerical methods:

### **General Relativity Implementation**
1. **Schwarzschild Metric**: Complete implementation with exact metric components `g_μν`
2. **Kerr Metric**: Full Boyer-Lindquist coordinates including frame-dragging term `g_tφ`
3. **Geodesic Integration**: 4th-order Runge-Kutta with adaptive step sizes near horizons
4. **Event Horizon Detection**: Exact calculation using `r < r_±` for Kerr, `r < 2M` for Schwarzschild

### **Physics Calculations**
1. **Orbital Mechanics**: Kepler's third law `f = (1/2π)√(GM/r³)` with exact GM computation
2. **Photon Spheres**: Newton's method for `r = 1.5rs` (Schwarzschild) and spin-dependent Kerr orbits  
3. **ISCO Computation**: `3rs` for Schwarzschild, complex spin-dependent formula for Kerr
4. **Frame Dragging**: Angular velocity `ω = 2Mar/((r²+a²)²-a²Δsin²θ)` from metric

### **Relativistic Jets**
1. **Synchrotron Emission**: Exact `f = (eB/2πm_e)γ²` with realistic magnetic field profiles
2. **Relativistic Beaming**: Doppler factor `δ = 1/(γ(1-β·n̂))` with proper Lorentz transformations
3. **Magnetic Field Evolution**: `B(r) ∝ (r_h/r)^1.5` with ergosphere amplification for rotating BHs
4. **Temperature Profiles**: Equipartition `T = 10⁹(B/B₀)^(2/3)` K with realistic cooling

### **Audio Sonification Accuracy**
1. **Physics Frequencies**: All derived from fundamental constants (ℏ, k_B, c, G, e, m_e)
2. **Logarithmic Mapping**: Proper frequency scaling from physics ranges to audible spectrum  
3. **Musical Quantization**: Equal temperament `f = A₄ × 2^((n-69)/12)` with consonant scale selection
4. **Relativistic Effects**: Gravitational redshift `√(-g_tt_obs/-g_tt_emit)` and Doppler corrections

### **Numerical Methods**
1. **Adaptive Integration**: Step size `h ∝ 1/√(|acceleration|)` near strong-field regions
2. **Coordinate Transformations**: Exact spherical ↔ Cartesian with singularity handling
3. **Metric Evaluation**: Cache-optimized computation of `Σ`, `Δ`, `A` functions for Kerr
4. **Professional Audio**: 44.1 kHz sampling with anti-aliasing and musical EQ processing

## References

### Physics & General Relativity
- Schwarzschild, K. (1916). "Über das Gravitationsfeld eines Massenpunktes nach der Einsteinschen Theorie"
- Kerr, R. P. (1963). "Gravitational field of a spinning mass as an example of algebraically special metrics"
- Chandrasekhar, S. (1983). "The Mathematical Theory of Black Holes"
- Luminet, J.-P. (1979). "Image of a spherical black hole with thin accretion disk"

### Sonification & Musical Theory
- Hermann, T. (2008). "Taxonomy and definitions for sonification and auditory display"
- Kramer, G. (1994). "Auditory Display: Sonification, Audification, and Auditory Interfaces"
- Dubus, G. & Bresin, R. (2013). "A systematic review of mapping strategies for the sonification of physical quantities"
- Milne, A. (2013). "A computational model of the cognition of tonality"

## Future Enhancements

### Visual Features
- [ ] Gravitational wave visualization ripples in spacetime
- [ ] Multiple black hole systems and binary mergers
- [ ] Charged black holes (Reissner-Nordström metric)
- [ ] Wormhole traversability visualization
- [ ] VR/AR support for immersive exploration

### Physics Extensions
- [ ] Tidal force visualization on extended objects
- [ ] Hawking radiation particle animation
- [ ] Frame-dragging field line visualization
- [ ] Plasma magnetohydrodynamics in jets

### Audio Features
- [ ] MIDI output for external synthesizers
- [ ] Real-time gravitational wave data sonification (LIGO/Virgo data)
- [ ] Interactive physics parameter control via audio feedback
- [ ] Binaural beats from orbital resonances and beat frequencies
- [ ] Machine learning-optimized pleasant frequency mapping
- [ ] Sonification of gravitational memory effects

## License

MIT License - See LICENSE file for details.

---

*"The most beautiful thing we can experience is the mysterious." - Albert Einstein*
