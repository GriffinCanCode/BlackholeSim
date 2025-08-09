# Black Hole Simulator

A physically accurate 3D black hole visualization and sonification system using C++ and OpenGL, implementing real general relativity physics with scientifically accurate audio generation.

## Features

### Physics Accuracy
- **Schwarzschild Metric**: Accurate spacetime curvature around non-rotating black holes
- **Geodesic Ray Tracing**: Light paths calculated using Einstein's field equations
- **4th-Order Runge-Kutta Integration**: Numerical solution of geodesic differential equations
- **Adaptive Step Sizes**: Smaller integration steps near the event horizon for accuracy

### Visual Effects
- **Gravitational Lensing**: Realistic light bending around the black hole
- **Event Horizon**: True black region where light cannot escape
- **Photon Sphere**: Critical orbital radius at 1.5 × Schwarzschild radius
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
- **Real-time Ray Marching**: GPU-accelerated fragment shader ray tracing
- **Camera System**: Free-look camera with orbital controls
- **Performance Optimized**: Adaptive algorithms designed for M4 Max MacBook
- **Cross-platform**: macOS native with OpenGL 3.3 Core Profile

## Physics Background

### Schwarzschild Radius
The event horizon radius is calculated as:
```
rs = 2GM/c²
```
Where G is the gravitational constant, M is the black hole mass, and c is the speed of light.

### Geodesic Equations
Light rays follow null geodesics in curved spacetime, described by:
```
d²xᵘ/dλ² + Γᵘₘᵥ (dxᵐ/dλ)(dxᵥ/dλ) = 0
```
Where Γᵘₘᵥ are the Christoffel symbols derived from the Schwarzschild metric.

### Sonification Physics
The audio system converts real physics calculations into sound using exact formulas:

**Orbital Dynamics**:
```
f_orbital = (1/2π) × √(GM/r³)     [Kepler's Third Law]
f_gravitational_wave = 2 × f_orbital   [Quadrupole Emission]
```

**Thermal Emission**:
```
f_thermal = k_B × T / h           [Planck's Energy Relation]
T_accretion(r) ∝ r^(-3/4)        [Shakura-Sunyaev Profile]
```

**Relativistic Effects**:
```
f_redshift = f_emit × √(-g_tt(r_obs)/-g_tt(r_emit))  [Gravitational Redshift]
f_doppler = f_emit × γ(1 + β·n̂)                      [Relativistic Doppler]
T_hawking = ℏc³/(8πGMk_B)                            [Hawking Temperature]
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

### Camera Controls
- **WASD**: Move camera forward/back/left/right
- **QE**: Move camera up/down  
- **Mouse**: Look around (free camera)
- **Scroll**: Adjust field of view

### Visual Controls
- **1-3**: Change render modes
- **ESC**: Exit

### Audio Controls
- **M**: Toggle music/sonification on/off
- **+/-**: Adjust audio volume
- **P**: Change musical scale (Pentatonic → Chromatic → Harmonic Series)
- **R**: Toggle relativistic audio effects
- **S**: Toggle spatial audio (stereo frame-dragging effects)

## Physics Parameters

The simulation uses a black hole with mass equal to our Sun by default:
- **Mass**: 1 Solar Mass (1.989 × 10³⁰ kg)
- **Schwarzschild Radius**: ~2,954 meters
- **Photon Sphere**: ~4,431 meters (exact calculation using Newton's method)
- **ISCO** (Innermost Stable Circular Orbit): ~8,862 meters

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

## Scientific Accuracy

This simulation implements:
1. **General Relativity**: Exact Schwarzschild and Kerr solutions to Einstein's field equations
2. **Null Geodesics**: Exact photon trajectories using 4th-order Runge-Kutta integration  
3. **Strong-Field Physics**: Darwin's exact deflection angles and Newton's method for photon spheres
4. **Sonification Accuracy**: All audio frequencies derived from real physics calculations
5. **Relativistic Audio**: Gravitational redshift, Doppler shifts, time dilation, and frame-dragging
6. **Musical Physics**: Frequencies mapped to consonant scales using harmonic theory

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
- [ ] Gravitational wave visualization  
- [ ] Multiple black hole systems
- [ ] Relativistic plasma dynamics
- [ ] VR/AR support

### Audio Features
- [ ] MIDI output for external synthesizers
- [ ] Real-time gravitational wave data sonification
- [ ] Interactive physics parameter control via audio
- [ ] Binaural beats from orbital resonances
- [ ] Machine learning-optimized pleasant frequency mapping

## License

MIT License - See LICENSE file for details.

---

*"The most beautiful thing we can experience is the mysterious." - Albert Einstein*
