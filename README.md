# Black Hole Simulator

A physically accurate 3D black hole visualization using C++ and OpenGL, implementing real general relativity physics.

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

## System Requirements

- macOS 10.15+ (optimized for Apple Silicon M4 Max)
- OpenGL 3.3 Core Profile support
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

- **WASD**: Move camera forward/back/left/right
- **QE**: Move camera up/down  
- **Mouse**: Look around (free camera)
- **Scroll**: Adjust field of view
- **1-3**: Change render modes
- **ESC**: Exit

## Physics Parameters

The simulation uses a black hole with mass equal to our Sun by default:
- **Mass**: 1 Solar Mass (1.989 × 10³⁰ kg)
- **Schwarzschild Radius**: ~2,954 meters
- **Photon Sphere**: ~4,431 meters
- **ISCO** (Innermost Stable Circular Orbit): ~8,862 meters

## Performance

Optimized for real-time rendering on M4 Max MacBook:
- **Adaptive Ray Marching**: Variable step sizes based on gravitational field strength
- **GPU Fragment Shaders**: Parallel ray computation
- **Efficient Integration**: 4th-order Runge-Kutta with optimized step selection

## Scientific Accuracy

This simulation implements:
1. **General Relativity**: Schwarzschild solution to Einstein's field equations
2. **Null Geodesics**: Exact photon trajectories in curved spacetime  
3. **Gravitational Redshift**: Frequency shifts due to gravitational potential
4. **Frame Dragging Effects**: (Planned for Kerr metric implementation)

## References

- Schwarzschild, K. (1916). "Über das Gravitationsfeld eines Massenpunktes nach der Einsteinschen Theorie"
- Chandrasekhar, S. (1983). "The Mathematical Theory of Black Holes"
- Luminet, J.-P. (1979). "Image of a spherical black hole with thin accretion disk"

## Future Enhancements

- [ ] Kerr metric for rotating black holes
- [ ] Gravitational wave visualization  
- [ ] Multiple black hole systems
- [ ] Relativistic plasma dynamics
- [ ] VR/AR support

## License

MIT License - See LICENSE file for details.

---

*"The most beautiful thing we can experience is the mysterious." - Albert Einstein*
