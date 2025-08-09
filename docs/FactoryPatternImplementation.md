# Factory Pattern Implementation for BlackHole Physics

### 1. Concrete Implementations

#### SchwarzchildPhysics
- Pure Schwarzschild (non-rotating) black hole implementation
- Optimized for non-rotating cases
- Direct mathematical calculations without Kerr complexity

#### KerrPhysicsEngine
- Rotating black hole implementation
- Delegates to existing KerrPhysics for calculations
- Handles complex frame-dragging and ergosphere effects

### 2. Factory Classes

#### BlackHolePhysicsFactory
```cpp
class BlackHolePhysicsFactory {
public:
    enum class PhysicsType { Schwarzschild, Kerr, Auto };
    
    static std::unique_ptr<IBlackHolePhysics> create(
        double mass, double spin = 0.0, PhysicsType type = Auto);
    static std::unique_ptr<IBlackHolePhysics> createSchwarzschild(double mass);
    static std::unique_ptr<IBlackHolePhysics> createKerr(double mass, double spin);
};
```

**Features:**
- Auto-selection based on spin parameter
- Explicit creation methods for specific types
- Type name utilities and metadata

#### BlackHolePhysicsRegistry (Extensibility)
```cpp
class BlackHolePhysicsRegistry {
public:
    using PhysicsCreator = std::function<std::unique_ptr<IBlackHolePhysics>(double, double)>;
    
    static void registerPhysicsType(const std::string& name, PhysicsCreator creator);
    static std::unique_ptr<IBlackHolePhysics> createByName(const std::string& name, double mass, double spin);
    static std::vector<std::string> getAvailableTypes();
};
```

**Benefits:**
- Extensible for future physics models (Reissner-Nordström, etc.)
- Plugin-like architecture for custom implementations
- Runtime discovery of available physics types

## Technical Details

### Memory Management
- Smart pointers (`std::unique_ptr`) for automatic cleanup
- RAII principles throughout
- No memory leaks or dangling pointers

### Performance Considerations
- Factory creation is O(1) with minimal overhead
- Method delegation is inlined by compiler (zero cost)
- Adaptive step sizes preserved for numerical integration
- Memory footprint comparable to original implementation

### Thread Safety
- Factory methods are thread-safe (stateless static methods)
- Physics engine instances are not thread-safe (same as before)
- Registry uses static data (initialization order guaranteed)

## Build Integration

### CMakeLists.txt Changes
```cmake
set(SOURCES
    src/main.cpp
    src/Physics.cpp
    src/PhysicsFactory.cpp  # Added new source file
    src/KerrPhysics.cpp
    # ... other sources
)
```

### File Structure
```
include/
├── Physics.h              # Updated facade class
├── PhysicsFactory.h       # New factory interfaces
└── KerrPhysics.h          # Unchanged

src/
├── Physics.cpp            # Updated to use factory internally
├── PhysicsFactory.cpp     # New factory implementations  
└── KerrPhysics.cpp        # Unchanged
```