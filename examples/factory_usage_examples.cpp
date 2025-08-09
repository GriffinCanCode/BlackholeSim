// Examples demonstrating how to use the new Factory Pattern in BlackholeSim

#include "PhysicsFactory.h"
#include "Physics.h"  // For backward compatibility
#include <iostream>
#include <memory>

void demonstrateBasicFactoryUsage() {
    std::cout << "=== Basic Factory Usage ===" << std::endl;
    
    // Method 1: Direct factory usage (new way)
    auto schwarzschild = BlackHolePhysicsFactory::createSchwarzschild(1.0);
    auto kerr = BlackHolePhysicsFactory::createKerr(1.0, 0.7);
    
    std::cout << "Schwarzschild physics type: " << schwarzschild->getPhysicsType() << std::endl;
    std::cout << "Kerr physics type: " << kerr->getPhysicsType() << std::endl;
    
    // Method 2: Auto-selection factory (recommended)
    auto auto_schwarzschild = BlackHolePhysicsFactory::create(1.0, 0.0);  // Small spin -> Schwarzschild
    auto auto_kerr = BlackHolePhysicsFactory::create(1.0, 0.5);           // Large spin -> Kerr
    
    std::cout << "Auto-selected for spin=0.0: " << auto_schwarzschild->getPhysicsType() << std::endl;
    std::cout << "Auto-selected for spin=0.5: " << auto_kerr->getPhysicsType() << std::endl;
}

void demonstrateBackwardCompatibility() {
    std::cout << "\n=== Backward Compatibility ===" << std::endl;
    
    // This still works exactly as before - no changes needed to existing code!
    BlackHolePhysics physics1(1.0, 0.0);  // Non-rotating black hole
    BlackHolePhysics physics2(1.0, 0.8);  // Rotating black hole
    
    std::cout << "Legacy physics1 type: " << physics1.getPhysicsType() << std::endl;
    std::cout << "Legacy physics1 horizon: " << physics1.schwarzschildRadius() << std::endl;
    
    std::cout << "Legacy physics2 type: " << physics2.getPhysicsType() << std::endl;
    std::cout << "Legacy physics2 outer horizon: " << physics2.outerHorizonRadius() << std::endl;
    std::cout << "Legacy physics2 inner horizon: " << physics2.innerHorizonRadius() << std::endl;
    
    // The facade class provides access to the underlying factory-created engine
    const IBlackHolePhysics* engine = physics2.getPhysicsEngine();
    std::cout << "Underlying engine spin: " << engine->getSpinParameter() << std::endl;
}

void demonstratePhysicsComparison() {
    std::cout << "\n=== Physics Comparison ===" << std::endl;
    
    double mass = 1.0;  // 1 solar mass
    double spin = 0.9;  // High spin
    
    auto schwarzschild = BlackHolePhysicsFactory::createSchwarzschild(mass);
    auto kerr = BlackHolePhysicsFactory::createKerr(mass, spin);
    
    std::cout << "Black hole mass: " << mass << " solar masses" << std::endl;
    std::cout << "Kerr spin parameter: " << spin << std::endl << std::endl;
    
    // Compare basic properties
    std::cout << "Property Comparison:" << std::endl;
    std::cout << "Schwarzschild radius: " << schwarzschild->schwarzschildRadius() << std::endl;
    std::cout << "Kerr outer horizon: " << kerr->outerHorizonRadius() << std::endl;
    std::cout << "Kerr inner horizon: " << kerr->innerHorizonRadius() << std::endl;
    
    std::cout << "\nPhoton Sphere:" << std::endl;
    std::cout << "Schwarzschild: " << schwarzschild->photonSphere() << std::endl;
    std::cout << "Kerr: " << kerr->photonSphere() << std::endl;
    
    std::cout << "\nISCO (Innermost Stable Circular Orbit):" << std::endl;
    std::cout << "Schwarzschild: " << schwarzschild->innerStableCircularOrbit() << std::endl;
    std::cout << "Kerr: " << kerr->innerStableCircularOrbit() << std::endl;
    
    // Kerr-specific properties
    std::cout << "\nKerr-specific properties:" << std::endl;
    std::cout << "Ergosphere radius (equator): " << kerr->ergosphereRadius(M_PI/2) << std::endl;
    std::cout << "Is point at r=1.5 in ergosphere? " << 
                 (kerr->isInErgosphere(1.5, M_PI/2) ? "Yes" : "No") << std::endl;
}

void demonstrateRegistryExtension() {
    std::cout << "\n=== Registry Extension (Future Physics Models) ===" << std::endl;
    
    // Show how to extend with custom physics models
    std::cout << "Available physics types: " << std::endl;
    auto types = BlackHolePhysicsRegistry::getAvailableTypes();
    for (const auto& type : types) {
        std::cout << "  - " << type << std::endl;
    }
    
    // Example of how to register a custom physics model (placeholder)
    std::cout << "\nNote: Custom physics models can be registered using:" << std::endl;
    std::cout << "BlackHolePhysicsRegistry::registerPhysicsType(\"Custom\", customCreator);" << std::endl;
}

int main() {
    std::cout << "BlackHole Physics Factory Pattern Demo" << std::endl;
    std::cout << "=====================================" << std::endl;
    
    try {
        demonstrateBasicFactoryUsage();
        demonstrateBackwardCompatibility();
        demonstratePhysicsComparison();
        demonstrateRegistryExtension();
        
        std::cout << "\n=== All Examples Completed Successfully! ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

/* 
Expected Output:
================

=== Basic Factory Usage ===
Schwarzschild physics type: Schwarzschild
Kerr physics type: Kerr
Auto-selected for spin=0.0: Schwarzschild
Auto-selected for spin=0.5: Kerr

=== Backward Compatibility ===
Legacy physics1 type: Schwarzschild
Legacy physics1 horizon: 2
Legacy physics2 type: Kerr
Legacy physics2 outer horizon: [Kerr outer horizon radius]
Legacy physics2 inner horizon: [Kerr inner horizon radius]
Underlying engine spin: 0.8

=== Physics Comparison ===
Black hole mass: 1 solar masses
Kerr spin parameter: 0.9

Property Comparison:
Schwarzschild radius: 2
Kerr outer horizon: [computed value]
Kerr inner horizon: [computed value]

Photon Sphere:
Schwarzschild: 3
Kerr: [Kerr photon sphere]

ISCO (Innermost Stable Circular Orbit):
Schwarzschild: 6
Kerr: [Kerr ISCO]

Kerr-specific properties:
Ergosphere radius (equator): [computed value]
Is point at r=1.5 in ergosphere? Yes

=== Registry Extension (Future Physics Models) ===
Available physics types: 
  - Schwarzschild
  - Kerr

Note: Custom physics models can be registered using:
BlackHolePhysicsRegistry::registerPhysicsType("Custom", customCreator);

=== All Examples Completed Successfully! ===
*/
