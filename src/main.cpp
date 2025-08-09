#include "../include/Renderer.h"
#include <iostream>
#include <memory>

int main() {
    std::cout << "=== Black Hole Simulator ===" << std::endl;
    std::cout << "A physically accurate 3D black hole visualization" << std::endl;
    std::cout << "Using Schwarzschild metric and geodesic ray tracing" << std::endl;
    std::cout << "=================================" << std::endl;
    
    // Create renderer
    auto renderer = std::make_unique<BlackHoleRenderer>(1200, 800);
    
    // Initialize the renderer
    if (!renderer->initialize()) {
        std::cerr << "Failed to initialize renderer!" << std::endl;
        return -1;
    }
    
    std::cout << "Initialization complete. Starting simulation..." << std::endl;
    
    // Main render loop
    while (!renderer->shouldClose()) {
        // Process input
        renderer->pollEvents();
        
        // Render frame
        renderer->render();
        
        // Swap buffers
        renderer->swapBuffers();
    }
    
    std::cout << "Shutting down..." << std::endl;
    return 0;
}
