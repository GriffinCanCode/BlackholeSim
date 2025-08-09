#!/bin/bash

echo "=== Black Hole Simulator Build Script ==="
echo "Building for macOS M4 Max with native OpenGL"

# Create build directory
mkdir -p build
cd build

# Configure with CMake
echo "Configuring project..."
cmake -DCMAKE_BUILD_TYPE=Release ..

if [ $? -ne 0 ]; then
    echo "CMake configuration failed!"
    exit 1
fi

# Build the project
echo "Building project..."
make -j$(sysctl -n hw.ncpu)

if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi

echo "Build completed successfully!"
echo "Executable: ./BlackHoleSimulator"
echo ""
echo "To run the simulation:"
echo "cd build && ./BlackHoleSimulator"
