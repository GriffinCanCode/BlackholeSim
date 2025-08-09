cmake_minimum_required(VERSION 3.20)
project(BlackHoleSimulator VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find packages
find_package(OpenGL REQUIRED)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/main.cpp
    src/RayTracer.cpp
    src/Renderer.cpp
    src/Physics.cpp
    src/Shader.cpp
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link libraries - macOS specific
if(APPLE)
    target_link_libraries(${PROJECT_NAME} 
        OpenGL::GL
        "-framework Cocoa"
        "-framework IOKit"
        "-framework CoreVideo"
    )
    target_compile_definitions(${PROJECT_NAME} PRIVATE GL_SILENCE_DEPRECATION)
    
    # Download and build GLFW from source for M4 Max compatibility
    include(FetchContent)
    FetchContent_Declare(
        glfw
        GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG 3.3.8
    )
    FetchContent_MakeAvailable(glfw)
    target_link_libraries(${PROJECT_NAME} glfw)
else()
    # Linux/Windows
    find_package(glfw3 REQUIRED)
    target_link_libraries(${PROJECT_NAME} OpenGL::GL glfw)
endif()

# Compiler flags for optimization
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(${PROJECT_NAME} PRIVATE -O3 -march=native -ffast-math)
endif()
