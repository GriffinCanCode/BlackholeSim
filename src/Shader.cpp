#include "Shader.h"
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#else
#include <GL/gl.h>
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

Shader::Shader() : program_id_(0) {}

Shader::~Shader() {
    if (program_id_ != 0) {
        glDeleteProgram(program_id_);
    }
}

bool Shader::loadFromFiles(const std::string& vertex_path, const std::string& fragment_path) {
    std::string vertex_source = loadFileToString(vertex_path);
    std::string fragment_source = loadFileToString(fragment_path);
    
    if (vertex_source.empty() || fragment_source.empty()) {
        std::cerr << "ERROR: Failed to load shader files" << std::endl;
        return false;
    }
    
    return loadFromStrings(vertex_source, fragment_source);
}

bool Shader::loadFromStrings(const std::string& vertex_source, const std::string& fragment_source) {
    // Compile vertex shader
    unsigned int vertex_shader = compileShader(vertex_source, GL_VERTEX_SHADER);
    if (vertex_shader == 0) return false;
    
    // Compile fragment shader
    unsigned int fragment_shader = compileShader(fragment_source, GL_FRAGMENT_SHADER);
    if (fragment_shader == 0) {
        glDeleteShader(vertex_shader);
        return false;
    }
    
    // Create program and link shaders
    program_id_ = glCreateProgram();
    glAttachShader(program_id_, vertex_shader);
    glAttachShader(program_id_, fragment_shader);
    glLinkProgram(program_id_);
    
    // Check for linking errors
    int success;
    char info_log[512];
    glGetProgramiv(program_id_, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(program_id_, 512, nullptr, info_log);
        std::cerr << "ERROR: Shader program linking failed\n" << info_log << std::endl;
        glDeleteProgram(program_id_);
        program_id_ = 0;
    }
    
    // Clean up individual shaders
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    
    return success;
}

void Shader::use() const {
    if (program_id_ != 0) {
        glUseProgram(program_id_);
    }
}

void Shader::setInt(const std::string& name, int value) {
    glUniform1i(getUniformLocation(name), value);
}

void Shader::setFloat(const std::string& name, float value) {
    glUniform1f(getUniformLocation(name), value);
}

void Shader::setVec2(const std::string& name, float x, float y) {
    glUniform2f(getUniformLocation(name), x, y);
}

void Shader::setVec3(const std::string& name, float x, float y, float z) {
    glUniform3f(getUniformLocation(name), x, y, z);
}

void Shader::setVec4(const std::string& name, float x, float y, float z, float w) {
    glUniform4f(getUniformLocation(name), x, y, z, w);
}

void Shader::setMat4(const std::string& name, const float* matrix) {
    glUniformMatrix4fv(getUniformLocation(name), 1, GL_FALSE, matrix);
}

int Shader::getUniformLocation(const std::string& name) {
    auto it = uniform_cache_.find(name);
    if (it != uniform_cache_.end()) {
        return it->second;
    }
    
    int location = glGetUniformLocation(program_id_, name.c_str());
    uniform_cache_[name] = location;
    
    if (location == -1) {
        std::cerr << "Warning: Uniform '" << name << "' not found in shader" << std::endl;
    }
    
    return location;
}

unsigned int Shader::compileShader(const std::string& source, unsigned int type) {
    unsigned int shader = glCreateShader(type);
    const char* src = source.c_str();
    glShaderSource(shader, 1, &src, nullptr);
    glCompileShader(shader);
    
    // Check for compilation errors
    int success;
    char info_log[512];
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(shader, 512, nullptr, info_log);
        std::string shader_type = (type == GL_VERTEX_SHADER) ? "VERTEX" : "FRAGMENT";
        std::cerr << "ERROR: " << shader_type << " shader compilation failed\n" << info_log << std::endl;
        glDeleteShader(shader);
        return 0;
    }
    
    return shader;
}

std::string Shader::loadFileToString(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open file: " << filepath << std::endl;
        return "";
    }
    
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}
