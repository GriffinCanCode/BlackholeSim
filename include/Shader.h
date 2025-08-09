#pragma once
#include <string>
#include <unordered_map>

class Shader {
public:
    Shader();
    ~Shader();
    
    // Load and compile shaders
    bool loadFromFiles(const std::string& vertex_path, const std::string& fragment_path);
    bool loadFromStrings(const std::string& vertex_source, const std::string& fragment_source);
    
    // Use the shader program
    void use() const;
    
    // Uniform setters
    void setInt(const std::string& name, int value);
    void setFloat(const std::string& name, float value);
    void setVec2(const std::string& name, float x, float y);
    void setVec3(const std::string& name, float x, float y, float z);
    void setVec4(const std::string& name, float x, float y, float z, float w);
    void setMat4(const std::string& name, const float* matrix);
    
    // Get uniform location
    int getUniformLocation(const std::string& name);
    
    // Get program ID
    unsigned int getProgramID() const { return program_id_; }
    
    // Check if shader is valid
    bool isValid() const { return program_id_ != 0; }
    
private:
    unsigned int program_id_;
    std::unordered_map<std::string, int> uniform_cache_;
    
    // Helper functions
    unsigned int compileShader(const std::string& source, unsigned int type);
    std::string loadFileToString(const std::string& filepath);
    void checkCompileErrors(unsigned int shader, const std::string& type);
};
