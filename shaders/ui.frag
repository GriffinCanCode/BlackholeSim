#version 330 core

in vec2 TexCoord;
out vec4 FragColor;

uniform vec4 color;
uniform float alpha;
uniform int renderMode;  // 0 = solid color, 1 = text rendering

void main() {
    if (renderMode == 0) {
        // Solid color panel background
        FragColor = vec4(color.rgb, color.a * alpha);
    } else {
        // Simple text rendering - create a basic bitmap-style font
        vec2 uv = TexCoord;
        vec2 grid = floor(uv * vec2(8.0, 16.0));
        
        // Simple procedural text - just create some basic shapes for now
        float text = 0.0;
        
        // Create a simple border and background for text areas
        if (uv.x < 0.05 || uv.x > 0.95 || uv.y < 0.05 || uv.y > 0.95) {
            text = 1.0; // Border
        } else {
            text = 0.3; // Background
        }
        
        FragColor = vec4(color.rgb * text, alpha);
    }
}
