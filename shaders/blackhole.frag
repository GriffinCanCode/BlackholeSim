#version 330 core
out vec4 FragColor;

in vec2 TexCoord;

uniform vec3 cameraPos;
uniform vec3 blackholePos;
uniform float schwarzschildRadius;
uniform vec2 resolution;
uniform int debugMode;


const float PI = 3.14159265359;

// Generate procedural stars
vec3 generateStars(vec3 direction) {
    vec3 d = normalize(direction);
    float hash = fract(sin(dot(d.xy, vec2(12.9898, 78.233)) + d.z * 37.719) * 43758.5453);
    
    if (hash > 0.998) {
        return vec3(1.0, 1.0, 0.9); // Bright star
    } else if (hash > 0.995) {
        return vec3(0.8, 0.8, 0.6); // Medium star
    } else if (hash > 0.99) {
        return vec3(0.6, 0.5, 0.4); // Dim star
    } else {
        // Space background
        return vec3(0.02, 0.02, 0.05);
    }
}

// Simple accretion disk
vec3 sampleAccretionDisk(vec3 rayPos, vec3 rayDir) {
    vec3 diskCenter = blackholePos;
    vec3 diskNormal = vec3(0, 0, 1); // Disk in XY plane
    
    // Ray-plane intersection
    float denom = dot(rayDir, diskNormal);
    if (abs(denom) < 0.001) return vec3(0.0); // Ray parallel to disk
    
    float t = dot(diskCenter - rayPos, diskNormal) / denom;
    if (t < 0.0) return vec3(0.0); // Behind ray
    
    vec3 hitPoint = rayPos + t * rayDir;
    float distFromCenter = length(hitPoint.xy - diskCenter.xy);
    
    float innerRadius = 3.0 * schwarzschildRadius; // ISCO
    float outerRadius = 12.0 * schwarzschildRadius;
    float diskThickness = 0.2 * schwarzschildRadius;
    
    if (distFromCenter >= innerRadius && distFromCenter <= outerRadius && 
        abs(hitPoint.z - diskCenter.z) <= diskThickness) {
        
        // Temperature gradient (hotter = bluer, cooler = redder)
        float temp = outerRadius / distFromCenter;
        vec3 diskColor = vec3(
            0.3 + 0.7 * (1.0 - temp),  // Red decreases with temperature
            0.4 + 0.6 * temp,          // Green increases with temperature  
            0.6 + 0.4 * temp           // Blue increases with temperature
        );
        
        // Add rotation and turbulence
        float angle = atan(hitPoint.y - diskCenter.y, hitPoint.x - diskCenter.x);
        float rotation = angle + 0.5 * sqrt(schwarzschildRadius / distFromCenter); // Static rotation
        float turbulence = 0.3 * sin(rotation * 8.0) * sin(distFromCenter * 0.5);
        
        return diskColor * (0.8 + turbulence);
    }
    
    return vec3(0.0);
}

// Simplified gravitational lensing
vec3 applyGravitationalLensing(vec3 rayPos, vec3 rayDir) {
    vec3 toBlackHole = blackholePos - rayPos;
    float distanceToBlackHole = length(toBlackHole);
    
    // Don't apply lensing if too close to avoid numerical issues
    if (distanceToBlackHole < 2.0 * schwarzschildRadius) {
        return rayDir;
    }
    
    // Calculate impact parameter (perpendicular distance)
    vec3 normalizedDir = normalize(rayDir);
    vec3 normalizedToBlackHole = normalize(toBlackHole);
    
    // Simple deflection based on distance
    float deflectionStrength = schwarzschildRadius / (distanceToBlackHole * distanceToBlackHole);
    deflectionStrength = min(deflectionStrength, 0.1); // Clamp for stability
    
    // Apply deflection perpendicular to ray direction
    vec3 perpendicular = normalize(cross(normalizedDir, cross(normalizedDir, normalizedToBlackHole)));
    vec3 deflectedDir = normalize(normalizedDir + perpendicular * deflectionStrength);
    
    return deflectedDir;
}

void main() {
    // Debug modes for testing
    if (debugMode == 1) {
        // Test pattern mode - show screen coordinates
        vec2 screenPos;
        screenPos.x = (gl_FragCoord.x / resolution.x) * 2.0 - 1.0;
        screenPos.y = (gl_FragCoord.y / resolution.y) * 2.0 - 1.0;
        vec3 color = vec3(0.1, 0.1, 0.2);
        
        // Show center as bright white
        if (abs(screenPos.x) < 0.05 && abs(screenPos.y) < 0.05) {
            color = vec3(1.0, 1.0, 1.0);
        }
        // Show coordinate grid
        color.r += abs(screenPos.x) * 0.5;
        color.g += abs(screenPos.y) * 0.5;
        FragColor = vec4(color, 1.0);
        return;
    } else if (debugMode == 2) {
        // Solid color mode
        FragColor = vec4(0.0, 1.0, 0.0, 1.0); // Green screen
        return;
    }
    
    // Normal rendering mode - Simple coordinate system
    vec2 screenCoords;
    screenCoords.x = (gl_FragCoord.x / resolution.x) * 2.0 - 1.0;
    screenCoords.y = (gl_FragCoord.y / resolution.y) * 2.0 - 1.0;
    
    // Simple ray tracing - camera at cameraPos looking at origin
    float aspectRatio = resolution.x / resolution.y;
    float fov = 45.0 * 3.14159265359 / 180.0;
    float tanHalfFov = tan(fov * 0.5);
    
    // Apply aspect ratio to screen coordinates
    screenCoords.x *= aspectRatio;
    
    // Create ray direction pointing from camera toward screen point at origin
    vec3 rayDir = normalize(vec3(
        screenCoords.x * tanHalfFov,
        screenCoords.y * tanHalfFov,
        -1.0
    ));
    
    // Distance from camera to black hole
    float distanceToBlackHole = length(cameraPos - blackholePos);
    
    // Simple ray-sphere intersection for event horizon
    vec3 oc = cameraPos - blackholePos;
    float a = dot(rayDir, rayDir);
    float b = 2.0 * dot(oc, rayDir);
    float c = dot(oc, oc) - schwarzschildRadius * schwarzschildRadius;
    float discriminant = b * b - 4.0 * a * c;
    
    // Start with star field background
    vec3 color = generateStars(rayDir);
    
    // Check if ray hits black hole
    if (discriminant >= 0.0) {
        float t = (-b - sqrt(discriminant)) / (2.0 * a);
        if (t > 0.0) {
            // Ray hits event horizon - pure black
            color = vec3(0.0, 0.0, 0.0);
        }
    } else {
        // Add orange glow around black hole for visibility
        vec3 toBlackHole = normalize(blackholePos - cameraPos);
        float alignment = dot(rayDir, toBlackHole);
        if (alignment > 0.95) { // Cone pointing toward black hole
            float glow = (alignment - 0.95) / 0.05; // Normalize to 0-1
            color += vec3(glow * 1.5, glow * 0.8, glow * 0.2); // Orange glow
        }
    }
    
    FragColor = vec4(color, 1.0);
}
