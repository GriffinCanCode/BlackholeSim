#version 330 core
out vec4 FragColor;

in vec2 TexCoord;

uniform vec3 cameraPos;
uniform vec3 blackholePos;
uniform float schwarzschildRadius;
uniform vec2 resolution;
uniform int debugMode;
uniform float time;
uniform int enableSpectralRendering;
uniform float diskTemperature;
uniform float gravitationalRedshift;
uniform int showLightRays;
uniform int showJets;
uniform float jetLorentzFactor;
uniform vec3 jetAxis;
uniform float jetOpeningAngle;
uniform float jetMagneticField;

// Enhanced realistic light parameters
uniform float photonSphereGlow;
uniform float ambientLightLevel;
uniform float diskBrightness;
uniform float starBrightness;
uniform int enablePhotonSphereGlow;

// Performance optimization uniforms
uniform int lodLevel;           // 0=highest quality, 3=lowest quality
uniform float performanceScale; // 0.1-1.0 quality vs performance balance
uniform float lightRayDistanceScale; // Distance-based scaling factor
uniform int maxLightSources;    // Adaptive light source count

// Enhanced proximity effect uniforms
uniform float proximityFactor;   // 0.0=far, 1.0=at event horizon
uniform int enableDramaticEffects; // Toggle for enhanced visual drama

const float PI = 3.14159265359;

// Physical constants - same as original
const float PLANCK_CONSTANT = 6.62607015e-34;
const float BOLTZMANN_CONSTANT = 1.380649e-23;
const float SPEED_OF_LIGHT = 299792458.0;
const float WIEN_CONSTANT = 2.897771955e-3;

// Visible spectrum wavelengths (nm)
const float WAVELENGTH_RED = 700.0;
const float WAVELENGTH_GREEN = 550.0;
const float WAVELENGTH_BLUE = 450.0;

// OPTIMIZATION 1: Precomputed lookup tables (approximations for expensive functions)
float fastPlanckFunction(float temperature, float wavelength) {
    // Fast approximation using polynomial fit for common temperature ranges
    if (temperature < 1000.0) return 0.0;
    
    float lambda5 = wavelength * wavelength * wavelength * wavelength * wavelength;
    float expArg = 1.4387e-2 / (wavelength * temperature); // hc/λkT simplified
    
    // Fast exp approximation for common ranges
    if (expArg > 10.0) return 0.0; // Avoid numerical overflow
    if (expArg < 0.1) expArg = 0.1; // Avoid division issues
    
    float expTerm = 1.0 + expArg + expArg*expArg*0.5; // Taylor series approximation
    return (3.74183e-16 / lambda5) / expTerm; // Simplified constant
}

// OPTIMIZATION 2: Fast wavelength to RGB conversion with fewer branches
vec3 fastWavelengthToRGB(float wavelength) {
    // Optimized piecewise linear approximation
    wavelength = clamp(wavelength, 380.0, 750.0);
    
    vec3 color = vec3(0.0);
    
    // Red component
    if (wavelength > 580.0) {
        color.r = (wavelength > 645.0) ? 1.0 : (wavelength - 580.0) / 65.0;
    }
    
    // Green component  
    if (wavelength > 490.0 && wavelength < 645.0) {
        color.g = (wavelength < 580.0) ? 1.0 : 1.0 - (wavelength - 580.0) / 65.0;
        if (wavelength < 510.0) color.g *= (wavelength - 490.0) / 20.0;
    }
    
    // Blue component
    if (wavelength < 510.0) {
        color.b = (wavelength < 440.0) ? 
                  0.3 + 0.7 * (wavelength - 380.0) / 60.0 : 
                  (wavelength < 490.0 ? 1.0 : 1.0 - (wavelength - 490.0) / 20.0);
    }
    
    return color;
}

// =============================================================================
// EXTREME PROXIMITY AND EVENT HORIZON DRAMATIC EFFECTS
// =============================================================================

// Enhanced spacetime curvature visualization with realistic physics
vec3 generateSpacetimeDistortion(vec2 pos, float proximity, float time) {
    if (proximity < 0.05) return vec3(0.0); // Start effects much earlier
    
    // Physics-based distortion scaling - follows 1/sqrt(r-rs) near event horizon
    float physicsScale = 1.0 / sqrt(max(0.01, 1.0 - proximity)); 
    vec2 p = pos * (1.0 + physicsScale * 0.3);
    
    // Tidal stretching effect - vertical stretching increases near event horizon
    float tidalStretch = proximity * proximity * 2.0;
    p.y *= (1.0 + tidalStretch);
    
    // Multi-scale fractal noise representing spacetime metric fluctuations
    float noise = 0.0;
    float amplitude = 1.0;
    float frequency = 1.0;
    
    // More iterations for closer approaches
    int iterations = int(2.0 + proximity * 6.0);
    for (int i = 0; i < iterations && i < 8; i++) {
        // Frame dragging rotation - spacetime rotation increases with proximity
        float frameDragRotation = time * 0.2 * proximity * proximity;
        vec2 rotated = mat2(cos(frameDragRotation), -sin(frameDragRotation), 
                           sin(frameDragRotation), cos(frameDragRotation)) * (p * frequency);
        
        // Spacetime curvature creates wave-like distortions
        noise += sin(rotated.x + cos(rotated.y * 1.3) + time * 0.5) * amplitude;
        amplitude *= 0.6;
        frequency *= 1.8;
    }
    
    // Time dilation effects - colors redshift dramatically near event horizon
    float redshift = 1.0 / sqrt(max(0.01, 1.0 - proximity * 0.9));
    float intensity = proximity * 2.0;  // MUCH MORE VISIBLE
    
    vec3 distortColor = vec3(
        sin(noise + time * 2.0 / redshift) * intensity,           // Redshifted plasma
        sin(noise * 1.5 + time * 1.5 / redshift) * intensity * 0.5, // Dimmed green
        cos(noise * 0.8 + time * 3.0 / redshift) * intensity * 0.2   // Heavily redshifted blue
    ) * redshift * 0.5;
    
    return abs(distortColor);
}

// Extreme chromatic aberration for dramatic color separation
vec3 applyChromaticAberration(vec3 rayDir, vec3 baseColor, float proximity) {
    if (proximity < 0.4) return baseColor; // Only active when very close
    
    // Different wavelengths bend differently in extreme gravitational fields
    float aberration = proximity * proximity * 0.1; // Increases dramatically with proximity
    
    // Separate color channels with different bending amounts
    vec3 redRay = rayDir + vec3(aberration * 0.8, 0.0, 0.0);
    vec3 greenRay = rayDir + vec3(0.0, aberration * 1.0, 0.0);  
    vec3 blueRay = rayDir + vec3(-aberration * 1.2, 0.0, 0.0);
    
    // Sample different rays for each color channel
    float redIntensity = baseColor.r * (1.0 + sin(length(redRay) * 10.0 + time) * proximity * 0.5);
    float greenIntensity = baseColor.g * (1.0 + sin(length(greenRay) * 12.0 + time * 1.2) * proximity * 0.3);
    float blueIntensity = baseColor.b * (1.0 + cos(length(blueRay) * 8.0 + time * 0.8) * proximity * 0.7);
    
    return vec3(redIntensity, greenIntensity, blueIntensity);
}

// Dramatic color shifting from extreme gravitational effects
vec3 applyExtremeColorShift(vec3 color, float proximity, vec3 position) {
    if (proximity < 0.2) return color;
    
    // Extreme gravitational redshift/blueshift effects
    float distortion = proximity * proximity;
    
    // Create rainbow-like prismatic effects from gravitational lensing
    float prism = sin(length(position) * 5.0 + time * 2.0) * distortion;
    
    // Shift colors dramatically - physics-inspired but amplified for drama
    vec3 shifted = color;
    shifted.r *= (1.0 + prism * 0.8);        // Enhanced red shift
    shifted.g *= (1.0 + sin(prism + PI/3) * 0.6); // Green modulation  
    shifted.b *= (1.0 + cos(prism + PI*2/3) * 1.2); // Enhanced blue shift
    
    // Add plasma-like energy colors at extreme proximity
    if (proximity > 0.7) {
        vec3 plasma = vec3(
            sin(time * 3.0 + position.x * 10.0) * 0.5 + 0.5,
            sin(time * 2.5 + position.y * 8.0) * 0.3 + 0.7,
            cos(time * 4.0 + position.z * 12.0) * 0.4 + 0.6
        );
        shifted = mix(shifted, plasma * 2.0, (proximity - 0.7) * 0.8);
    }
    
    return shifted;
}

// Event horizon "firewall" effect - theoretical physics visualization
vec3 generateFirewallEffect(vec2 screenPos, float proximity, float time) {
    if (proximity < 0.8) return vec3(0.0); // Only very close to event horizon
    
    // Generate plasma-like "firewall" patterns
    vec2 p = screenPos * 20.0;
    
    // Turbulent plasma simulation
    float turbulence = 0.0;
    for (int i = 0; i < 3; i++) {
        float scale = pow(2.0, float(i));
        turbulence += sin(p.x * scale + time * 2.0) * cos(p.y * scale + time * 1.5) / scale;
    }
    
    // Firewall intensity increases dramatically near event horizon
    float intensity = pow(proximity - 0.8, 2.0) * 5.0; // Sharp increase at horizon
    
    // Hot plasma colors - white-hot at center, cooler at edges
    vec3 firewall = vec3(
        1.0,                                    // White-hot core
        0.8 + turbulence * 0.3,                // Slightly cooler
        0.3 + max(0.0, turbulence) * 0.4       // Blue-white plasma
    ) * intensity;
    
    return firewall;
}

// Enhanced extreme gravitational lensing with realistic spacetime curvature
vec3 applyExtremeLensing(vec3 rayPos, vec3 rayDir, float proximity) {
    vec3 toBlackHole = normalize(blackholePos - rayPos);
    float distanceToBlackHole = length(blackholePos - rayPos);
    
    if (distanceToBlackHole < 1.5 * schwarzschildRadius) {
        return rayDir; // Don't apply if inside photon sphere
    }
    
    vec3 normalizedDir = normalize(rayDir);
    vec3 normalizedToBlackHole = normalize(toBlackHole);
    
    float cosAngle = dot(normalizedDir, normalizedToBlackHole);
    float impactParameter = distanceToBlackHole * sqrt(max(0.0, 1.0 - cosAngle * cosAngle));
    
    // Physics-based deflection using proper strong-field formula
    // For Schwarzschild: δφ ≈ 4Rs/b + 15πRs²/(4b²) for strong fields
    float rs_over_b = schwarzschildRadius / max(impactParameter, 0.01 * schwarzschildRadius);
    float deflectionAngle = 4.0 * rs_over_b + 15.0 * PI * rs_over_b * rs_over_b * 0.25;
    
    // Add proximity-based amplification for visual drama
    deflectionAngle *= (1.0 + proximity * proximity * 3.0);
    deflectionAngle = min(deflectionAngle, PI * 0.9); // Prevent complete reversal
    
    // Apply the deflection with proper vector mathematics
    vec3 perpendicular = normalize(cross(normalizedDir, normalizedToBlackHole));
    if (length(perpendicular) < 0.01) {
        // Ray pointing directly at black hole - use arbitrary perpendicular
        perpendicular = vec3(0.0, 1.0, 0.0);
    }
    vec3 deflectionAxis = normalize(cross(normalizedDir, perpendicular));
    
    // Frame-dragging effects for rotating black holes
    float frameDragAngle = proximity * proximity * 0.2 * sin(time * 1.5);
    
    vec3 deflected = normalizedDir * cos(deflectionAngle) + 
                     deflectionAxis * sin(deflectionAngle);
    
    // Apply frame dragging rotation around the black hole axis
    float cosFrame = cos(frameDragAngle);
    float sinFrame = sin(frameDragAngle);
    vec3 frameDragged = deflected * cosFrame + 
                       cross(normalizedToBlackHole, deflected) * sinFrame +
                       normalizedToBlackHole * dot(normalizedToBlackHole, deflected) * (1.0 - cosFrame);
    
    return normalize(frameDragged);
}

// New function: Tidal stretching effects (spaghettification)
vec3 applyTidalDistortion(vec3 rayDir, vec3 rayPos, float proximity) {
    if (proximity < 0.1) return rayDir;
    
    // Calculate radial and tangential components relative to black hole
    vec3 radial = normalize(blackholePos - rayPos);
    vec3 tangential = rayDir - dot(rayDir, radial) * radial;
    
    // Tidal force: differential gravity creates stretching
    // F_tidal ∝ (2GM/r³) for radial direction, -(GM/r³) for tangential
    float tidalStrength = proximity * proximity * proximity; // Increases rapidly near horizon
    
    // Stretch vertically (away from radial direction)
    vec3 vertical = normalize(cross(radial, normalize(tangential + vec3(0.001, 0.001, 0.001))));
    float verticalComponent = dot(rayDir, vertical);
    float stretchFactor = 1.0 + tidalStrength * 4.0;
    
    // Compress horizontally
    vec3 horizontal = normalize(cross(radial, vertical));
    float horizontalComponent = dot(rayDir, horizontal);
    float compressFactor = 1.0 / (1.0 + tidalStrength * 2.0);
    
    // Reconstruct ray direction with tidal effects
    vec3 tidalRay = dot(rayDir, radial) * radial + 
                   verticalComponent * stretchFactor * vertical +
                   horizontalComponent * compressFactor * horizontal;
    
    return normalize(tidalRay);
}

// Enhanced time dilation visual effects
vec3 applyTimeDilationEffects(vec3 color, vec2 screenPos, float proximity, float time) {
    if (proximity < 0.05) return color;
    
    // Time dilation factor: proper time / coordinate time = sqrt(1 - rs/r)
    // As proximity → 1, time dilation → ∞
    float dilationFactor = 1.0 / sqrt(max(0.01, 1.0 - proximity * 0.95));
    
    // Slow down time-dependent effects visually
    float dilatedTime = time / dilationFactor;
    
    // Create ripple effects representing spacetime distortion
    float rippleFreq = 5.0 / dilationFactor;  // Slower ripples near horizon
    float ripple1 = sin(length(screenPos) * rippleFreq + dilatedTime * 2.0) * 0.1;
    float ripple2 = cos(screenPos.x * rippleFreq * 1.3 + dilatedTime * 1.5) * 0.05;
    float ripple3 = sin(screenPos.y * rippleFreq * 0.8 + dilatedTime * 2.5) * 0.07;
    
    float totalRipple = (ripple1 + ripple2 + ripple3) * proximity;
    
    // Time dilation creates a "slowing" visual effect with trailing patterns
    vec2 trailOffset = screenPos * proximity * 0.1 * sin(dilatedTime * 0.5);
    float trail = exp(-length(trailOffset) * 10.0) * proximity;
    
    // Color effects: extreme time dilation causes redshift and distortion
    vec3 dilatedColor = color;
    if (proximity > 0.3) {
        // Gravitational redshift becomes extreme
        float redshift = dilationFactor;
        dilatedColor.r *= (1.0 + redshift * 0.3);
        dilatedColor.g /= (1.0 + redshift * 0.2); 
        dilatedColor.b /= (1.0 + redshift * 0.5);
    }
    
    // Add temporal distortion patterns
    vec3 temporalPattern = vec3(
        totalRipple * 0.8,
        totalRipple * 0.4, 
        totalRipple * 0.2
    );
    
    dilatedColor += temporalPattern + vec3(trail * 0.5);
    
    return dilatedColor;
}

// OPTIMIZATION 3: Fast temperature to RGB using lookup approximation
vec3 fastTemperatureToRGB(float temperature) {
    // Fast blackbody approximation without expensive loops
    if (temperature < 1000.0) return vec3(0.0);
    
    // Use Wien's law for dominant wavelength
    float peakWavelength = 2.898e6 / temperature; // nm
    vec3 peakColor = fastWavelengthToRGB(peakWavelength);
    
    // Temperature-based color adjustment
    if (temperature < 3000.0) {
        return vec3(1.0, 0.3, 0.1) * (temperature / 3000.0);
    } else if (temperature < 5000.0) {
        float t = (temperature - 3000.0) / 2000.0;
        return mix(vec3(1.0, 0.3, 0.1), vec3(1.0, 0.8, 0.4), t);
    } else if (temperature < 6500.0) {
        float t = (temperature - 5000.0) / 1500.0;
        return mix(vec3(1.0, 0.8, 0.4), vec3(1.0, 1.0, 0.8), t);
    } else {
        float t = min(1.0, (temperature - 6500.0) / 3500.0);
        return mix(vec3(1.0, 1.0, 0.8), vec3(0.8, 0.9, 1.0), t);
    }
}

// OPTIMIZATION 4: LOD-based star generation
vec3 generateStarsLOD(vec3 direction, vec3 observerPos, int lod) {
    vec3 d = normalize(direction);
    float hash = fract(sin(dot(d.xy, vec2(12.9898, 78.233)) + d.z * 37.719) * 43758.5453);
    
    // Adjust star density based on LOD
    float starThreshold = 0.995 + (1.0 - performanceScale) * 0.003; // Fewer stars at low performance
    
    if (hash > starThreshold) {
        float starTemp = 8000.0 + hash * 10000.0;
        vec3 starColor = (enableSpectralRendering == 1) ? 
                        fastTemperatureToRGB(starTemp) : vec3(1.0, 1.0, 0.9);
        
        // Simplified lensing calculation
        float distanceToBlackHole = length(observerPos - blackholePos);
        float lensingFactor = 1.0 + schwarzschildRadius * 0.5 / (distanceToBlackHole * distanceToBlackHole);
        
        return starColor * lensingFactor * starBrightness;
    } else if (hash > starThreshold - 0.002) {
        float starTemp = 5000.0 + hash * 3000.0;
        vec3 starColor = (enableSpectralRendering == 1) ? 
                        fastTemperatureToRGB(starTemp) : vec3(0.8, 0.8, 0.6);
        
        float distanceToBlackHole = length(observerPos - blackholePos);
        float lensingFactor = 1.0 + schwarzschildRadius * 0.25 / (distanceToBlackHole * distanceToBlackHole);
        
        return starColor * lensingFactor * starBrightness;
    } else if (hash > starThreshold - 0.003 && lod <= 1) {
        // Only show dim stars at high quality
        float starTemp = 3000.0 + hash * 2000.0;
        return ((enableSpectralRendering == 1) ? fastTemperatureToRGB(starTemp) : vec3(0.6, 0.5, 0.4)) * starBrightness * 0.8;
    } else {
        return vec3(0.02, 0.02, 0.05) + vec3(ambientLightLevel);
    }
}

// OPTIMIZATION 5: Simplified accretion disk with early exit
vec3 sampleAccretionDiskOptimized(vec3 rayPos, vec3 rayDir, vec3 observerPos) {
    vec3 diskCenter = blackholePos;
    vec3 diskNormal = vec3(0, 0, 1);
    
    // Early exit if ray doesn't intersect disk plane
    float denom = dot(rayDir, diskNormal);
    if (abs(denom) < 0.001) return vec3(0.0);
    
    float t = dot(diskCenter - rayPos, diskNormal) / denom;
    if (t < 0.0) return vec3(0.0);
    
    vec3 hitPoint = rayPos + t * rayDir;
    float distFromCenter = length(hitPoint.xy - diskCenter.xy);
    
    float innerRadius = 3.0 * schwarzschildRadius;
    float outerRadius = 12.0 * schwarzschildRadius;
    float diskThickness = 0.2 * schwarzschildRadius;
    
    // Early exit if outside disk bounds
    if (distFromCenter < innerRadius || distFromCenter > outerRadius || 
        abs(hitPoint.z - diskCenter.z) > diskThickness) {
        return vec3(0.0);
    }
    
    // Simplified temperature calculation
    float localTemp = diskTemperature * pow(innerRadius / distFromCenter, 0.75);
    
    // Simplified relativistic correction
    if (distFromCenter < 6.0 * schwarzschildRadius) {
        localTemp *= sqrt(1.0 - schwarzschildRadius / distFromCenter);
    }
    
    vec3 diskColor = (enableSpectralRendering == 1) ? 
                    fastTemperatureToRGB(localTemp) : 
                    vec3(0.3 + 0.7 * (outerRadius / distFromCenter - 1.0),
                        0.4 + 0.6 * (outerRadius / distFromCenter),
                        0.6 + 0.4 * (outerRadius / distFromCenter));
    
    // Simplified Doppler and redshift effects
    if (enableSpectralRendering == 1) {
        float orbitalVelocity = sqrt(schwarzschildRadius / (2.0 * distFromCenter));
        float angle = atan(hitPoint.y - diskCenter.y, hitPoint.x - diskCenter.x);
        
        vec3 velocityDir = vec3(-sin(angle), cos(angle), 0.0);
        vec3 photonDir = normalize(observerPos - hitPoint);
        float dopplerShift = 1.0 + orbitalVelocity * dot(velocityDir, photonDir) * 3.33e-9; // c^-1 simplified
        
        float gravRedshift = gravitationalRedshift * sqrt(1.0 - schwarzschildRadius / distFromCenter) /
                            sqrt(1.0 - schwarzschildRadius / length(observerPos - diskCenter));
        
        float totalShift = dopplerShift * gravRedshift;
        if (abs(totalShift - 1.0) > 0.01) {
            diskColor = fastTemperatureToRGB(localTemp * totalShift);
        }
    }
    
    // Simplified turbulence and falloff
    float rotation = atan(hitPoint.y - diskCenter.y, hitPoint.x - diskCenter.x) + 
                    time * sqrt(schwarzschildRadius / distFromCenter) * 0.1;
    float turbulence = 0.8 + 0.3 * sin(rotation * 8.0 + time * 2.0);
    
    float heightFactor = exp(-hitPoint.z * hitPoint.z / (diskThickness * diskThickness));
    float intensityFalloff = 1.0 / (1.0 + pow(distFromCenter / schwarzschildRadius, 1.5));
    
    return diskColor * turbulence * heightFactor * intensityFalloff * diskBrightness;
}

// OPTIMIZATION 6: Simplified gravitational lensing
vec3 applyGravitationalLensingFast(vec3 rayPos, vec3 rayDir) {
    vec3 toBlackHole = blackholePos - rayPos;
    float distanceToBlackHole = length(toBlackHole);
    
    if (distanceToBlackHole < 2.0 * schwarzschildRadius) {
        return rayDir;
    }
    
    vec3 normalizedDir = normalize(rayDir);
    vec3 normalizedToBlackHole = normalize(toBlackHole);
    
    // Simplified deflection calculation
    float cosAngle = dot(normalizedDir, normalizedToBlackHole);
    float impactParameter = distanceToBlackHole * sqrt(1.0 - cosAngle * cosAngle);
    
    float deflectionAngle = 2.0 * schwarzschildRadius / max(impactParameter, 0.1 * schwarzschildRadius);
    deflectionAngle = min(deflectionAngle, PI * 0.25);
    
    // Simplified deflection application
    vec3 perpendicular = normalize(cross(normalizedDir, cross(normalizedDir, normalizedToBlackHole)));
    return normalize(normalizedDir + perpendicular * sin(deflectionAngle));
}

// OPTIMIZATION 7: ULTRA-FAST ray visualization with precomputed patterns
vec3 generateLightRayVisualizationLOD(vec2 screenCoords, vec3 cameraPos, vec3 rayDir) {
    if (showLightRays == 0) return vec3(0.0);
    
    // EMERGENCY EXIT: Turn off entirely for ANY performance issues
    if (performanceScale < 0.7) {
        return vec3(0.0);
    }
    
    // Distance-based ray density - fewer rays when farther away
    float cameraDistance = length(cameraPos - blackholePos);
    
    // CRITICAL: Turn off for distant views (major GPU saver)
    if (cameraDistance > 30.0 * schwarzschildRadius) {
        return vec3(0.0);
    }
    
    // ULTRA-SIMPLE light ring effect only
    vec3 toBlackHole = normalize(blackholePos - cameraPos);
    float alignment = dot(rayDir, toBlackHole);
    
    // Single simple calculation - no complex math
    if (alignment > 0.75) {
        float intensity = (alignment - 0.75) * 4.0; // Precomputed (1/0.25)
        vec3 rayColor = vec3(0.8, 0.6, 0.3) * intensity * performanceScale * 0.2;
        return rayColor;
    }
    
    return vec3(0.0);
}

// OPTIMIZATION 8: Simplified relativistic jets
vec3 calculateJetEmissionOptimized(vec3 rayPos, vec3 rayDir, vec3 observerPos) {
    if (showJets == 0) return vec3(0.0);
    
    // EMERGENCY EXIT: Turn off entirely for ANY performance issues
    if (performanceScale < 0.8) {
        return vec3(0.0);
    }
    
    // ULTRA-SIMPLE jet visualization - single calculation only
    vec3 jetDirection = normalize(jetAxis);
    float jetAlignment = dot(rayDir, jetDirection);
    
    // Single simple jet glow - no complex sampling
    if (jetAlignment > 0.85) {
        float intensity = (jetAlignment - 0.85) * 6.67; // Precomputed (1/0.15)
        vec3 jetColor = vec3(0.6, 0.8, 1.0) * intensity * performanceScale * 0.3;
        return jetColor;
    }
    
    return vec3(0.0);
}

void main() {
    // Debug modes
    if (debugMode == 1) {
        vec2 screenPos = (gl_FragCoord.xy / resolution) * 2.0 - 1.0;
        vec3 color = vec3(0.1, 0.1, 0.2);
        if (abs(screenPos.x) < 0.05 && abs(screenPos.y) < 0.05) color = vec3(1.0);
        color.r += abs(screenPos.x) * 0.5;
        color.g += abs(screenPos.y) * 0.5;
        FragColor = vec4(color, 1.0);
        return;
    } else if (debugMode == 2) {
        FragColor = vec4(0.0, 1.0, 0.0, 1.0);
        return;
    }
    
    // EMERGENCY MODE: Ultra-minimal rendering for critical performance
    // CRITICAL FIX: Lower threshold from 0.2 to 0.4 for earlier intervention
    if (performanceScale < 0.4) {
        vec2 screenPos = (gl_FragCoord.xy / resolution) * 2.0 - 1.0;
        float dist = length(screenPos);
        
        // ULTRA-SIMPLE rendering - no complex math at all
        if (dist < 0.08) {
            FragColor = vec4(0.0, 0.0, 0.0, 1.0); // Event horizon - pure black
        } else if (dist < 0.25) {
            // Simple gradient - no complex calculations
            float fade = (dist - 0.08) * 5.88; // Precomputed 1/(0.25-0.08)
            FragColor = vec4(fade * 0.08, fade * 0.04, 0.0, 1.0); // Simple orange glow
        } else {
            // ULTRA-SIMPLE star field - single hash, minimal branching
            float hash = fract(sin(dot(gl_FragCoord.xy, vec2(12.9898, 78.233))) * 43758.5453);
            // Single branch decision
            if (hash > 0.997) {
                FragColor = vec4(0.6, 0.6, 0.8, 1.0); // Single star brightness
            } else {
                FragColor = vec4(0.01, 0.01, 0.03, 1.0); // Dark space
            }
        }
        return;
    }
    
    // SEVERE PERFORMANCE MODE: Still very simple but slightly better
    if (performanceScale < 0.6) {
        vec2 screenPos = (gl_FragCoord.xy / resolution) * 2.0 - 1.0;
        float dist = length(screenPos);
        
        if (dist < 0.1) {
            FragColor = vec4(0.0, 0.0, 0.0, 1.0); // Event horizon
        } else if (dist < 0.3) {
            float fade = (dist - 0.1) * 5.0; // Precomputed
            FragColor = vec4(fade * 0.1, fade * 0.05, fade * 0.02, 1.0); // Simple glow
        } else {
            // Simple stars with two levels only
            float hash = fract(sin(dot(gl_FragCoord.xy, vec2(12.9898, 78.233))) * 43758.5453);
            if (hash > 0.996) {
                FragColor = vec4(0.8, 0.8, 1.0, 1.0);
            } else {
                FragColor = vec4(0.02, 0.02, 0.05, 1.0);
            }
        }
        return;
    }
    
    // Normal rendering with optimizations
    vec2 screenCoords = (gl_FragCoord.xy / resolution) * 2.0 - 1.0;
    
    float aspectRatio = resolution.x / resolution.y;
    float fov = 45.0 * PI / 180.0;
    float tanHalfFov = tan(fov * 0.5);
    
    screenCoords.x *= aspectRatio;
    
    vec3 rayDir = normalize(vec3(
        screenCoords.x * tanHalfFov,
        screenCoords.y * tanHalfFov,
        -1.0
    ));
    
    // Event horizon check (simplified)
    vec3 oc = cameraPos - blackholePos;
    float a = dot(rayDir, rayDir);
    float b = 2.0 * dot(oc, rayDir);
    float c = dot(oc, oc) - schwarzschildRadius * schwarzschildRadius;
    float discriminant = b * b - 4.0 * a * c;
    
    // Use proximity factor calculated in C++ and passed as uniform
    // (avoid recalculating to prevent inconsistencies)
    float currentCameraDistance = length(cameraPos - blackholePos);
    
    // Apply enhanced gravitational lensing (normal or extreme based on proximity)
    vec3 lensedRayDir;
    if (enableDramaticEffects == 1 && proximityFactor > 0.05) {
        // Apply tidal distortion first (affects the view itself)
        vec3 tidalDistortedRayDir = applyTidalDistortion(rayDir, cameraPos, proximityFactor);
        // Then apply gravitational lensing to the distorted rays
        lensedRayDir = applyExtremeLensing(cameraPos, tidalDistortedRayDir, proximityFactor);
    } else if (performanceScale > 0.5) {
        lensedRayDir = applyGravitationalLensingFast(cameraPos, rayDir);
    } else {
        lensedRayDir = rayDir;
    }
    
    // Start with optimized background (simplified at low performance)
    vec3 color = (performanceScale > 0.3) ? 
                 generateStarsLOD(lensedRayDir, cameraPos, lodLevel) :
                 vec3(0.02, 0.02, 0.05); // Simple space color
    
    // Add optimized light ray visualization with early exit
    if (showLightRays == 1) {
        // Use already calculated camera distance to avoid duplicate calculations
        
        // CRITICAL PERFORMANCE: Turn off light rays entirely for severe issues
        if (performanceScale < 0.15) {
            // Skip light rays completely
        }
        // Super aggressive culling for very distant views or poor performance
        else if (currentCameraDistance > 80.0 * schwarzschildRadius && performanceScale < 0.5) {
            // Show minimal "fake" light ring only
            vec3 toBlackHole = normalize(blackholePos - cameraPos);
            float alignment = dot(lensedRayDir, toBlackHole);
            if (alignment > 0.8) {
                color += vec3(0.1, 0.08, 0.05) * (alignment - 0.8) * 5.0 * performanceScale;
            }
        }
        // Aggressive culling for moderate performance issues
        else if (performanceScale < 0.3 || currentCameraDistance > 50.0 * schwarzschildRadius) {
            // Ultra-simple light rays only
            vec3 toBlackHole = normalize(blackholePos - cameraPos);
            float alignment = dot(lensedRayDir, toBlackHole);
            if (alignment > 0.7) {
                color += vec3(0.2, 0.15, 0.1) * (alignment - 0.7) * 3.0 * performanceScale;
            }
        }
        else {
            color += generateLightRayVisualizationLOD(screenCoords, cameraPos, rayDir);
        }
    }
    
    // Add optimized jet emission
    color += calculateJetEmissionOptimized(cameraPos, lensedRayDir, cameraPos);
    
    // =========================================================================
    // APPLY PROGRESSIVE PHYSICS-BASED PROXIMITY EFFECTS
    // =========================================================================
    
    if (enableDramaticEffects == 1 && proximityFactor > 0.01) {
        
        // Stage 1: VISIBLE spacetime curvature (starts very early for extended experience)
        if (proximityFactor > 0.01) {
            vec3 spacetimeEffect = generateSpacetimeDistortion(screenCoords, proximityFactor, time);
            color += spacetimeEffect * 5.0 * smoothstep(0.01, 0.1, proximityFactor); // MUCH MORE VISIBLE
            
            // Add subtle color warping even at large distances
            if (proximityFactor > 0.02) {
                color.r += sin(time * 2.0 + screenCoords.x * 10.0) * proximityFactor * 0.3;
                color.b += cos(time * 1.5 + screenCoords.y * 8.0) * proximityFactor * 0.2;
            }
        }
        
        // Stage 2: DRAMATIC time dilation redshift effects 
        if (proximityFactor > 0.05) {
            float redshiftFactor = 1.0 / sqrt(max(0.05, 1.0 - proximityFactor * 0.9));
            color.r *= (1.0 + (redshiftFactor - 1.0) * 1.5); // MUCH more red
            color.g *= (1.0 + (redshiftFactor - 1.0) * 0.8); // More green shift
            color.b /= (1.0 + (redshiftFactor - 1.0) * 2.0);  // Much less blue
            
            // Add pulsing red glow for extreme redshift
            if (proximityFactor > 0.3) {
                float pulse = sin(time * 3.0) * 0.5 + 0.5;
                color += vec3(pulse * proximityFactor * 0.8, 0.0, 0.0);
            }
        }
        
        // Stage 3: Color shifting and chromatic aberration (medium proximity)
        if (proximityFactor > 0.1) {
            color = applyExtremeColorShift(color, proximityFactor, cameraPos);
            color = applyChromaticAberration(lensedRayDir, color, proximityFactor);
        }
        
        // Stage 4: Intense spacetime warping (close proximity)
        if (proximityFactor > 0.3) {
            // Add swirling patterns representing frame dragging
            vec2 swirl = screenCoords;
            float swirlAngle = proximityFactor * proximityFactor * time * 2.0;
            mat2 swirlMatrix = mat2(cos(swirlAngle), -sin(swirlAngle), 
                                   sin(swirlAngle), cos(swirlAngle));
            vec2 swirled = swirlMatrix * swirl;
            
            float swirlIntensity = proximityFactor * 1.5; // MUCH MORE INTENSE
            vec3 swirlColor = vec3(
                sin(swirled.x * 10.0 + time) * swirlIntensity,
                sin(swirled.y * 10.0 + time * 1.5) * swirlIntensity * 0.6,
                cos(length(swirled) * 15.0 + time * 2.0) * swirlIntensity * 0.4
            );
            color += abs(swirlColor);
            
            // Add bright flashing for extreme warping
            if (proximityFactor > 0.5) {
                float flash = sin(time * 10.0) * sin(time * 7.0) * proximityFactor;
                color += vec3(abs(flash) * 0.5, abs(flash) * 0.3, abs(flash) * 0.1);
            }
        }
        
        // Stage 5: EXTREME gravitational effects (very close)
        if (proximityFactor > 0.5) {
            // Add dramatic energy distortion patterns across the screen
            float energyPattern = sin(screenCoords.x * 20.0 + time * 3.0) * 
                                 cos(screenCoords.y * 15.0 + time * 2.0) * 
                                 proximityFactor * 1.0;  // 10x MORE INTENSE
            color += vec3(abs(energyPattern), abs(energyPattern) * 0.5, abs(energyPattern) * 0.2);
            
            // Add screen-wide electric field effects
            float electricField = sin(screenCoords.x * 50.0 + time * 8.0) *
                                 sin(screenCoords.y * 30.0 + time * 5.0) *
                                 proximityFactor * 0.8;
            color += vec3(0.2, 0.4, 0.8) * abs(electricField);
        }
        
        // Stage 6: DRAMATIC Event horizon firewall effect (extreme proximity)
        if (proximityFactor > 0.8) {
            vec3 firewall = generateFirewallEffect(screenCoords, proximityFactor, time);
            float firewallMix = smoothstep(0.8, 1.0, proximityFactor);
            color = mix(color, firewall * 3.0, firewallMix); // MUCH MORE INTENSE
            
            // Add chaotic screen distortion at event horizon
            if (proximityFactor > 0.95) {
                vec2 distortedCoords = screenCoords + vec2(
                    sin(time * 15.0 + screenCoords.y * 20.0) * 0.1,
                    cos(time * 12.0 + screenCoords.x * 18.0) * 0.1
                ) * proximityFactor;
                
                float chaos = sin(distortedCoords.x * 100.0) * cos(distortedCoords.y * 80.0);
                color += vec3(abs(chaos) * 2.0, abs(chaos) * 1.5, abs(chaos) * 0.5);
            }
        }
        
        // Apply time dilation effects throughout (affects visual perception of time)
        color = applyTimeDilationEffects(color, screenCoords, proximityFactor, time);
    }
    
    // Enhanced event horizon rendering with physics-based approach
    if (discriminant >= 0.0) {
        float t = (-b - sqrt(discriminant)) / (2.0 * a);
        if (t > 0.0) {
            if (enableDramaticEffects == 1) {
                // Distance to event horizon surface for smooth transition
                vec3 rayHitPoint = cameraPos + rayDir * t;
                float distanceToSurface = abs(length(rayHitPoint - blackholePos) - schwarzschildRadius);
                float surfaceProximity = 1.0 - clamp(distanceToSurface / (schwarzschildRadius * 0.1), 0.0, 1.0);
                
                if (proximityFactor > 0.95) {
                    // Ultra-close to event horizon - show theoretical physics effects
                    vec2 vortexPos = screenCoords * (1.0 + surfaceProximity * 3.0);
                    float vortexAngle = atan(vortexPos.y, vortexPos.x) + time * 8.0 * surfaceProximity;
                    float vortexRadius = length(vortexPos);
                    
                    // Information paradox visualization - quantum fluctuations at the horizon
                    float quantumNoise = fract(sin(dot(vortexPos * 100.0, vec2(12.9898, 78.233))) * 43758.5453);
                    quantumNoise = smoothstep(0.3, 0.7, quantumNoise) * surfaceProximity;
                    
                    vec3 vortexColor = vec3(
                        sin(vortexAngle * 4.0 + vortexRadius * 15.0 + quantumNoise * 10.0) * 0.6 + 0.4,
                        cos(vortexAngle * 3.0 + time * 4.0 + quantumNoise * 8.0) * 0.4 + 0.6,
                        sin(vortexAngle * 5.0 - time * 3.0 + quantumNoise * 12.0) * 0.5 + 0.5
                    );
                    
                    // Hawking radiation-like glow at the very edge
                    float hawkingGlow = exp(-vortexRadius * 5.0) * surfaceProximity * 0.8;
                    vortexColor += vec3(hawkingGlow * 1.2, hawkingGlow * 0.8, hawkingGlow * 0.4);
                    
                    color = mix(vec3(0.05, 0.02, 0.0), vortexColor, surfaceProximity * 0.8);
                } else {
                    // Gradual fade to black with some remaining light effects
                    float fadeOut = 1.0 - smoothstep(0.8, 1.0, proximityFactor);
                    color *= fadeOut;
                    
                    // Add some residual "last light" effects
                    if (fadeOut > 0.1) {
                        vec3 lastLight = vec3(0.3, 0.15, 0.05) * fadeOut * fadeOut;
                        color += lastLight;
                    }
                }
            } else {
                color = vec3(0.0); // Traditional black hole interior
            }
        }
    } else {
        // Sample optimized accretion disk
        vec3 diskContribution = sampleAccretionDiskOptimized(cameraPos, lensedRayDir, cameraPos);
        
        if (length(diskContribution) > 0.01) {
            color = diskContribution;
        } else {
            // Simplified photon sphere glow
            vec3 toBlackHole = normalize(blackholePos - cameraPos);
            float alignment = dot(lensedRayDir, toBlackHole);
            
            if (alignment > 0.85) {
                float glow = pow((alignment - 0.85) / 0.15, 2.0);
                vec3 glowColor = vec3(1.5, 0.8, 0.2);
                
                // Add jet enhancement
                if (showJets == 1) {
                    vec3 jetDirection = normalize(jetAxis);
                    float jetAlignment = dot(lensedRayDir, jetDirection);
                    if (jetAlignment > 0.9) {
                        glowColor += vec3(0.6, 0.8, 1.2) * jetLorentzFactor * 0.1;
                    }
                }
                
                color += glowColor * glow * photonSphereGlow;
            }
        }
    }
    
    // Add perpetual photon sphere glow for enhanced realism (always visible)
    if (enablePhotonSphereGlow == 1) {
        float distanceFromBH = length(cameraPos - blackholePos);
        float photonSphereRadius = 1.5 * schwarzschildRadius;
        
        // Create a subtle but persistent glow near the photon sphere
        if (distanceFromBH > schwarzschildRadius && distanceFromBH < photonSphereRadius * 4.0) {
            float glowIntensity = photonSphereGlow * 0.5 * (1.0 / (1.0 + pow(distanceFromBH / photonSphereRadius, 2.0)));
            vec3 persistentGlow = vec3(0.2, 0.4, 0.8) * glowIntensity; // Blue-ish glow
            color += persistentGlow;
        }
    }
    
    // Final dramatic enhancement for extreme proximity
    if (enableDramaticEffects == 1 && proximityFactor > 0.5) {
        // Add energy distortion around the entire screen when very close
        float screenDist = length(screenCoords);
        float energyPulse = sin(time * 4.0 + screenDist * 8.0) * proximityFactor * 0.1;
        color *= (1.0 + energyPulse);
        
        // Boost overall intensity for dramatic effect
        color *= (1.0 + proximityFactor * 0.3);
    }
    
    // Ensure colors don't exceed reasonable bounds
    color = clamp(color, 0.0, 3.0);
    
    FragColor = vec4(color, 1.0);
}
