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

// Enhanced proximity effect uniforms
uniform float proximityFactor;   // 0.0=far, 1.0=at event horizon
uniform int enableDramaticEffects; // Toggle for enhanced visual drama

const float PI = 3.14159265359;

// Physical constants for spectral rendering
const float PLANCK_CONSTANT = 6.62607015e-34;
const float BOLTZMANN_CONSTANT = 1.380649e-23;
const float SPEED_OF_LIGHT = 299792458.0;
const float WIEN_CONSTANT = 2.897771955e-3;

// Visible spectrum wavelengths (nm)
const float WAVELENGTH_RED = 700.0;
const float WAVELENGTH_GREEN = 550.0;
const float WAVELENGTH_BLUE = 450.0;

// Forward declaration
float planckFunction(float temperature, float wavelength);

// Convert wavelength to RGB color
vec3 wavelengthToRGB(float wavelength) {
    float r = 0.0, g = 0.0, b = 0.0;
    
    if (wavelength >= 580.0 && wavelength < 645.0) {
        r = (wavelength - 580.0) / (645.0 - 580.0);
    } else if (wavelength >= 645.0 && wavelength <= 750.0) {
        r = 1.0;
    } else if (wavelength > 750.0) {
        r = 1.0 - 0.8 * (wavelength - 750.0) / (800.0 - 750.0);
    }
    
    if (wavelength >= 490.0 && wavelength < 510.0) {
        g = (wavelength - 490.0) / (510.0 - 490.0);
    } else if (wavelength >= 510.0 && wavelength < 580.0) {
        g = 1.0;
    } else if (wavelength >= 580.0 && wavelength <= 645.0) {
        g = 1.0 - (wavelength - 580.0) / (645.0 - 580.0);
    }
    
    if (wavelength >= 380.0 && wavelength < 440.0) {
        b = 0.3 + 0.7 * (wavelength - 380.0) / (440.0 - 380.0);
    } else if (wavelength >= 440.0 && wavelength <= 490.0) {
        b = 1.0;
    } else if (wavelength > 490.0 && wavelength < 510.0) {
        b = 1.0 - (wavelength - 490.0) / (510.0 - 490.0);
    }
    
    // Intensity falloff at spectrum edges
    float intensity = 1.0;
    if (wavelength < 420.0) {
        intensity = 0.3 + 0.7 * (wavelength - 380.0) / (420.0 - 380.0);
    } else if (wavelength > 700.0) {
        intensity = 0.3 + 0.7 * (750.0 - wavelength) / (750.0 - 700.0);
    }
    
    return vec3(r * intensity, g * intensity, b * intensity);
}

// Blackbody temperature to RGB conversion
vec3 temperatureToRGB(float temperature) {
    // Wien's displacement law for peak wavelength
    float peakWavelength = WIEN_CONSTANT / temperature * 1e9; // Convert to nm
    
    // Sample blackbody spectrum and integrate
    vec3 color = vec3(0.0);
    for (int i = 0; i < 32; ++i) {
        float lambda = 380.0 + float(i) * (750.0 - 380.0) / 32.0;
        float planck = planckFunction(temperature, lambda * 1e-9);
        vec3 spectral = wavelengthToRGB(lambda);
        color += spectral * planck * (750.0 - 380.0) / 32.0;
    }
    
    // Normalize
    float maxComponent = max(max(color.r, color.g), color.b);
    if (maxComponent > 1.0) {
        color /= maxComponent;
    }
    
    return color;
}

// Planck blackbody function
float planckFunction(float temperature, float wavelength) {
    float lambda5 = pow(wavelength, 5.0);
    float expTerm = exp(PLANCK_CONSTANT * SPEED_OF_LIGHT / (wavelength * BOLTZMANN_CONSTANT * temperature));
    
    if (expTerm > 1e10) return 0.0; // Avoid overflow
    
    return (2.0 * PLANCK_CONSTANT * SPEED_OF_LIGHT * SPEED_OF_LIGHT / lambda5) / (expTerm - 1.0);
}

// Enhanced stellar spectrum generation
vec3 generateStars(vec3 direction, vec3 observerPos) {
    vec3 d = normalize(direction);
    float hash = fract(sin(dot(d.xy, vec2(12.9898, 78.233)) + d.z * 37.719) * 43758.5453);
    
    if (hash > 0.998) {
        // Hot blue-white star
        float starTemp = 8000.0 + hash * 10000.0;
        vec3 starColor = (enableSpectralRendering == 1) ? temperatureToRGB(starTemp) : vec3(1.0, 1.0, 0.9);
        
        // Apply gravitational lensing brightness
        float distanceToBlackHole = length(observerPos - blackholePos);
        float lensingFactor = 1.0 + schwarzschildRadius / (distanceToBlackHole * distanceToBlackHole);
        
        return starColor * lensingFactor * starBrightness;
    } else if (hash > 0.995) {
        // Sun-like yellow star
        float starTemp = 5000.0 + hash * 3000.0;
        vec3 starColor = (enableSpectralRendering == 1) ? temperatureToRGB(starTemp) : vec3(0.8, 0.8, 0.6);
        
        float distanceToBlackHole = length(observerPos - blackholePos);
        float lensingFactor = 1.0 + 0.5 * schwarzschildRadius / (distanceToBlackHole * distanceToBlackHole);
        
        return starColor * lensingFactor * starBrightness;
    } else if (hash > 0.99) {
        // Red dwarf
        float starTemp = 3000.0 + hash * 2000.0;
        vec3 starColor = (enableSpectralRendering == 1) ? temperatureToRGB(starTemp) : vec3(0.6, 0.5, 0.4);
        
        return starColor * starBrightness * 0.8; // Slightly dimmer than main stars
    } else {
        // Deep space background with cosmic microwave background and ambient light
        float cmbTemp = 2.7; // Kelvin
        vec3 cmbColor = (enableSpectralRendering == 1) ? temperatureToRGB(cmbTemp) * 0.001 : vec3(0.02, 0.02, 0.05);
        return cmbColor + vec3(ambientLightLevel); // Add ambient light for better visibility
    }
}

// Enhanced accretion disk with realistic physics
vec3 sampleAccretionDisk(vec3 rayPos, vec3 rayDir, vec3 observerPos) {
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
        
        // Enhanced physics-based temperature calculation
        float baseDiskTemp = diskTemperature;
        float localTemp = baseDiskTemp * pow(innerRadius / distFromCenter, 0.75);
        
        // Relativistic correction near ISCO
        if (distFromCenter < 6.0 * schwarzschildRadius) {
            float relFactor = sqrt(1.0 - schwarzschildRadius / distFromCenter);
            localTemp *= relFactor;
        }
        
        vec3 diskColor;
        if (enableSpectralRendering == 1) {
            // Use physics-based blackbody spectrum
            diskColor = temperatureToRGB(localTemp);
        } else {
            // Fallback temperature gradient
            float temp = outerRadius / distFromCenter;
            diskColor = vec3(
                0.3 + 0.7 * (1.0 - temp),
                0.4 + 0.6 * temp,
                0.6 + 0.4 * temp
            );
        }
        
        // Orbital motion and Doppler effects
        float angle = atan(hitPoint.y - diskCenter.y, hitPoint.x - diskCenter.x);
        float orbitalVelocity = sqrt(schwarzschildRadius / (2.0 * distFromCenter));
        
        // Calculate Doppler shift
        vec3 velocityDir = vec3(-sin(angle), cos(angle), 0.0);
        vec3 photonDir = normalize(observerPos - hitPoint);
        float dopplerShift = 1.0 + orbitalVelocity * dot(velocityDir, photonDir) / SPEED_OF_LIGHT;
        
        // Apply gravitational redshift
        float gravRedshift = gravitationalRedshift * sqrt(1.0 - schwarzschildRadius / distFromCenter) /
                            sqrt(1.0 - schwarzschildRadius / length(observerPos - diskCenter));
        
        float totalShift = dopplerShift * gravRedshift;
        
        // Shift color based on redshift/blueshift
        if (enableSpectralRendering == 1 && totalShift != 1.0) {
            float shiftedTemp = localTemp * totalShift;
            diskColor = temperatureToRGB(shiftedTemp);
        }
        
        // Dynamic turbulence and rotation
        float rotation = angle + time * sqrt(schwarzschildRadius / distFromCenter) * 0.1;
        float turbulence = 0.3 * sin(rotation * 8.0 + time * 2.0) * 
                          sin(distFromCenter * 0.5 + time * 3.0) *
                          cos(hitPoint.z * 5.0);
        
        // Height-based intensity falloff
        float heightFactor = exp(-hitPoint.z * hitPoint.z / (diskThickness * diskThickness));
        
        // Distance-based intensity
        float intensityFalloff = 1.0 / (1.0 + pow(distFromCenter / schwarzschildRadius, 1.5));
        
        // Apply enhanced brightness for more realistic light display
        return diskColor * (0.8 + turbulence) * heightFactor * intensityFalloff * diskBrightness;
    }
    
    return vec3(0.0);
}

// Enhanced gravitational lensing with proper geodesic approximation
vec3 applyGravitationalLensing(vec3 rayPos, vec3 rayDir) {
    vec3 toBlackHole = blackholePos - rayPos;
    float distanceToBlackHole = length(toBlackHole);
    
    // Don't apply lensing if too close to avoid numerical issues
    if (distanceToBlackHole < 2.0 * schwarzschildRadius) {
        return rayDir;
    }
    
    // Calculate impact parameter (perpendicular distance from ray to black hole)
    vec3 normalizedDir = normalize(rayDir);
    vec3 normalizedToBlackHole = normalize(toBlackHole);
    
    // More accurate deflection using the proper Einstein deflection formula
    float impactParameter = distanceToBlackHole * sin(acos(dot(normalizedDir, normalizedToBlackHole)));
    
    // Einstein deflection angle: α ≈ 4GM/(c²b) for weak field
    // In geometric units where GM=rs/2: α ≈ 2rs/b
    float deflectionAngle = 2.0 * schwarzschildRadius / max(impactParameter, 0.1 * schwarzschildRadius);
    deflectionAngle = min(deflectionAngle, PI / 4.0); // Clamp for stability
    
    // Apply deflection in the plane containing ray and black hole
    vec3 perpendicular = normalize(cross(normalizedDir, cross(normalizedDir, normalizedToBlackHole)));
    vec3 deflectedDir = normalize(normalizedDir + perpendicular * sin(deflectionAngle));
    
    return deflectedDir;
}

// Trace a light ray through curved spacetime with multiple deflections
struct RayTraceResult {
    vec3 finalDirection;
    vec3 accumulatedColor;
    float totalDeflection;
    bool hitEventHorizon;
};

RayTraceResult traceGeodesicRay(vec3 startPos, vec3 startDir, int maxSteps) {
    RayTraceResult result;
    result.finalDirection = startDir;
    result.accumulatedColor = vec3(0.0);
    result.totalDeflection = 0.0;
    result.hitEventHorizon = false;
    
    vec3 currentPos = startPos;
    vec3 currentDir = normalize(startDir);
    float stepSize = 0.5 * schwarzschildRadius;
    
    for (int i = 0; i < maxSteps; i++) {
        vec3 prevDir = currentDir;
        
        // Apply gravitational lensing at each step
        currentDir = applyGravitationalLensing(currentPos, currentDir);
        
        // Calculate deflection this step
        float stepDeflection = acos(clamp(dot(prevDir, currentDir), -1.0, 1.0));
        result.totalDeflection += stepDeflection;
        
        // Move along the ray
        currentPos += currentDir * stepSize;
        
        // Check if we hit the event horizon
        float distToBlackHole = length(currentPos - blackholePos);
        if (distToBlackHole < schwarzschildRadius) {
            result.hitEventHorizon = true;
            break;
        }
        
        // Sample any light sources along the way (accretion disk, stars)
        vec3 diskLight = sampleAccretionDisk(currentPos, currentDir, startPos);
        if (length(diskLight) > 0.01) {
            result.accumulatedColor += diskLight * 0.1; // Reduced contribution for multiple sampling
        }
        
        // If we've moved far from the black hole, we can stop
        if (distToBlackHole > 20.0 * schwarzschildRadius) {
            break;
        }
        
        // Adaptive step sizing - smaller steps closer to the black hole
        stepSize = max(0.1 * schwarzschildRadius, 
                      min(1.0 * schwarzschildRadius, distToBlackHole * 0.1));
    }
    
    result.finalDirection = currentDir;
    return result;
}

// Generate visible light ray effects
vec3 generateLightRayVisualization(vec2 screenCoords, vec3 cameraPos, vec3 rayDir) {
    vec3 rayColor = vec3(0.0);
    
    // Create multiple light sources positioned around the black hole
    vec3 lightSources[6];
    lightSources[0] = blackholePos + vec3(15.0 * schwarzschildRadius, 0.0, 0.0);
    lightSources[1] = blackholePos + vec3(-15.0 * schwarzschildRadius, 0.0, 0.0);
    lightSources[2] = blackholePos + vec3(0.0, 15.0 * schwarzschildRadius, 0.0);
    lightSources[3] = blackholePos + vec3(0.0, -15.0 * schwarzschildRadius, 0.0);
    lightSources[4] = blackholePos + vec3(0.0, 0.0, 15.0 * schwarzschildRadius);
    lightSources[5] = blackholePos + vec3(0.0, 0.0, -15.0 * schwarzschildRadius);
    
    // Generate rays from multiple sources
    for (int source = 0; source < 6; source++) {
        vec3 sourcePos = lightSources[source];
        
        // Generate rays in a denser pattern
        for (float i = -3.0; i <= 3.0; i += 0.4) {
            for (float j = -3.0; j <= 3.0; j += 0.4) {
                vec3 testRayStart = sourcePos + vec3(i * schwarzschildRadius * 0.5, 
                                                   j * schwarzschildRadius * 0.5, 0.0);
                vec3 testRayDir = normalize(blackholePos - testRayStart);
                
                // Trace this test ray
                RayTraceResult traceResult = traceGeodesicRay(testRayStart, testRayDir, 15);
                
                // Check if this ray passes close to our viewing ray (more generous threshold)
                vec3 rayStart = testRayStart;
                vec3 rayEnd = testRayStart + traceResult.finalDirection * 30.0 * schwarzschildRadius;
                
                // Calculate closest approach between our view ray and this light ray
                vec3 viewRayStart = cameraPos;
                vec3 viewRayEnd = cameraPos + rayDir * 50.0 * schwarzschildRadius;
                
                // Line-to-line distance calculation
                vec3 w = rayStart - viewRayStart;
                vec3 u = normalize(rayEnd - rayStart);
                vec3 v = normalize(viewRayEnd - viewRayStart);
                
                float a = dot(u, u);
                float b = dot(u, v);
                float c = dot(v, v);
                float d = dot(u, w);
                float e = dot(v, w);
                
                float denominator = a * c - b * b;
                float distanceToRay = length(w + u * d / a - v * e / c);
                
                // Increased visibility threshold
                if (distanceToRay < 2.0 * schwarzschildRadius) {
                    // Enhanced ray visibility with better intensity falloff
                    float rayIntensity = exp(-distanceToRay / schwarzschildRadius);
                    
                    // Color the ray based on deflection and source position
                    vec3 spectralColor = vec3(1.2, 0.9, 0.6); // Brighter yellowish light
                    
                    if (enableSpectralRendering == 1) {
                        // Apply spectral shifting based on deflection
                        float redshiftFactor = 1.0 + traceResult.totalDeflection * 0.2;
                        float baseWavelength = 550.0; // Green light
                        float shiftedWavelength = baseWavelength * redshiftFactor;
                        spectralColor = wavelengthToRGB(shiftedWavelength);
                        
                        // Add some blue shift for rays coming from different directions
                        if (source % 2 == 1) {
                            shiftedWavelength = baseWavelength * 0.8; // Blue shift
                            spectralColor = mix(spectralColor, wavelengthToRGB(shiftedWavelength), 0.3);
                        }
                    }
                    
                    // Much higher intensity multiplier for visibility
                    rayColor += spectralColor * rayIntensity * 1.5;
                }
            }
        }
    }
    
    return rayColor;
}

// Relativistic jet emission calculation
vec3 calculateJetEmission(vec3 rayPos, vec3 rayDir, vec3 observerPos) {
    if (showJets == 0) return vec3(0.0);
    
    vec3 jetDirection = normalize(jetAxis);
    vec3 jetBase = blackholePos + jetDirection * (2.0 * schwarzschildRadius);
    float jetHeight = 20.0 * schwarzschildRadius;
    
    vec3 jetColor = vec3(0.0);
    
    // Sample along jet axis
    for (int i = 0; i < 20; i++) {
        float heightFraction = float(i) / 19.0;
        float height = 2.0 * schwarzschildRadius + heightFraction * jetHeight;
        vec3 jetPoint = blackholePos + jetDirection * height;
        
        // Check if ray passes near this jet element
        vec3 toJetPoint = jetPoint - rayPos;
        vec3 closestPoint = rayPos + rayDir * dot(toJetPoint, rayDir);
        float distanceToJet = length(closestPoint - jetPoint);
        
        // Jet radius at this height (conical expansion)
        float jetRadius = 0.1 * schwarzschildRadius * (1.0 + height * tan(jetOpeningAngle) / (10.0 * schwarzschildRadius));
        
        if (distanceToJet < jetRadius) {
            // Calculate relativistic beaming
            float beta = sqrt(1.0 - 1.0/(jetLorentzFactor * jetLorentzFactor));
            vec3 jetVelocity = jetDirection * beta;
            vec3 observerDirection = normalize(observerPos - jetPoint);
            
            // Doppler beaming factor
            float cosTheta = dot(jetVelocity, observerDirection);
            float beamingFactor = jetLorentzFactor * (1.0 + beta * cosTheta);
            
            // Check if within beaming cone
            float beamingConeAngle = 1.0 / jetLorentzFactor;
            vec3 jetToObserver = normalize(observerPos - jetPoint);
            float viewingAngle = acos(dot(jetDirection, jetToObserver));
            
            if (viewingAngle < 2.0 * beamingConeAngle) {
                // Synchrotron emission color - magnetic field dependent
                float magneticFieldStrength = jetMagneticField * pow(2.0 * schwarzschildRadius / height, 1.5);
                float syncFreq = magneticFieldStrength * 1e15; // Simplified
                
                vec3 spectralColor;
                if (enableSpectralRendering == 1) {
                    // Physics-based spectrum
                    float temperature = 1e9 * pow(magneticFieldStrength / jetMagneticField, 2.0/3.0);
                    spectralColor = temperatureToRGB(temperature);
                } else {
                    // Simplified blue-white jet color
                    spectralColor = vec3(0.8, 0.9, 1.2);
                }
                
                // Apply relativistic beaming enhancement
                float intensity = 1.0 / (1.0 + pow(distanceToJet / jetRadius, 2.0));
                
                if (viewingAngle < beamingConeAngle) {
                    // Strong beaming - intensity scales as δ³
                    spectralColor *= pow(beamingFactor, 3.0) * intensity;
                } else {
                    // Transition region with exponential suppression
                    float suppression = exp(-(viewingAngle - beamingConeAngle) / beamingConeAngle);
                    spectralColor *= beamingFactor * suppression * intensity;
                }
                
                // Add variability based on turbulence
                float turbulence = 0.5 + 0.5 * sin(time * 3.0 + height * 0.1) * cos(time * 2.0 + viewingAngle * 10.0);
                spectralColor *= turbulence;
                
                jetColor += spectralColor * 0.1; // Scale factor
            }
        }
    }
    
    return jetColor;
}

// Enhanced photon sphere with jet interactions  
vec3 calculatePhotonSphereWithJets(vec3 rayDir, vec3 observerPos) {
    vec3 toBlackHole = normalize(blackholePos - observerPos);
    float alignment = dot(rayDir, toBlackHole);
    
    if (alignment < 0.85) return vec3(0.0);
    
    float glow = pow((alignment - 0.85) / 0.15, 2.0);
    vec3 baseGlow = vec3(1.5, 0.8, 0.2);
    
    // Add jet-enhanced emission if jets are active
    if (showJets == 1) {
        vec3 jetDirection = normalize(jetAxis);
        float jetAlignment = dot(rayDir, jetDirection);
        
        if (jetAlignment > 0.9) {
            // Looking along jet axis - strong beaming enhancement
            float jetBeaming = pow(jetLorentzFactor, 2.0) * (1.0 + jetAlignment);
            vec3 jetGlow = vec3(0.6, 0.8, 1.2) * jetBeaming * 0.1;
            baseGlow += jetGlow;
        }
    }
    
    return baseGlow * glow;
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
    
    // Apply enhanced gravitational lensing to our view ray
    vec3 lensedRayDir = applyGravitationalLensing(cameraPos, rayDir);
    
    // Start with enhanced star field background using lensed ray direction
    vec3 color = generateStars(lensedRayDir, cameraPos);
    
    // Add light ray visualization if enabled
    if (showLightRays == 1) {
        vec3 rayVisualization = generateLightRayVisualization(screenCoords, cameraPos, rayDir);
        color += rayVisualization;
    }
    
    // Add relativistic jet emission
    vec3 jetEmission = calculateJetEmission(cameraPos, lensedRayDir, cameraPos);
    color += jetEmission;
    
    // Check if ray hits black hole (using original ray for intersection)
    if (discriminant >= 0.0) {
        float t = (-b - sqrt(discriminant)) / (2.0 * a);
        if (t > 0.0) {
            // Ray hits event horizon - pure black
            color = vec3(0.0, 0.0, 0.0);
        }
    } else {
        // Trace geodesic for enhanced accuracy
        RayTraceResult geodesicResult = traceGeodesicRay(cameraPos, rayDir, 15);
        
        if (geodesicResult.hitEventHorizon) {
            color = vec3(0.0, 0.0, 0.0);
        } else {
            // Sample accretion disk using both original and deflected rays
            vec3 diskContribution = sampleAccretionDisk(cameraPos, lensedRayDir, cameraPos);
            vec3 geodesicDiskContribution = geodesicResult.accumulatedColor;
            
            if (length(diskContribution) > 0.01 || length(geodesicDiskContribution) > 0.01) {
                // Combine both contributions with enhanced spectral effects
                vec3 totalDisk = diskContribution + geodesicDiskContribution * 0.5;
                
                // Apply chromatic dispersion - different wavelengths bend differently
                if (enableSpectralRendering == 1 && geodesicResult.totalDeflection > 0.01) {
                    // Simulate chromatic aberration from gravitational lensing
                    float redDeflection = geodesicResult.totalDeflection * 0.95;  // Red bends less
                    float greenDeflection = geodesicResult.totalDeflection;       // Green baseline
                    float blueDeflection = geodesicResult.totalDeflection * 1.05; // Blue bends more
                    
                    // Slightly separate the color channels
                    vec2 redOffset = screenCoords + vec2(redDeflection * 0.001, 0.0);
                    vec2 blueOffset = screenCoords - vec2(blueDeflection * 0.001, 0.0);
                    
                    totalDisk.r *= (1.0 + sin(redOffset.x * 50.0) * 0.1);
                    totalDisk.b *= (1.0 + sin(blueOffset.x * 50.0) * 0.1);
                }
                
                color = totalDisk;
            } else {
                // Enhanced photon sphere effects with jet interactions
                vec3 photonSphereGlow = calculatePhotonSphereWithJets(lensedRayDir, cameraPos);
                
                float photonSphereRadius = 1.5 * schwarzschildRadius;
                float distanceToCenter = length(cameraPos - blackholePos);
                
                if (length(photonSphereGlow) > 0.0 && distanceToCenter > photonSphereRadius) {
                    // Enhanced physics-based glow with Einstein ring effects
                    vec3 glowColor = photonSphereGlow;
                    
                    if (enableSpectralRendering == 1) {
                        float ambientTemp = 3000.0;
                        float gravitationalBlueshift = sqrt(1.0 - schwarzschildRadius / distanceToCenter);
                        float lensingIntensification = 1.0 + geodesicResult.totalDeflection;
                        
                        float shiftedTemp = ambientTemp * gravitationalBlueshift * lensingIntensification;
                        vec3 tempColor = temperatureToRGB(shiftedTemp);
                        glowColor = mix(glowColor, tempColor, 0.5); // Blend with temperature
                    }
                    
                    // Einstein ring pattern with jet modulation
                    float ringPattern = sin(distanceToCenter / schwarzschildRadius * 10.0 + time) * 0.1 + 0.9;
                    
                    // Apply enhanced photon sphere glow from uniform parameter
                    color += glowColor * ringPattern * (1.0 + geodesicResult.totalDeflection) * photonSphereGlow;
                }
                
                // Add gravitational lensing magnification to background stars
                if (geodesicResult.totalDeflection > 0.01) {
                    float magnification = 1.0 + geodesicResult.totalDeflection * 2.0;
                    color *= magnification;
                }
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
    
    // =========================================================================
    // DRAMATIC PHYSICS-BASED PROXIMITY EFFECTS (Standard Shader Version)
    // =========================================================================
    
    if (enableDramaticEffects == 1 && proximityFactor > 0.01) {
        vec2 screenCoords = (gl_FragCoord.xy / resolution) * 2.0 - 1.0;
        screenCoords.x *= resolution.x / resolution.y;
        
        // Stage 1: Early spacetime distortion
        if (proximityFactor > 0.01) {
            float intensity = proximityFactor * 1.5;
            vec2 p = screenCoords * (1.0 + proximityFactor * 2.0);
            
            float distortion = sin(p.x * 10.0 + time * 2.0) * cos(p.y * 8.0 + time * 1.5) * intensity;
            color.r += abs(distortion) * 0.8;
            color.g += abs(distortion) * 0.4;
            color.b += abs(distortion) * 0.2;
        }
        
        // Stage 2: Gravitational redshift
        if (proximityFactor > 0.1) {
            float redshift = 1.0 / sqrt(max(0.05, 1.0 - proximityFactor * 0.8));
            color.r *= (1.0 + (redshift - 1.0) * 1.2);
            color.g *= (1.0 + (redshift - 1.0) * 0.6);
            color.b /= (1.0 + (redshift - 1.0) * 1.5);
            
            // Pulsing effects
            if (proximityFactor > 0.3) {
                float pulse = sin(time * 3.0) * 0.5 + 0.5;
                color += vec3(pulse * proximityFactor * 0.6, 0.0, 0.0);
            }
        }
        
        // Stage 3: Intense warping
        if (proximityFactor > 0.4) {
            vec2 warp = screenCoords * proximityFactor * 2.0;
            float warpAngle = proximityFactor * time * 3.0;
            vec2 warped = mat2(cos(warpAngle), -sin(warpAngle), sin(warpAngle), cos(warpAngle)) * warp;
            
            float warpIntensity = sin(warped.x * 15.0) * cos(warped.y * 12.0) * proximityFactor * 1.0;
            color += vec3(abs(warpIntensity), abs(warpIntensity) * 0.6, abs(warpIntensity) * 0.3);
        }
        
        // Stage 4: Extreme effects near event horizon
        if (proximityFactor > 0.7) {
            float chaos = sin(screenCoords.x * 30.0 + time * 8.0) * cos(screenCoords.y * 25.0 + time * 6.0);
            float firewall = proximityFactor * proximityFactor * abs(chaos) * 2.0;
            color = mix(color, vec3(1.0, 0.6, 0.2) * firewall, smoothstep(0.7, 1.0, proximityFactor) * 0.5);
        }
        
        // Final intensity boost
        if (proximityFactor > 0.5) {
            color *= (1.0 + proximityFactor * 0.4);
        }
    }
    
    FragColor = vec4(color, 1.0);
}
