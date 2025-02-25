#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <cmath>
#include <iostream>
#include <memory>
#include "raylib.h"

struct Body {
    double m;
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    Color color;
};

class Simulation {
public:
    Simulation(bool useSolarSystem = true) {
        if (useSolarSystem) {
            initSolarSystem();
        } else {
            initTripleStarSystem();
        }
    }

    void reset(bool useSolarSystem) {
        std::lock_guard<std::mutex> lock(mutex_);
        bodies_.clear();
        
        if (useSolarSystem) {
            initSolarSystem();
        } else {
            initTripleStarSystem();
        }
    }

    void initSolarSystem() {
        const double G = 6.67430e-11;  // Gravitational constant
        const double AU = 1.496e11;    // Astronomical unit in meters
        
        // Create Sun
        Body sun;
        sun.m = 1.989e30;
        sun.x = 0;
        sun.y = 0;
        sun.z = 0;
        sun.vx = 0;
        sun.vy = 0;
        sun.vz = 0;
        sun.color = YELLOW;
        bodies_.push_back(sun);

        // Helper function to create planets
        auto createPlanet = [&](double mass, double au_distance, double orbital_speed, Color color) {
            Body planet;
            planet.m = mass;
            planet.x = au_distance * AU;  // Place along x-axis
            planet.y = 0;
            planet.z = 0;
            planet.vx = 0;
            planet.vy = orbital_speed;    // Initial tangential velocity
            planet.vz = 0;
            planet.color = color;
            return planet;
        };

        // Add inner planets (approximate values)
        bodies_.push_back(createPlanet(
            3.301e23,   // Mercury mass
            0.387,      // AU distance
            std::sqrt(G * sun.m / (0.387 * AU)),  // Orbital velocity
            GRAY
        ));

        bodies_.push_back(createPlanet(
            4.867e24,   // Venus mass
            0.723,      // AU distance
            std::sqrt(G * sun.m / (0.723 * AU)),  // ~35,000 m/s
            BEIGE
        ));

        // Add Earth
        bodies_.push_back(createPlanet(
            5.972e24,   // Earth mass
            1.0,        // AU distance
            std::sqrt(G * sun.m / (1.0 * AU)),    // ~29,780 m/s
            BLUE
        ));

        bodies_.push_back(createPlanet(
            6.417e23,   // Mars mass
            1.524,      // AU distance
            std::sqrt(G * sun.m / (1.524 * AU)),  // ~24,100 m/s
            RED
        ));

        // Add outer gas giants
        bodies_.push_back(createPlanet(
            1.898e27,   // Jupiter mass
            5.203,      // AU distance
            std::sqrt(G * sun.m / (5.203 * AU)),  // ~13,060 m/s
            ORANGE
        ));

        bodies_.push_back(createPlanet(
            5.683e26,   // Saturn mass
            9.537,      // AU distance
            std::sqrt(G * sun.m / (9.537 * AU)),  // ~9,680 m/s
            GOLD
        ));
    }

    void initTripleStarSystem() {
        const double G = 6.67430e-11;  // Gravitational constant
        const double AU = 1.496e11;    // 1 astronomical unit in meters
        const double mass_base = 1e28; // Base mass (similar to small stars)
        
        // Common velocity magnitude calculation
        auto calc_velocity = [&](double m1, double m2, double distance) {
            return std::sqrt(G * (m1 + m2) / distance);
        };
    
        // Triple star system with chaotic initial conditions
        Body star1, star2, star3;
        
        // Mass configuration (similar but not identical)
        star1.m = mass_base * 0.5;
        star2.m = mass_base * 0.6;
        star3.m = mass_base * 0.7;
    
        // Initial positions (triangular configuration)
        star1.x =  AU * 0.7;
        star1.y = -AU * 0.3;
        star2.x = -AU * 0.4;
        star2.y =  AU * 0.6;
        star3.x =  AU * 0.2;
        star3.y =  AU * 0.1;
    
        // Velocity configuration (non-symmetric, non-orbital)
        const double speed_base = 1e3; // 10 km/s base speed
        star1.vx =  speed_base * 0.01;
        star1.vy = -speed_base * 0.01;
        star2.vx = -speed_base * 1.9;
        star2.vy =  speed_base * 1.9;
        star3.vx =  speed_base * 1.5;
        star3.vy =  speed_base * 1.5;
    
        // Center-of-mass correction
        const double total_mass = star1.m + star2.m + star3.m;
        const double com_vx = (star1.m*star1.vx + star2.m*star2.vx + star3.m*star3.vx) / total_mass;
        const double com_vy = (star1.m*star1.vy + star2.m*star2.vy + star3.m*star3.vy) / total_mass;
        
        // Adjust velocities to maintain system in view
        star1.vx -= com_vx;
        star1.vy -= com_vy;
        star2.vx -= com_vx;
        star2.vy -= com_vy;
        star3.vx -= com_vx;
        star3.vy -= com_vy;
    
        bodies_.push_back(star1);
        bodies_.push_back(star2);
        bodies_.push_back(star3);
    }

    void addRandomBody(bool randomOrbit = false) {
        std::lock_guard<std::mutex> lock(mutex_);
        
        // Generate positions between -5.0 to +4.9 AU
        double x_au = (rand() % 100) / 10.0 - 5.0; 
        double y_au = (rand() % 100) / 10.0 - 5.0;

        // Convert AU to meters
        double x = x_au * 1.496e11;
        double y = y_au * 1.496e11;

        // Calculate distance from the Sun (assuming Sun is at origin)
        double r = std::sqrt(x*x + y*y);

        // Gravitational parameters (SI units)
        const double G = 6.67430e-11;      // Gravitational constant
        double M_sun = 0;
        
        // Get mass of central body (first in the array)
        if (!bodies_.empty()) {
            M_sun = bodies_[0].m;
        } else {
            M_sun = 1.989e30;  // Sun's mass as default
        }

        // Calculate required orbital velocity for circular orbit
        double v_mag = std::sqrt(G * M_sun / r);

        // Optional: Add slight randomness to velocity magnitude (0.9-1.1x)
        v_mag *= 0.9 + (rand() % 2001)/10000.0;  // Random factor between 0.9-1.1

        double vx, vy;
        
        if (randomOrbit) {
            // For random chaotic orbits
            vx = 10.0 * (-10 + (rand() % 20));
            vy = 10.0 * (-10 + (rand() % 20));
        } else {
            // For mostly circular/elliptical orbits
            // Set velocity direction (tangential, counter-clockwise)
            vx = -v_mag * y / r;
            vy = v_mag * x / r;
        }

        // Create and add body
        double earth_mass = 0.5 * 1e24;
        Body body;
        body.m = earth_mass;
        body.x = x;
        body.y = y;
        body.z = 0.0;
        body.vx = vx;
        body.vy = vy;
        body.vz = 0.0;
        body.ax = body.ay = body.az = 0.0;
        
        // Random color
        body.color = Color{
            255,
            static_cast<unsigned char>(rand() % 200 + 55),
            static_cast<unsigned char>(rand() % 200 + 55),
            static_cast<unsigned char>(rand() % 200 + 55)            
        };
        
        bodies_.push_back(body);
    }

    void removeMostDistant() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (bodies_.size() > 1) {
            double max_dist = 0;
            size_t max_id = 0;
            for (size_t i = 1; i < bodies_.size(); ++i) {
                const auto& body = bodies_[i];
                double dist = std::sqrt(body.x*body.x + body.y*body.y + body.z*body.z);
                if (dist > max_dist) {
                    max_dist = dist;
                    max_id = i;
                }
            }
            bodies_.erase(bodies_.begin() + max_id);
        }
    }

    void step(double dt) {
        std::lock_guard<std::mutex> lock(mutex_);
        const double G = 6.67430e-11;
        
        // Reset accelerations
        for (auto& body : bodies_) {
            body.ax = body.ay = body.az = 0.0;
        }

        // Calculate forces
        for (size_t i = 0; i < bodies_.size(); ++i) {
            for (size_t j = i + 1; j < bodies_.size(); ++j) {
                Body& a = bodies_[i];
                Body& b = bodies_[j];
                double dx = b.x - a.x;
                double dy = b.y - a.y;
                double dz = b.z - a.z;
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                // Avoid division by zero with softening
                double F = G * a.m * b.m / (dist*dist + 1e-10);  
                
                // Force components (direction from a to b)
                double Fx = F * dx / dist;
                double Fy = F * dy / dist;
                double Fz = F * dz / dist;

                // Apply equal and opposite forces
                a.ax += Fx / a.m;  // a is pulled toward b
                a.ay += Fy / a.m;
                a.az += Fz / a.m;
                
                b.ax -= Fx / b.m;  // b is pulled toward a
                b.ay -= Fy / b.m;
                b.az -= Fz / b.m;
            }
        }

        // Update velocities and positions
        for (auto& body : bodies_) {
            body.vx += body.ax * dt;
            body.vy += body.ay * dt;
            body.vz += body.az * dt;
            body.x += body.vx * dt + 0.5 * body.ax * dt * dt;
            body.y += body.vy * dt + 0.5 * body.ay * dt * dt;
            body.z += body.vz * dt + 0.5 * body.az * dt * dt;
        }
        
        // Check for collisions and handle mergers
        handleCollisions();
    }
    
    // NEW FUNCTION: Handle collisions between bodies
    void handleCollisions() {
        const double AU = 1.496e11;    // Astronomical unit in meters
        // Collision threshold: fraction of an AU (adjustable)
        const double collisionThreshold = 0.05 * AU; 
        
        // Use a flag to track if any collisions occurred
        bool collisionOccurred = true;
        
        // Continue checking for collisions until no more are found in a pass
        while (collisionOccurred && bodies_.size() > 1) {
            collisionOccurred = false;
            
            for (size_t i = 0; i < bodies_.size(); ++i) {
                for (size_t j = i + 1; j < bodies_.size(); ++j) {
                    Body& a = bodies_[i];
                    Body& b = bodies_[j];
                    
                    // Calculate distance between bodies
                    double dx = b.x - a.x;
                    double dy = b.y - a.y;
                    double dz = b.z - a.z;
                    double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                    
                    // Check if bodies are close enough to collide
                    if (dist < collisionThreshold) {
                        // Merge the bodies - conserve momentum (p = mv)
                        double totalMass = a.m + b.m;
                        
                        // Calculate new velocities based on conservation of momentum
                        double newVx = (a.m * a.vx + b.m * b.vx) / totalMass;
                        double newVy = (a.m * a.vy + b.m * b.vy) / totalMass;
                        double newVz = (a.m * a.vz + b.m * b.vz) / totalMass;
                        
                        // Calculate new position (mass-weighted center)
                        double newX = (a.m * a.x + b.m * b.x) / totalMass;
                        double newY = (a.m * a.y + b.m * b.y) / totalMass;
                        double newZ = (a.m * a.z + b.m * b.z) / totalMass;
                        
                        // Update the first body with new properties
                        a.m = totalMass;
                        a.x = newX;
                        a.y = newY;
                        a.z = newZ;
                        a.vx = newVx;
                        a.vy = newVy;
                        a.vz = newVz;
                        
                        // Mix colors (simple average)
                        a.color.r = (a.color.r + b.color.r) / 2;
                        a.color.g = (a.color.g + b.color.g) / 2;
                        a.color.b = (a.color.b + b.color.b) / 2;
                        
                        // Remove the second body
                        bodies_.erase(bodies_.begin() + j);
                        
                        // Flag that a collision occurred
                        collisionOccurred = true;
                        
                        // Break out of inner loop since bodies_ size changed
                        break;
                    }
                }
                
                // If a collision occurred, break out of outer loop too
                if (collisionOccurred) {
                    break;
                }
            }
        }
    }
    
    std::vector<Body> getBodies() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return bodies_;
    }

private:
    std::vector<Body> bodies_;
    mutable std::mutex mutex_;
};

// Global variables
bool pause = false;
float zoom = 1.0f;
float offsetX = 0.0f;
float offsetY = 0.0f;
bool showTrails = true;
bool showInfo = true;
int currentSystem = 0; // 0 = Solar System, 1 = Triple Star

// For tracking body positions for trails
const int MAX_TRAIL_POINTS = 300;
struct TrailPoint {
    Vector2 position;
    float alpha;
};
std::vector<std::vector<TrailPoint>> trails;

int main() {
    // Initialize window
    const int screenWidth = 1200;
    const int screenHeight = 800;
    InitWindow(screenWidth, screenHeight, "N-Body Gravity Simulation with Collisions");
    SetTargetFPS(60);

    // Create simulation
    Simulation simulation(true); // Start with solar system
    
    // Initialize random seed
    srand(time(NULL));
    
    // Time step for simulation
    const double timeStep = 3600.0 * 24; // One day in seconds
    const int stepsPerFrame = 5;         // Simulation steps per frame
    
    // Scale factor for visualization
    // 1 AU = 100 pixels
    const double AU = 1.496e11;    // 1 astronomical unit in meters
    const double visualScale = 100.0 / AU;
    
    // Initialize trails vector
    auto bodies = simulation.getBodies();
    trails.resize(bodies.size());
    
    // Collision counter
    int collisionCount = 0;
    int previousBodyCount = bodies.size();
    
    // Main game loop
    while (!WindowShouldClose()) {
        // Handle input
        if (IsKeyPressed(KEY_SPACE)) pause = !pause;
        if (IsKeyPressed(KEY_T)) showTrails = !showTrails;
        if (IsKeyPressed(KEY_I)) showInfo = !showInfo;
        
        // Zoom controls
        if (IsKeyDown(KEY_UP)) zoom *= 1.02f;
        if (IsKeyDown(KEY_DOWN)) zoom /= 1.02f;
        
        // Pan controls
        if (IsKeyDown(KEY_RIGHT)) offsetX -= 10.0f / zoom;
        if (IsKeyDown(KEY_LEFT)) offsetX += 10.0f / zoom;
        if (IsKeyDown(KEY_S)) offsetY -= 10.0f / zoom;
        if (IsKeyDown(KEY_W)) offsetY += 10.0f / zoom;
        
        // Add/remove bodies
        if (IsKeyPressed(KEY_A)) simulation.addRandomBody(false);
        if (IsKeyPressed(KEY_R)) simulation.addRandomBody(true); // Random chaotic orbit
        if (IsKeyPressed(KEY_D)) simulation.removeMostDistant();
        
        // Switch between solar system and triple star system
        if (IsKeyPressed(KEY_TAB)) {
            currentSystem = 1 - currentSystem; // Toggle between 0 and 1
            
            // Reset simulation with new system
            simulation.reset(currentSystem == 0);
            
            // Reset trails
            trails.clear();
            bodies = simulation.getBodies();
            trails.resize(bodies.size());
            
            // Reset view
            zoom = 1.0f;
            offsetX = 0.0f;
            offsetY = 0.0f;
            
            // Reset collision counter
            collisionCount = 0;
            previousBodyCount = bodies.size();
        }
        
        // Update simulation if not paused
        if (!pause) {
            for (int i = 0; i < stepsPerFrame; i++) {
                simulation.step(timeStep);
            }
        }
        
        // Get updated body positions
        bodies = simulation.getBodies();
        
        // Check if a collision occurred by comparing body count
        if (bodies.size() < previousBodyCount) {
            collisionCount += (previousBodyCount - bodies.size());
            previousBodyCount = bodies.size();
        }
        
        // Ensure trails vector size matches bodies
        if (trails.size() != bodies.size()) {
            // When bodies are merged, trails need to be adjusted
            if (trails.size() > bodies.size()) {
                // Bodies were removed - truncate trails
                trails.resize(bodies.size());
            } else {
                // New bodies were added
                trails.resize(bodies.size());
            }
        }
        
        // Update trails
        for (size_t i = 0; i < bodies.size(); i++) {
            Vector2 pos = {
                static_cast<float>(bodies[i].x * visualScale),
                static_cast<float>(bodies[i].y * visualScale)
            };
            
            if (showTrails) {
                // Add new point to trail
                if (trails[i].size() >= MAX_TRAIL_POINTS) {
                    trails[i].erase(trails[i].begin());
                }
                
                TrailPoint newPoint;
                newPoint.position = pos;
                newPoint.alpha = 1.0f;
                trails[i].push_back(newPoint);
                
                // Age all points in the trail
                for (auto& point : trails[i]) {
                    point.alpha -= 0.01f;
                    if (point.alpha < 0.0f) point.alpha = 0.0f;
                }
            } else {
                // Clear trails if not showing
                trails[i].clear();
            }
        }
        
        // Begin drawing
        BeginDrawing();
        ClearBackground(BLACK);
        
        // Calculate center of screen
        float centerX = screenWidth / 2.0f;
        float centerY = screenHeight / 2.0f;
        
        // Draw trails first (so they're behind the bodies)
        if (showTrails) {
            for (size_t i = 0; i < bodies.size() && i < trails.size(); i++) {
                Color trailColor = bodies[i].color;
                
                for (size_t j = 0; j < trails[i].size(); j++) {
                    // Fade out trail points based on age
                    Color pointColor = trailColor;
                    pointColor.a = static_cast<unsigned char>(255 * trails[i][j].alpha);
                    
                    DrawCircle(
                        centerX + (trails[i][j].position.x + offsetX) * zoom, 
                        centerY + (trails[i][j].position.y + offsetY) * zoom, 
                        1.0f, 
                        pointColor
                    );
                }
            }
        }
        
        // Draw bodies
        for (size_t i = 0; i < bodies.size(); i++) {
            // Calculate radius based on mass
            float baseSize = 2.0f;
            float radiusFactor = 0.0f;
            
            if (i == 0) {
                // Sun or central body
                radiusFactor = log10(bodies[i].m) - 20; // Adjusted for visibility
            } else {
                // Planets or smaller bodies
                radiusFactor = log10(bodies[i].m) - 22; // Adjusted for visibility
            }
            
            float radius = baseSize + radiusFactor * 2.0f;
            if (radius < 2.0f) radius = 2.0f;
            
            // Draw the body
            DrawCircle(
                centerX + (bodies[i].x * visualScale + offsetX) * zoom,
                centerY + (bodies[i].y * visualScale + offsetY) * zoom,
                radius * zoom,
                bodies[i].color
            );
        }
        
        // Draw UI information
        if (showInfo) {
            DrawText(TextFormat("Bodies: %d", (int)bodies.size()), 10, 10, 20, WHITE);
            DrawText(TextFormat("Zoom: %.2f", zoom), 10, 35, 20, WHITE);
            DrawText(TextFormat("System: %s", currentSystem == 0 ? "Solar System" : "Triple Star"), 10, 60, 20, WHITE);
            DrawText(TextFormat("Collisions: %d", collisionCount), 10, 85, 20, WHITE);
            
            DrawText("Controls:", 10, 120, 20, WHITE);
            DrawText("SPACE: Pause/Resume", 10, 145, 20, WHITE);
            DrawText("A: Add Stable Body", 10, 170, 20, WHITE);
            DrawText("R: Add Random Body", 10, 195, 20, WHITE);
            DrawText("D: Remove Distant Body", 10, 220, 20, WHITE);
            DrawText("T: Toggle Trails", 10, 245, 20, WHITE);
            DrawText("I: Toggle Info", 10, 270, 20, WHITE);
            DrawText("TAB: Switch System", 10, 295, 20, WHITE);
            DrawText("Arrows/WASD: Pan", 10, 320, 20, WHITE);
            DrawText("UP/DOWN: Zoom", 10, 345, 20, WHITE);
        }
        
        EndDrawing();
    }
    
    CloseWindow();
    
    return 0;
}