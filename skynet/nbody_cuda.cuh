// nbody_cuda.cuh - Header file
#ifndef NBODY_CUDA_H
#define NBODY_CUDA_H

typedef struct my_Color {
    unsigned char r, g, b,a;
} my_Color;

struct Body {
    double x, y, z;        // Position
    double vx, vy, vz;     // Velocity
    double ax, ay, az;     // Acceleration
    double m;              // Mass
    double radius;         // For collision detection
    my_Color color;           // Color
};

// Function declarations
void cudaStep(Body* bodies, int& numBodies, double dt);

#endif
