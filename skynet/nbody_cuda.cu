// nbody_cuda.cu - CUDA implementation file
#include "nbody_cuda.cuh"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <math.h>

// Constant memory for physical constants
__constant__ double d_G;

// CUDA kernel to calculate forces between bodies
__global__ void calculateForces(Body* bodies, int numBodies) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < numBodies) {
        // Reset acceleration
        bodies[idx].ax = bodies[idx].ay = bodies[idx].az = 0.0;
        
        // Calculate forces from all other bodies
        for (int j = 0; j < numBodies; j++) {
            if (idx != j) {
                double dx = bodies[j].x - bodies[idx].x;
                double dy = bodies[j].y - bodies[idx].y;
                double dz = bodies[j].z - bodies[idx].z;
                
                double distSq = dx*dx + dy*dy + dz*dz;
                double dist = sqrt(distSq + 1e-10); // Softening to avoid division by zero
                
                double F = d_G * bodies[idx].m * bodies[j].m / (distSq + 1e-10);
                
                // Force components
                double Fx = F * dx / dist;
                double Fy = F * dy / dist;
                double Fz = F * dz / dist;
                
                // Accumulate acceleration
                bodies[idx].ax += Fx / bodies[idx].m;
                bodies[idx].ay += Fy / bodies[idx].m;
                bodies[idx].az += Fz / bodies[idx].m;
            }
        }
    }
}

// CUDA kernel to update positions and velocities
__global__ void integratePositions(Body* bodies, int numBodies, double dt) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < numBodies) {
        // Update velocity
        bodies[idx].vx += bodies[idx].ax * dt;
        bodies[idx].vy += bodies[idx].ay * dt;
        bodies[idx].vz += bodies[idx].az * dt;
        
        // Update position (with semi-implicit Euler integration)
        bodies[idx].x += bodies[idx].vx * dt + 0.5 * bodies[idx].ax * dt * dt;
        bodies[idx].y += bodies[idx].vy * dt + 0.5 * bodies[idx].ay * dt * dt;
        bodies[idx].z += bodies[idx].vz * dt + 0.5 * bodies[idx].az * dt * dt;
    }
}

// Helper function to check for collisions on GPU and mark bodies for merging
__global__ void detectCollisions(Body* bodies, int* collisionMap, int numBodies) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < numBodies && collisionMap[idx] == -1) {
        for (int j = idx + 1; j < numBodies; j++) {
            if (collisionMap[j] == -1) { // Only check non-merged bodies
                double dx = bodies[idx].x - bodies[j].x;
                double dy = bodies[idx].y - bodies[j].y;
                double dz = bodies[idx].z - bodies[j].z;
                
                double distSq = dx*dx + dy*dy + dz*dz;
                double radiusSum = bodies[idx].radius + bodies[j].radius;
                
                if (distSq < radiusSum * radiusSum) {
                    // Mark the smaller body for merging into the larger one
                    if (bodies[idx].m < bodies[j].m) {
                        collisionMap[idx] = j;
                    } else {
                        collisionMap[j] = idx;
                    }
                }
            }
        }
    }
}

// Implementation of the CUDA step function
void cudaStep(Body* h_bodies, int& numBodies, double dt) {
    // Allocate device memory
    Body* d_bodies;
    int* d_collisionMap;
    cudaMalloc((void**)&d_bodies, numBodies * sizeof(Body));
    cudaMalloc((void**)&d_collisionMap, numBodies * sizeof(int));
    
    // Copy data to device
    cudaMemcpy(d_bodies, h_bodies, numBodies * sizeof(Body), cudaMemcpyHostToDevice);
    
    // Set physical constants in constant memory
    double h_G = 6.67430e-11;
    cudaMemcpyToSymbol(d_G, &h_G, sizeof(double));
    
    // Initialize collision map to -1 (no collisions)
    int* h_collisionMap = new int[numBodies];
    for (int i = 0; i < numBodies; i++) {
        h_collisionMap[i] = -1;
    }
    cudaMemcpy(d_collisionMap, h_collisionMap, numBodies * sizeof(int), cudaMemcpyHostToDevice);
    
    // Calculate number of CUDA blocks and threads
    int threadsPerBlock = 256;
    int blocksPerGrid = (numBodies + threadsPerBlock - 1) / threadsPerBlock;
    
    // Execute kernels
    calculateForces<<<blocksPerGrid, threadsPerBlock>>>(d_bodies, numBodies);
    integratePositions<<<blocksPerGrid, threadsPerBlock>>>(d_bodies, numBodies, dt);
    detectCollisions<<<blocksPerGrid, threadsPerBlock>>>(d_bodies, d_collisionMap, numBodies);
    
    // Synchronize to ensure all kernels have completed
    cudaDeviceSynchronize();
    
    // Copy results back to host
    cudaMemcpy(h_bodies, d_bodies, numBodies * sizeof(Body), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_collisionMap, d_collisionMap, numBodies * sizeof(int), cudaMemcpyDeviceToHost);
    
    // Handle collisions and mergers on CPU
    int newNumBodies = numBodies;
    for (int i = 0; i < numBodies; i++) {
        if (h_collisionMap[i] != -1) {
            int targetIdx = h_collisionMap[i];
            
            // Conservation of momentum
            h_bodies[targetIdx].vx = (h_bodies[targetIdx].m * h_bodies[targetIdx].vx + h_bodies[i].m * h_bodies[i].vx) / 
                                     (h_bodies[targetIdx].m + h_bodies[i].m);
            h_bodies[targetIdx].vy = (h_bodies[targetIdx].m * h_bodies[targetIdx].vy + h_bodies[i].m * h_bodies[i].vy) / 
                                     (h_bodies[targetIdx].m + h_bodies[i].m);
            h_bodies[targetIdx].vz = (h_bodies[targetIdx].m * h_bodies[targetIdx].vz + h_bodies[i].m * h_bodies[i].vz) / 
                                     (h_bodies[targetIdx].m + h_bodies[i].m);
            
            // Combine masses
            h_bodies[targetIdx].m += h_bodies[i].m;
            
            // Update radius (assuming volume conservation)
            h_bodies[targetIdx].radius = pow(pow(h_bodies[targetIdx].radius, 3) + pow(h_bodies[i].radius, 3), 1.0/3.0);
            
            // Mark this body for removal
            h_bodies[i].m = 0;
            newNumBodies--;
        }
    }
    
    // Compact the array to remove merged bodies
    if (newNumBodies < numBodies) {
        int writeIdx = 0;
        for (int i = 0; i < numBodies; i++) {
            if (h_bodies[i].m > 0) {
                if (writeIdx != i) {
                    h_bodies[writeIdx] = h_bodies[i];
                }
                writeIdx++;
            }
        }
        numBodies = newNumBodies;
    }
    
    // Clean up
    cudaFree(d_bodies);
    cudaFree(d_collisionMap);
    delete[] h_collisionMap;
}
