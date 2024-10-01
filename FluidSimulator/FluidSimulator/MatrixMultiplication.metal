//
//  MatrixMultiplication.metal
//  FluidSimulator
//
//  Created by Nick Venanzi on 9/29/24.
//

#include <metal_stdlib>
using namespace metal;

struct MatrixDimensions {
    uint n;
};

kernel void matrix_multiplication(const device float* A [[buffer(0)]],
                                  const device float* b [[buffer(1)]],
                                  device float* result [[buffer(2)]],
                                  constant MatrixDimensions &dims [[buffer(3)]],
                                  uint row [[thread_position_in_grid]]) {
    // Initialize the result to zero
    float sum = 0.0;

    // Perform the dot product of row from matrix A and column from matrix B
    for (uint i = 0; i < dims.n; i++) {
        sum += A[row * dims.n + i] * b[i];
    }

    // Write the result
    result[row] = sum;
}
