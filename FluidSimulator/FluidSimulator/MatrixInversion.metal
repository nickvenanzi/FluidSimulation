//
//  MatrixInversion.metal
//  FluidSimulator
//
//  Created by Nick Venanzi on 10/11/24.
//

#include <metal_stdlib>
using namespace metal;

kernel void cg_initialization(const device float* A [[buffer(0)]],
                              const device float* AplusI [[buffer(1)]],
                              const device float* AplusJ [[buffer(2)]],
                              const device float* AplusK [[buffer(3)]],
                              const device float* x [[buffer(4)]],
                              const device float* b [[buffer(5)]],
                              device float* r [[buffer(6)]],
                              device float* rr [[buffer(7)]],
                              device float* p [[buffer(8)]],
                              uint3 gid [[thread_position_in_grid]],
                              uint3 g_size [[threads_per_grid]]) {
    uint x_jmp = g_size.y * g_size.z; // eventually replace with buffer pointing to dimens
    uint y_jmp = g_size.z;
    uint z_jmp = 1;
    uint index = gid.x * x_jmp + gid.y * y_jmp + gid.z;
    
    if (
        gid.x == 0 || gid.x == g_size.x - 1 ||
        gid.y == 0 || gid.y == g_size.y - 1 ||
        gid.z == 0 || gid.z == g_size.z - 1
    ) {
        r[index] = 0;
        rr[index] = 0;
        p[index] = 0;
        return;
    }
    
    float tmp = A[index]*x[index] + AplusI[index]*x[index + x_jmp] + AplusJ[index]*x[index + y_jmp] + AplusK[index]*x[index + z_jmp] + AplusI[index - x_jmp]*x[index - x_jmp] + AplusJ[index - y_jmp]*x[index - y_jmp] + AplusK[index - z_jmp]*x[index - z_jmp];
    float result = b[index] - tmp;
    r[index] = result;
    rr[index] = result * result;
    p[index] = result;
}

kernel void cg_iterationStep1(const device float* A [[buffer(0)]],
                              const device float* AplusI [[buffer(1)]],
                              const device float* AplusJ [[buffer(2)]],
                              const device float* AplusK [[buffer(3)]],
                              const device float* p [[buffer(4)]],
                              device float* Ap [[buffer(5)]],
                              device float* pAp [[buffer(6)]],
                              uint3 gid [[thread_position_in_grid]],
                              uint3 g_size [[threads_per_grid]]) {
    uint x_jmp = g_size.y * g_size.z;
    uint y_jmp = g_size.z;
    uint z_jmp = 1;
    uint index = gid.x * x_jmp + gid.y * y_jmp + gid.z;
    
    if (
        gid.x == 0 || gid.x == g_size.x - 1 ||
        gid.y == 0 || gid.y == g_size.y - 1 ||
        gid.z == 0 || gid.z == g_size.z - 1
    ) {
        Ap[index] = 0;
        pAp[index] = 0;
        return;
    }
    
    float tmp = A[index]*p[index] + AplusI[index]*p[index + x_jmp] + AplusJ[index]*p[index + y_jmp] + AplusK[index]*p[index + z_jmp] + AplusI[index - x_jmp]*p[index - x_jmp] + AplusJ[index - y_jmp]*p[index - y_jmp] + AplusK[index - z_jmp]*p[index - z_jmp];
    Ap[index] = tmp;
    pAp[index] = p[index] * tmp;
}

kernel void cg_iterationStep3(const device float* p [[buffer(0)]],
                              const device float* Ap [[buffer(1)]],
                              const device float* pAp [[buffer(2)]],
                              const device float* rrOld [[buffer(3)]],
                              device float* x [[buffer(4)]],
                              device float* r [[buffer(5)]],
                              device float* rr [[buffer(6)]],
                              uint index [[thread_position_in_grid]]) {
    float alpha = rrOld[0] / pAp[0];
    x[index] = x[index] + alpha * p[index];
    float tmp = r[index] - alpha * Ap[index];
    r[index] = tmp;
    rr[index] = tmp*tmp;
}

kernel void cg_iterationStep5(const device float* r [[buffer(0)]],
                              const device float* rrOld [[buffer(1)]],
                              const device float* rr [[buffer(2)]],
                              device float* p [[buffer(3)]],
                              uint index [[thread_position_in_grid]]) {
    p[index] = r[index] + (rr[0] / rrOld[0]) * p[index];
}

kernel void cg_iterationStep6(const device float* rr [[buffer(0)]],
                              device float* rrOld [[buffer(1)]]) {
    rrOld[0] = rr[0];
}
