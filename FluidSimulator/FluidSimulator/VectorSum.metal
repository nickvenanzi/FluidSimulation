//
//  VectorSum.metal
//  FluidSimulator
//
//  Created by Nick Venanzi on 10/10/24.
//

#include <metal_stdlib>
using namespace metal;
constant int MAX_THREADGROUP_SIZE = 1024;

kernel void vector_sum(
    device const float* input [[buffer(0)]],
    device float* partialSums [[buffer(1)]],
    uint gid [[thread_position_in_grid]],
    uint tid [[thread_position_in_threadgroup]],
    uint tg_size [[threads_per_threadgroup]]
) {
    threadgroup float localData[MAX_THREADGROUP_SIZE];
    localData[tid] = input[gid];
    if (tid + tg_size < MAX_THREADGROUP_SIZE) {
        localData[tid+tg_size] = 0.0;
    }
    uint start_stride;
    if (tg_size <= 2) {
        start_stride = 1;
    } else if (tg_size <= 4) {
        start_stride = 2;
    } else if (tg_size <= 8) {
        start_stride = 4;
    } else if (tg_size <= 16) {
        start_stride = 8;
    } else if (tg_size <= 32) {
        start_stride = 16;
    } else if (tg_size <= 64) {
        start_stride = 32;
    } else if (tg_size <= 128) {
        start_stride = 64;
    } else if (tg_size <= 256) {
        start_stride = 128;
    } else if (tg_size <= 512) {
        start_stride = 256;
    } else {
        start_stride = 512;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint stride = start_stride; stride > 0; stride /= 2) {
        if (tid < stride) {
            localData[tid] += localData[tid + stride];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    if (tid == 0) {
        partialSums[gid / tg_size] = localData[0];
    }
}
