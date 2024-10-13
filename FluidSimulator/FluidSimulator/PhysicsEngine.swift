//
//  MatrixMultiplyView.swift
//  FluidSimulator
//
//  Created by Nick Venanzi on 9/29/24.
//

import SwiftUI
import Metal

class PhysicsEngine {
    
    let n: Int = 100
    
    var A: [Float]
    var AplusI: [Float]
    var AplusJ: [Float]
    var AplusK: [Float]
    
    var x: [Float]
    var b: [Float]
    
    var device: MTLDevice!
    var commandQueue: MTLCommandQueue!
    var cgInitPipeline: MTLComputePipelineState!
    var sumPipeline: MTLComputePipelineState!
    var cgIter1Pipeline: MTLComputePipelineState!
    var cgIter3Pipeline: MTLComputePipelineState!
    var cgIter5Pipeline: MTLComputePipelineState!
    var cgIter6Pipeline: MTLComputePipelineState!
    
    var aBuffer: MTLBuffer!
    var aPlusIBuffer: MTLBuffer!
    var aPlusJBuffer: MTLBuffer!
    var aPlusKBuffer: MTLBuffer!
    var xBuffer: MTLBuffer!
    var bBuffer: MTLBuffer!
    var rBuffer: MTLBuffer!
    var rrBuffer: MTLBuffer!
    var rrOldBuffer: MTLBuffer!
    var pBuffer: MTLBuffer!
    var ApBuffer: MTLBuffer!
    var pApBuffer: MTLBuffer!
    var partialSumsBuffer: MTLBuffer!
    var partialSums2Buffer: MTLBuffer!
    
    init() {
        A = Array(repeating: 0.0, count: n*n*n)
        AplusI = Array(repeating: 0.0, count: n*n*n)
        AplusJ = Array(repeating: 0.0, count: n*n*n)
        AplusK = Array(repeating: 0.0, count: n*n*n)
        x = Array(repeating: 0.0, count: n*n*n)
        b = Array(repeating: 0.0, count: n*n*n)
        generateData()
        
        // Get the Metal device and command queue
        device = MTLCreateSystemDefaultDevice()
        commandQueue = device.makeCommandQueue()
        
        guard let library = device.makeDefaultLibrary(), let cgInitializationFunc = library.makeFunction(name: "cg_initialization"), let sumFunc = library.makeFunction(name: "vector_sum"), let cgIter1 = library.makeFunction(name: "cg_iterationStep1"), let cgIter3 = library.makeFunction(name: "cg_iterationStep3"), let cgIter5 = library.makeFunction(name: "cg_iterationStep5"), let cgIter6 = library.makeFunction(name: "cg_iterationStep6"),
              let pipe1 = try? device.makeComputePipelineState(function: cgInitializationFunc),
              let pipe2 = try? device.makeComputePipelineState(function: sumFunc),
              let pipe3 = try? device.makeComputePipelineState(function: cgIter1),
              let pipe4 = try? device.makeComputePipelineState(function: cgIter3),
              let pipe5 = try? device.makeComputePipelineState(function: cgIter5),
              let pipe6 = try? device.makeComputePipelineState(function: cgIter6)
        else {
            print("Failed to set up Metal")
            return
        }
        cgInitPipeline = pipe1
        sumPipeline = pipe2
        cgIter1Pipeline = pipe3
        cgIter3Pipeline = pipe4
        cgIter5Pipeline = pipe5
        cgIter6Pipeline = pipe6
        
        computeX()
        
        /*
         Get time in seconds: let startTime = CFAbsoluteTimeGetCurrent()
         */
    }
    
    func printRR() {
        let ptr = rBuffer.contents().bindMemory(to: Float.self, capacity: n*n*n)
        var max: Float = 0
        for i in 0..<n*n*n {
            if abs(ptr[i]) > max {
                max = abs(ptr[i])
            }
        }
        print("Max element: \(max)")
        let b_ptr = rrOldBuffer.contents().bindMemory(to: Float.self, capacity: 1)
        print("RR: \(b_ptr[0])\n--------------")
    }
    
    func iterate() {
        let commandBuffer = commandQueue.makeCommandBuffer()!
        configureCGIteration(commandBuffer)
        commandBuffer.commit()
        commandBuffer.waitUntilCompleted()
        printRR()
    }
    
    func generateData() {
        var index: Int
        for i in 1..<n-1 {
            for j in 1..<n-1 {
                for k in 1..<n-1 {
                    index = i*n*n + j*n + k
                    
                    A[index] = Float.random(in: 5..<6)
                    AplusI[index] = Float.random(in: 0..<1)
                    AplusJ[index] = Float.random(in: 0..<1)
                    AplusK[index] = Float.random(in: 0..<1)
                    
                    x[index] = Float.random(in: 0..<1)
                }
            }
        }
        
        var b_val: Float = 0
        for i in 1..<n-1 {
            for j in 1..<n-1 {
                for k in 1..<n-1 {
                    index = idx(i, j, k)
                    b_val = A[index] * x[index]
                    b_val += AplusI[index] * x[idx(i+1, j, k)]
                    b_val += AplusJ[index] * x[idx(i, j+1, k)]
                    b_val += AplusK[index] * x[idx(i, j, k+1)]
                    b_val += AplusI[idx(i-1, j, k)] * x[idx(i-1, j, k)]
                    b_val += AplusJ[idx(i, j-1, k)] * x[idx(i, j-1, k)]
                    b_val += AplusK[idx(i, j, k-1)] * x[idx(i, j, k-1)]
                    
                    b[index] = b_val
                }
            }
        }
    }
    
    func idx(_ i: Int, _ j: Int, _ k: Int) -> Int {
        return i * n * n + j * n + k
    }
    
    func computeX() {
        // set buffers
        setBufferData()
        let startTime = CFAbsoluteTimeGetCurrent()

        // Create a command buffer
        let commandBuffer = commandQueue.makeCommandBuffer()!
        
        // CG Initialization
        configureCGInit(commandBuffer)
        
        // CG Iterations
        for _ in 0..<50 {
            configureCGIteration(commandBuffer)
        }
        
        // Commit the command buffer and wait for completion
        commandBuffer.commit()
        commandBuffer.waitUntilCompleted()
        
        let endTime = CFAbsoluteTimeGetCurrent()
        print("Time: \(endTime - startTime) seconds")
        
        printRR()
        return
    }
    
    func configureCGInit(_ commandBuffer: MTLCommandBuffer) {
        //initialization pipeline
        let computeEncoder = commandBuffer.makeComputeCommandEncoder()!
        computeEncoder.setComputePipelineState(cgInitPipeline)
        computeEncoder.setBuffer(aBuffer, offset: 0, index: 0)
        computeEncoder.setBuffer(aPlusIBuffer, offset: 0, index: 1)
        computeEncoder.setBuffer(aPlusJBuffer, offset: 0, index: 2)
        computeEncoder.setBuffer(aPlusKBuffer, offset: 0, index: 3)
        computeEncoder.setBuffer(xBuffer, offset: 0, index: 4)
        computeEncoder.setBuffer(bBuffer, offset: 0, index: 5)
        computeEncoder.setBuffer(rBuffer, offset: 0, index: 6)
        computeEncoder.setBuffer(rrBuffer, offset: 0, index: 7)
        computeEncoder.setBuffer(pBuffer, offset: 0, index: 8)
        
        let gridSize = MTLSize(width: n, height: n, depth: n)
        let threadgroupSize = MTLSize(width: 10, height: 10, depth: 10)
        computeEncoder.dispatchThreads(gridSize, threadsPerThreadgroup: threadgroupSize)
        computeEncoder.endEncoding()
        
        // initial RR compute
        configureSum(commandBuffer, rrBuffer, rrOldBuffer)
    }
    
    func configureCGIteration(_ commandBuffer: MTLCommandBuffer) {
        // Step 1: Compute A*p and p.A*p
        let computeEncoder = commandBuffer.makeComputeCommandEncoder()!
        computeEncoder.setComputePipelineState(cgIter1Pipeline)
        computeEncoder.setBuffer(aBuffer, offset: 0, index: 0)
        computeEncoder.setBuffer(aPlusIBuffer, offset: 0, index: 1)
        computeEncoder.setBuffer(aPlusJBuffer, offset: 0, index: 2)
        computeEncoder.setBuffer(aPlusKBuffer, offset: 0, index: 3)
        computeEncoder.setBuffer(pBuffer, offset: 0, index: 4)
        computeEncoder.setBuffer(ApBuffer, offset: 0, index: 5)
        computeEncoder.setBuffer(pApBuffer, offset: 0, index: 6)
        
        let gridSize = MTLSize(width: n, height: n, depth: n)
        let threadgroupSize = MTLSize(width: 10, height: 10, depth: 10)
        computeEncoder.dispatchThreads(gridSize, threadsPerThreadgroup: threadgroupSize)
        computeEncoder.endEncoding()
        
        //Step 2: Sum pAp (to implcitly get alpha)
        configureSum(commandBuffer, pApBuffer, pApBuffer)
        
        //Step 3: Compute x, r, and rr
        let computeEncoder3 = commandBuffer.makeComputeCommandEncoder()!
        computeEncoder3.setComputePipelineState(cgIter3Pipeline)
        computeEncoder3.setBuffer(pBuffer, offset: 0, index: 0)
        computeEncoder3.setBuffer(ApBuffer, offset: 0, index: 1)
        computeEncoder3.setBuffer(pApBuffer, offset: 0, index: 2)
        computeEncoder3.setBuffer(rrOldBuffer, offset: 0, index: 3)
        computeEncoder3.setBuffer(xBuffer, offset: 0, index: 4)
        computeEncoder3.setBuffer(rBuffer, offset: 0, index: 5)
        computeEncoder3.setBuffer(rrBuffer, offset: 0, index: 6)
        
        let gridSize3 = MTLSize(width: n*n*n, height: 1, depth: 1)
        let threadgroupSize3 = MTLSize(width: 1024, height: 1, depth: 1)
        computeEncoder3.dispatchThreads(gridSize3, threadsPerThreadgroup: threadgroupSize3)
        computeEncoder3.endEncoding()
        
        //Step 4: Compute r_new * r_new
        configureSum(commandBuffer, rrBuffer, rrBuffer)
        
        //Step 5: Compute p
        let computeEncoder5 = commandBuffer.makeComputeCommandEncoder()!
        computeEncoder5.setComputePipelineState(cgIter5Pipeline)
        computeEncoder5.setBuffer(rBuffer, offset: 0, index: 0)
        computeEncoder5.setBuffer(rrOldBuffer, offset: 0, index: 1)
        computeEncoder5.setBuffer(rrBuffer, offset: 0, index: 2)
        computeEncoder5.setBuffer(pBuffer, offset: 0, index: 3)

        let gridSize5 = MTLSize(width: n*n*n, height: 1, depth: 1)
        let threadgroupSize5 = MTLSize(width: 1024, height: 1, depth: 1)
        computeEncoder5.dispatchThreads(gridSize5, threadsPerThreadgroup: threadgroupSize5)
        computeEncoder5.endEncoding()
        
        //Step 6: Set rr_new to rr_old
        let computeEncoder6 = commandBuffer.makeComputeCommandEncoder()!
        computeEncoder6.setComputePipelineState(cgIter6Pipeline)
        computeEncoder6.setBuffer(rrBuffer, offset: 0, index: 0)
        computeEncoder6.setBuffer(rrOldBuffer, offset: 0, index: 1)

        let gridSize6 = MTLSize(width: 1, height: 1, depth: 1)
        let threadgroupSize6 = MTLSize(width: 1, height: 1, depth: 1)
        computeEncoder6.dispatchThreads(gridSize6, threadsPerThreadgroup: threadgroupSize6)
        computeEncoder6.endEncoding()
    }
    
    func configureSum(_ commandBuffer: MTLCommandBuffer, _ inputBuffer: MTLBuffer, _ outputBuffer: MTLBuffer) {
        let numElements = n*n*n
        let groupSize = 1024
        
        let numThreadGroups1 = (numElements + groupSize - 1) / groupSize
        let numThreadGroups2 = (numThreadGroups1 + groupSize - 1) / groupSize

        // First kernel for parallel reduction
        let computeEncoder = commandBuffer.makeComputeCommandEncoder()!
        computeEncoder.setComputePipelineState(sumPipeline)
        computeEncoder.setBuffer(inputBuffer, offset: 0, index: 0)
        computeEncoder.setBuffer(partialSumsBuffer, offset: 0, index: 1)
        
        let gridSize = MTLSize(width: numElements, height: 1, depth: 1)
        let threadgroupSize = MTLSize(width: groupSize, height: 1, depth: 1)
        computeEncoder.dispatchThreads(gridSize, threadsPerThreadgroup: threadgroupSize)
        computeEncoder.endEncoding()
        
        // Second kernel for parallel reduction, if needed
        var finalInputBuffer = partialSumsBuffer
        var finalThreadCount = numThreadGroups1
        if numThreadGroups2 > 1 {
            let computeEncoder2 = commandBuffer.makeComputeCommandEncoder()!
            computeEncoder2.setComputePipelineState(sumPipeline)
            computeEncoder2.setBuffer(partialSumsBuffer, offset: 0, index: 0)
            computeEncoder2.setBuffer(partialSums2Buffer, offset: 0, index: 1)
            let gridSize2 = MTLSize(width: numThreadGroups1, height: 1, depth: 1)
            let threadgroupSize2 = MTLSize(width: groupSize, height: 1, depth: 1)
            computeEncoder2.dispatchThreads(gridSize2, threadsPerThreadgroup: threadgroupSize2)
            computeEncoder2.endEncoding()
            
            finalInputBuffer = partialSums2Buffer
            finalThreadCount = numThreadGroups2
        }
        
        // final kernel to sum partial sums
        let finalEncoder = commandBuffer.makeComputeCommandEncoder()!
        finalEncoder.setComputePipelineState(sumPipeline)
        finalEncoder.setBuffer(finalInputBuffer, offset: 0, index: 0)
        finalEncoder.setBuffer(outputBuffer, offset: 0, index: 1)
        
        let finalGridSize = MTLSize(width: finalThreadCount, height: 1, depth: 1)
        finalEncoder.dispatchThreads(finalGridSize, threadsPerThreadgroup: MTLSize(width: finalThreadCount, height: 1, depth: 1))
        finalEncoder.endEncoding()
    }
    
    func setBufferData() {
        let numElements = n*n*n
        let groupSize = 1024
        
        let numThreadGroups1 = (numElements + groupSize - 1) / groupSize
        let numThreadGroups2 = (numThreadGroups1 + groupSize - 1) / groupSize

        // Create buffers
        aBuffer = device.makeBuffer(bytes: A, length: MemoryLayout<Float>.size * numElements, options: .storageModeShared)!
        aPlusIBuffer = device.makeBuffer(bytes: AplusI, length: MemoryLayout<Float>.size * numElements, options: .storageModeShared)!
        aPlusJBuffer = device.makeBuffer(bytes: AplusJ, length: MemoryLayout<Float>.size * numElements, options: .storageModeShared)!
        aPlusKBuffer = device.makeBuffer(bytes: AplusK, length: MemoryLayout<Float>.size * numElements, options: .storageModeShared)!
        xBuffer = device.makeBuffer(bytes: Array(repeating: 0, count: n*n*n), length: MemoryLayout<Float>.size * numElements, options: .storageModeShared)!
        bBuffer = device.makeBuffer(bytes: b, length: MemoryLayout<Float>.size * numElements, options: .storageModeShared)!
        rBuffer = device.makeBuffer(length: MemoryLayout<Float>.size * numElements, options: .storageModeShared)!
        rrBuffer = device.makeBuffer(length: MemoryLayout<Float>.size * numElements, options: .storageModeShared)!
        rrOldBuffer = device.makeBuffer(length: MemoryLayout<Float>.size, options: .storageModeShared)!
        pBuffer = device.makeBuffer(length: MemoryLayout<Float>.size * numElements, options: .storageModeShared)!
        ApBuffer = device.makeBuffer(length: MemoryLayout<Float>.size * numElements, options: .storageModeShared)!
        pApBuffer = device.makeBuffer(length: MemoryLayout<Float>.size * numElements, options: .storageModeShared)!

        partialSumsBuffer = device.makeBuffer(length: MemoryLayout<Float>.size * numThreadGroups1, options: .storageModeShared)!
        partialSums2Buffer = device.makeBuffer(length: MemoryLayout<Float>.size * numThreadGroups2, options: .storageModeShared)!
    }
}
