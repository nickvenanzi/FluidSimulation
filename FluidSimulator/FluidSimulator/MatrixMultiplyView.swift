//
//  MatrixMultiplyView.swift
//  FluidSimulator
//
//  Created by Nick Venanzi on 9/29/24.
//

import SwiftUI
import Metal

struct MatrixMultiplyView: View {
    @State private var result: [Float] = []
    
    var body: some View {
        VStack {
            if result.isEmpty {
                Text("Tap to compute matrix multiplication")
                    .padding()
                    .onAppear(perform: computeMatrixMultiplication)
            } else {
                Text("Matrix Multiplication Result:")
                ForEach(0..<result.count, id: \.self) { val in
                    HStack {
                        Text(String(format: "%.2f", result[val]))
                            .frame(width: 50, height: 50)
                    }
                }
            }
        }
    }
    
    func generateMatrix(_ n: Int) -> [Float] {
        var matrix: [Float] = Array(repeating: 0, count: n*n)
        for i in 0..<(n*n) {
            matrix[i] = Float.random(in: 0..<1)
        }
        return matrix
    }
    
    func generateArray(_ n: Int) -> [Float] {
        var array: [Float] = Array(repeating: 0, count: n)
        for i in 0..<n {
            array[i] = Float.random(in: 0..<1)
        }
        return array
    }
    
    func computeMatrixMultiplicationCPU(_ A: [Float], _ b: [Float]) -> [Float] {
        var result: [Float] = Array(repeating: 0, count: b.count)
        for row in 0..<b.count {
            var sum: Float = 0
            for i in 0..<b.count {
                sum += A[row * b.count + i] * b[i]
            }
            result[row] = sum
        }
        return result
    }
    
    func computeMatrixMultiplication() {
        let n = 1000
        let A: [Float] = generateMatrix(n)
        let b: [Float] = generateArray(n)
        
        var startTime = CFAbsoluteTimeGetCurrent()
        let _ = computeMatrixMultiplicationCPU(A, b)
        var endTime = CFAbsoluteTimeGetCurrent()  // Start the timer
        print("CPU Time elapsed: \(endTime - startTime) seconds")


        // GPU-based Matrix Multiplication
        startTime = CFAbsoluteTimeGetCurrent()  // Start the timer
        if let result = performMatrixMultiplication(A, b) {
            endTime = CFAbsoluteTimeGetCurrent()  // Start the timer
            print("GPU Time elapsed: \(endTime - startTime) seconds")
            self.result = result
        }
    }
    
    func performMatrixMultiplication(_ A: [Float], _ b: [Float]) -> [Float]? {
        // Get the Metal device and command queue
        guard let device = MTLCreateSystemDefaultDevice(),
              let commandQueue = device.makeCommandQueue(),
              let library = device.makeDefaultLibrary(),
              let kernelFunction = library.makeFunction(name: "matrix_multiplication"),
              let pipelineState = try? device.makeComputePipelineState(function: kernelFunction)
        else {
            print("Failed to set up Metal")
            return nil
        }
        
        // Create buffers
        let matrixADataSize = A.count * MemoryLayout<Float>.size
        let arrayBDataSize = b.count * MemoryLayout<Float>.size
        let resultDataSize = b.count * MemoryLayout<Float>.size
        
        guard let bufferA = device.makeBuffer(bytes: A, length: matrixADataSize, options: []),
              let bufferB = device.makeBuffer(bytes: b, length: arrayBDataSize, options: []),
              let resultBuffer = device.makeBuffer(length: resultDataSize, options: [])
        else {
            print("Failed to create buffers")
            return nil
        }
        
        // Create command buffer and compute command encoder
        guard let commandBuffer = commandQueue.makeCommandBuffer(),
              let computeEncoder = commandBuffer.makeComputeCommandEncoder() else {
            print("Failed to create command buffer/encoder")
            return nil
        }
        
        computeEncoder.setComputePipelineState(pipelineState)
        computeEncoder.setBuffer(bufferA, offset: 0, index: 0)
        computeEncoder.setBuffer(bufferB, offset: 0, index: 1)
        computeEncoder.setBuffer(resultBuffer, offset: 0, index: 2)
        
        // Send matrix dimensions to the GPU
        // Create the buffer for the struct
        var dimensions = MatrixDimensions(n: UInt32(b.count))
        let dimensionsBuffer = device.makeBuffer(bytes: &dimensions, length: MemoryLayout<MatrixDimensions>.stride, options: [])

        // Set the buffer on the Metal command encoder
        computeEncoder.setBuffer(dimensionsBuffer, offset: 0, index: 3)
        
        // Set threadgroup size and dispatch
        let gridSize = MTLSize(width: b.count, height: 1, depth: 1)
        let threadGroupSize = MTLSize(width: b.count, height: 1, depth: 1)
        
        computeEncoder.dispatchThreads(gridSize, threadsPerThreadgroup: threadGroupSize)
        computeEncoder.endEncoding()
        
        // Commit command buffer
        commandBuffer.commit()
        commandBuffer.waitUntilCompleted()
        
        // Get the result
        let resultPointer = resultBuffer.contents().bindMemory(to: Float.self, capacity: b.count)
        return Array(UnsafeBufferPointer(start: resultPointer, count: b.count))
    }
}

// Define a Swift struct that matches the Metal struct
struct MatrixDimensions {
    var n: UInt32
}
