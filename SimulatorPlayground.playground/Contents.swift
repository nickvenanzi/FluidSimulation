import PlaygroundSupport

var A: [[Float]] = [[6, 2, 1, 3], [2, 5, 2, 1], [1, 2, 4, 2], [3, 1, 2, 7]]
var b: [Float] = [7, -8000, 9, 10]
var x: [Float] = [0, 0, 0, 0]
let RESIDUAL_THRESHOLD: Float = 1e-6

func mul(_ A: [[Float]], _ b: [Float]) -> [Float] {
    var result: [Float] = Array(repeating: 0, count: A.count)
    for i in 0..<A.count {
        var sum: Float = 0
        for j in 0..<b.count {
            sum += A[i][j] * b[j]
        }
        result[i] = sum
    }
    return result
}

func dot(_ a: [Float], _ b: [Float]) -> Float {
    var sum: Float = 0
    for i in 0..<a.count {
        sum += a[i] * b[i]
    }
    return sum
}

func sum(_ a: [Float], _ b: [Float]) -> [Float] {
    var result: [Float] = Array(repeating: 0.0, count: a.count)
    for i in 0..<a.count {
        result[i] = a[i] + b[i]
    }
    return result
}

func sub(_ a: [Float], _ b: [Float]) -> [Float] {
    var result: [Float] = Array(repeating: 0.0, count: a.count)
    for i in 0..<a.count {
        result[i] = a[i] - b[i]
    }
    return result
}

func sumC(_ a: [Float], _ b: [Float], _ C: Float) -> [Float] {
    var result: [Float] = Array(repeating: 0.0, count: a.count)
    for i in 0..<a.count {
        result[i] = a[i] + (C * b[i])
    }
    return result
}

func subC(_ a: [Float], _ b: [Float], _ C: Float) -> [Float] {
    var result: [Float] = Array(repeating: 0.0, count: a.count)
    for i in 0..<a.count {
        result[i] = a[i] - (C * b[i])
    }
    return result
}

func checkResidual(_ r: [Float]) -> Bool {
    var max: Float = 0
    for i in 0..<r.count {
        if abs(r[i]) > max {
            max = abs(r[i])
        }
    }
    return max < RESIDUAL_THRESHOLD
}

func check(_ r: [Float], _ iterations: Int) -> Bool {
    if checkResidual(r) {
        print("Threshold met after \(iterations) iterations.  r: \(r)")
        return true
    } else {
        print("Iteration \(iterations): \nr: \(r)")
        return false
    }
}

var iterations: Int = 0
var r = sub(b, mul(A, x))
var r_dot_r = dot(r, r)

var p = r

if check(r, iterations) { PlaygroundPage.current.finishExecution() }

while iterations < 100 {
    let Ap = mul(A, p)
    let alpha = r_dot_r / dot(p, Ap)
    x = sumC(x, p, alpha)
    r = subC(r, Ap, alpha)

    iterations += 1
    if check(r, iterations) { PlaygroundPage.current.finishExecution() }

    let new_r_dot_r = dot(r, r)
    let beta = new_r_dot_r / r_dot_r

    r_dot_r = new_r_dot_r

    p = sumC(r, p, beta)
}

print("Reached max iterations without converging")
