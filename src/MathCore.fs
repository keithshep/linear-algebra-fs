module MathCore

// determines if the given number is near 0
let nearZero = 1e-8
let inline isNearZero x = x < nearZero && x > -(nearZero)
let inline isNear num x = isNearZero (x - num)

let mean vec = Array.sum vec / double (Array.length vec)

let variance vec =
    let mu = mean vec
    let sumSqDiffs = ref 0.0
    for x in vec do sumSqDiffs := !sumSqDiffs + (mu - x) ** 2.0
    !sumSqDiffs

let stdDev vec = sqrt (variance vec)


