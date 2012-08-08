module LATest

module LA = LinearAlgebra
module U = Utilities

(*
Note: basis definition from wikipedia
In linear algebra, a basis is a set of vectors that, in a linear combination,
can represent every vector in a given vector space or free module, and such
that no element of the set can be represented as a linear combination of the
others. In other words, a basis is a linearly independent spanning set.
*)

// The fourier coefficients are the scalar values that we have to
// multiply each member of the `orthonormalSet` by in order
// to recover `vec`
let fourierCoefficients orthonormalBasis v =
    let coefForRow row = LA.dotProd (U.matRow orthonormalBasis row) v
    Array.init (Array2D.length1 orthonormalBasis) coefForRow

// uses the fourier expansion method to calculate the vector. in theory if
// orthonormalBasis is really what it claims to be than this function
// should just return vec (with some precision error)
let fourierExpand orthonormalBasis vec =
    let coefs = fourierCoefficients orthonormalBasis vec
    let expandCol col =
        let mutable sum = 0.0
        for row = 0 to Array2D.length1 orthonormalBasis - 1 do
            sum <- sum + coefs.[row] * orthonormalBasis.[row, col]
        sum
    Array.init (Array2D.length2 orthonormalBasis) expandCol

let testMat1 =
    Array2D.map double
        (array2D
            [[2; 0; -1; 1];
             [1; 2; 0; 1]])

let testMat2 =
    Array2D.map double
        (array2D
            [[1; 5; -7];
             [1; 1; 0];
             [0; -1; 1];
             [2; 0; 0]])

let x =
    Array2D.map double
        (array2D
            [[1; 1];
             [1; 2];
             [1; 3];
             [1; 4]])

let y = Array.map double [|6; 5; 7; 10|]

let lCoefMatrix =
    array2D
        [[2.0;  1.0; 1.0];
         [4.0;  3.0; 1.0];
         [-2.0; 2.0; 1.0]]

let lRHSVector = [|1.0; -1.0; 7.0|]

let mCoefMatrix =
    array2D
        [[2.0;  1.0; 1.0];
         [6.0;  3.0; 1.0];
         [-2.0; 2.0; 1.0]]

let mRHSVector = [|1.0; -1.0; 7.0|]

let nRHSVector = [|33.3; 2.0; 1.11|]

// see http://en.wikipedia.org/wiki/Ordinary_least_squares for data
let heights = [|1.47; 1.50; 1.52; 1.55; 1.57; 1.60; 1.63; 1.65; 1.68; 1.70; 1.73; 1.75; 1.78; 1.80; 1.83|]
let heightsSq = Array.map (fun x -> x ** 2.0) heights
let weights = [|52.21; 53.12; 54.48; 55.84; 57.20; 58.57; 59.93; 61.29; 63.11; 64.47; 66.28; 68.10; 69.92; 72.19; 74.46|]
let heightWeightPred =
    Array2D.init heights.Length 3 (
        fun (row : int) (col : int) ->
            match col with
            | 0 -> heights.[row]
            | 1 -> heightsSq.[row]
            | 2 -> 1.0
            | _ -> failwithf "bad column index %i" col
    )

[<EntryPoint>]
let main _ =
    let println() = System.Console.Out.WriteLine()
    
    //printfn "Matrix Addition:"
    //printfn "%A" (LA.matAdd (array2D [[1.0; 2.0; 3.0]]) (array2D [[0.1; 0.2; 0.3]]))
    printfn "=== Matrix Multiplication (calculate C for AB=C) ==="
    printfn "A="
    printfn "%A" testMat1
    printfn "B="
    printfn "%A" testMat2
    printfn "C="
    printfn "%A" (LA.matMult testMat1 testMat2)
    println()
    printfn "=== Identity Mat I_5 ==="
    printfn "%A" (LA.identityMatrix 5)
    println()
    printfn "=== Gaussian/Back-substitution solution (solve for x in for Mx=m) ==="
    printfn "%A" (LA.solveWithGaussAndBackSub mCoefMatrix mRHSVector)
    println()
    (*
    printfn "=== Gaussian/Back-substitution solution for Mx=m, Mx=n ==="
    printfn "%A" (LA.solveManyWithGaussAndBackSub mCoefMatrix [mRHSVector; nRHSVector])
    println()
    printfn "=== Gaussian/Back-substitution solution for Mx=n ==="
    printfn "%A" (LA.solveWithGaussAndBackSub mCoefMatrix nRHSVector)
    println()
    *)
    printfn "=== Least Squares (y=XB + e, solve for B that minimizes e^2) ==="
    printfn "X="
    printfn "%A" x
    printfn "y="
    printfn "%A" y
    printfn "B="
    let b = LA.solveLeastSquares x y
    printfn "%A" b
    printfn "XB="
    let xb = LA.matMult x (U.columnVector b)
    printfn "%A" xb
    println()
    printfn "=== Least Squares for Height and Weight from http://en.wikipedia.org/wiki/Ordinary_least_squares ==="
    printfn "X="
    printfn "%A" heightWeightPred
    printfn "y="
    printfn "%A" weights
    printfn "B="
    let b = LA.solveLeastSquares heightWeightPred weights
    printfn "%A" b
    printfn "XB="
    let xb = LA.matMult heightWeightPred (U.columnVector b)
    printfn "%A" xb
    println()
    printfn "=== Gaussian/Back-substitution solution for Lx=l ==="
    printfn "L="
    printfn "%A" lCoefMatrix
    printfn "l="
    printfn "%A" lRHSVector
    printfn "x="
    printfn "%A" (LA.solveWithGaussAndBackSub lCoefMatrix lRHSVector)

    // exit success
    0

