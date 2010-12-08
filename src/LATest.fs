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
        let sum = ref 0.0
        for row = 0 to Array2D.length1 orthonormalBasis - 1 do
            sum := !sum + coefs.[row] * orthonormalBasis.[row, col]
        !sum
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

[<EntryPoint>]
let main _ =
    printfn "Matrix Addition:"
    printfn "%A" (LA.matAdd (array2D [[1.0; 2.0; 3.0]]) (array2D [[0.1; 0.2; 0.3]]))
    printfn "Matrix Multiplication:"
    printfn "%A" (LA.matMult testMat1 testMat2)
    printfn "Identity Mat 10:"
    printfn "%A" (LA.identityMatrix 10)
    printfn "Least Squares:"
    printfn "%A" (LA.solveLeastSquares x y)
    printfn "Gaussian/Back-substitution solution for Lx=l"
    printfn "%A" (LA.solveWithGaussAndBackSub lCoefMatrix lRHSVector)
    printfn "Gaussian/Back-substitution solution for Mx=m"
    printfn "%A" (LA.solveWithGaussAndBackSub mCoefMatrix mRHSVector)
    printfn "Gaussian/Back-substitution solution for Mx=m, Mx=n"
    printfn "%A" (LA.solveManyWithGaussAndBackSub mCoefMatrix [mRHSVector; nRHSVector])
    printfn "Gaussian/Back-substitution solution for Mx=n"
    printfn "%A" (LA.solveWithGaussAndBackSub mCoefMatrix nRHSVector)

    // exit success
    0

