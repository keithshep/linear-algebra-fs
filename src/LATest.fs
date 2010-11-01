module LATest

module LA = LinearAlgebra

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

let main =
    printfn "Matrix Addition:"
    printfn "%A" (LA.matAdd (array2D [[1.0; 2.0; 3.0]]) (array2D [[0.1; 0.2; 0.3]]))
    printfn "Matrix Multiplication:"
    printfn "%A" (LA.matMult testMat1 testMat2)
    printfn "Identity Mat 10:"
    printfn "%A" (LA.identityMatrix 10)
    printfn "Least Squares:"
    printfn "%A" (LA.solveLeastSquares x y)
    printfn "Gaussian/Back-substitution solution for L"
    printfn "%A" (LA.solveWithGaussAndBackSub lCoefMatrix lRHSVector)
    printfn "Gaussian/Back-substitution solution for M"
    printfn "%A" (LA.solveWithGaussAndBackSub mCoefMatrix mRHSVector)

