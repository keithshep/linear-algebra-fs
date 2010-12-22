module LinearAlgebra

module MC = MathCore
module U = Utilities

////////////////////////////////////////////////////////////////////////////////
// BASIC LINEAR ALGEBRA FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

// the identity matrix has 1's along the diagonal and 0's elsewhere
let inline identityMatrix n = Array2D.init n n (fun i j -> if i = j then 1.0 else 0.0)

// adds two matrices
let inline matAdd m1 m2 = U.matZipWith (+) m1 m2

// Multiplies the two given matrices
// The wikipedia page on matrix multiplication is great: http://en.wikipedia.org/wiki/Matrix_multiplication
// 
// The most important points are:
// * The row count of the resulting matrix is determined by the row count of matrix1
// * The column count of the resulting matrix is determined by the column count of matrix2
// * element[i, j] of the resulting matrix is calculated by taking the dot product
//   of the ith row vector from matrix1 and the jth column from matrix2. Also note
//   that because these two vectors must be the same length, the column count of
//   matrix1 must equal the row count of matrix2.
let matMult m1 m2 =
    let elemCount = Array2D.length2 m1
    assert (elemCount = Array2D.length1 m2)
    
    if elemCount = 0 then
        Array2D.zeroCreate 0 0
    else
        // dot product of a row from m1 and column from m2
        let dotProdRowCol row col =
            let rowColDotProd = ref (m1.[row, 0] * m2.[0, col])
            for i = 1 to elemCount - 1 do
                rowColDotProd := m1.[row, i] * m2.[i, col] + !rowColDotProd
            !rowColDotProd
        
        Array2D.init (Array2D.length1 m1) (Array2D.length2 m2) dotProdRowCol

// Calculate the norm of the given vector.
// Math notation for norm is ||X||
// TODO: check for complex numbers
let vecPNorm p vec =
    let norm = ref 0.0
    for v in vec do
        norm := v ** p + !norm
    !norm ** (1.0 / p)

// Calculate the norm of the given vector. TODO check for complex numbers
let vecNorm vec = vecPNorm 2.0 vec

// As P approaches infinity the vector p-norm becomes the max absolute value in
// the vector
let vecInfNorm vec =
    let infNorm = ref 0.0
    for v in vec do
        let absV = abs v
        if absV > !infNorm then infNorm := absV
    !infNorm

// Normalize `vector`. Geometricly this means have it point in the same
// direction, but scale its magnitude to 1
let normalizeVector vec =
    let norm = vecNorm vec
    Array.map (fun x -> x / norm) vec

// standardize the vector by mean centering and dividing by
// standard deviation
let standardize vec =
    let mu = MC.mean vec
    let vecStdDev = MC.stdDev vec
    [| for x in vec -> (x - mu) / vecStdDev |]

/// The dot product is the inner product for euclidean spaces. Notations include
/// <x,y> <x|y> and x . y
let dotProd vec1 vec2 =
    let elemCount = Array.length vec1
    assert (elemCount = Array.length vec2)
    let result = ref 0.0
    for i = 0 to elemCount - 1 do
        result := !result + (vec1.[i] * vec2.[i])
    !result

let vecCorrelation vec1 vec2 =
    let stdVec1 = standardize vec1
    let stdVec2 = standardize vec2
    dotProd stdVec1 stdVec2 / (vecNorm stdVec1 * vecNorm stdVec2)

let trace m =
    assert (U.isSquareMat m)
    let diagSum = ref m.[0, 0]
    for i = 1 to Array2D.length1 m - 1 do
        diagSum := !diagSum + m.[i, i]
    !diagSum

let matInnerProd m1 m2 =
    // TODO: look further into
    //       http://en.wikipedia.org/wiki/Trace_%28linear_algebra%29#Inner_product
    //       for info about the conjugate transpose (in the case of imaginaries)
    trace (matMult (U.transpose m1) m2)

// calculate cos theta using the given vectors
let vecCosine v1 v2 = dotProd v1 v2 / (vecNorm v1 * vecNorm v2)

// coefficient of correlation is the same thing as `vecCosine`
let coefOfCorr = vecCosine

// two vectors are orthogonal if thier dot product is zero
let inline orthogonal v1 v2 = MC.isNearZero (dotProd v1 v2)

let isOrthonormal m =
    let nCol = Array2D.length2 m
    U.all <| seq {
        for i = 0 to nCol do
            // all column norms should be 1
            let innerCol = U.matColumn m i
            yield (MC.isNear 1.0 (vecNorm innerCol))

            // all column pairs should be orthogonal
            for j = i + 1 to nCol do
                let outerCol = U.matColumn m j
                yield (orthogonal innerCol outerCol)
    }

let inline isOrthonormalBasis m = U.isSquareMat m && isOrthonormal m

////////////////////////////////////////////////////////////////////////////////
// Guassian Elimination Etc.
////////////////////////////////////////////////////////////////////////////////

// solve a matrix that is zeroed out below the diagonal using
// back substitution
let backSubstituteUpper (upperTriangleMat : double [,]) (rhsVec : double []) =
    let n = Array2D.length1 upperTriangleMat
    assert (Array2D.length2 upperTriangleMat = n)
    
    let solution = (Array.zeroCreate n : double array)
    for i = n - 1 downto 0 do
        let mutable sum = 0.0
        for j = n - 1 downto i + 1 do
            sum <- sum + (upperTriangleMat.[i, j] * solution.[j])
        solution.[i] <- ((rhsVec.[i] - sum) / upperTriangleMat.[i, i] : double)
    solution

// solve a matrix that is zeroed out above the diagonal using
// back substitution
let backSubstituteLower (lowerTriangleMat : double [,]) (rhsVec : double []) =
    let n = Array2D.length1 lowerTriangleMat
    assert (Array2D.length2 lowerTriangleMat = n)
    
    let solution = (Array.zeroCreate n : double array)
    for i = 0 to n - 1 do
        let sum = ref 0.0
        for j = 0 to i - 1 do
            sum := !sum + (lowerTriangleMat.[i, j] * solution.[j])
        solution.[i] <- (rhsVec.[i] - !sum) / lowerTriangleMat.[i, i]
    solution

/// <summary>
/// Use guassian elimination with pivoting to calculate LU decomposition for the
/// given matrix. Since pivoting changes row ordering the reordered indices are
/// also returned.
/// </summary>
/// <param name="m">the coefficient matrix (must be square)</param>
/// <returns>
/// A tuple containing: (luMatrix, rowOrder). The upper triangle of luMatrix
/// (including diagonal) contains U while the lower triangle of luMatrix
/// contains L
/// </returns>
let luDecompose matrix =
    let matrix = Array2D.copy matrix
    let size = Array2D.length1 matrix

    // do some sanity checking
    if size = 0 then failwith "matrix cannot be empty"
    if not (U.isSquareMat matrix) then
        invalidArg "matrix" "the coefficient matrix must be a square matrix"

    let rowOrder = [|0 .. size - 1|]
    for diagIndex = 0 to size - 2 do
        // perform "partial pivoting" which bubbles the max absolute coefficient
        // to the top (improves algorithm's accuracy)
        let maxIndex = ref diagIndex
        let maxAbsCoef = ref (abs matrix.[rowOrder.[diagIndex], diagIndex])
        for j = diagIndex + 1 to size - 1 do
            let currAbsVal = abs matrix.[rowOrder.[j], diagIndex]
            if currAbsVal > !maxAbsCoef then
                maxAbsCoef := currAbsVal
                maxIndex := j
        if MC.isNearZero !maxAbsCoef then
            let errorMsg =
                sprintf
                    "the matrix is singular: could not find a non-zero \
                     coefficient under the diagonal at position %i"
                    rowOrder.[diagIndex]
            failwith errorMsg

        // now swap the max row with the current
        if !maxIndex <> diagIndex then
            let tmp = rowOrder.[diagIndex]
            rowOrder.[diagIndex] <- rowOrder.[!maxIndex]
            rowOrder.[!maxIndex] <- tmp

        // now "zero out" the coefficients below the diagonal
        let diagCoef = matrix.[rowOrder.[diagIndex], diagIndex]
        for row = diagIndex + 1 to size - 1 do
            let orderedRow = rowOrder.[row]
            let currCoef = matrix.[orderedRow, diagIndex]
            let zeroFactor = currCoef / diagCoef
            if not (MC.isNearZero zeroFactor) then
                for col = diagIndex + 1 to size - 1 do
                    matrix.[orderedRow, col] <-
                        matrix.[orderedRow, col] -
                        (zeroFactor * matrix.[rowOrder.[diagIndex], col])
            matrix.[orderedRow, diagIndex] <- zeroFactor
    
    (U.reorderMatrixRows matrix rowOrder, rowOrder)

let solveLU (luMatrix : double [,]) (rhsVec : double []) =
    let size = Array.length rhsVec
    for col = 0 to size - 2 do
        for row = col + 1 to size - 1 do
            rhsVec.[row] <- rhsVec.[row] - (rhsVec.[col] * luMatrix.[row, col])
    backSubstituteUpper luMatrix rhsVec

let solveManyWithGaussAndBackSub coefMat rhsVecs =
    let (luMatrix, rowOrder) = luDecompose coefMat
    List.map (fun rhs -> U.reorderArray rhs rowOrder |> solveLU luMatrix) rhsVecs

// Solve the linear equations using gaussian elimination and back
// substitution
let inline solveWithGaussAndBackSub coefMat rhsVec =
    solveManyWithGaussAndBackSub coefMat [rhsVec] |> List.head

////////////////////////////////////////////////////////////////////////////////
// Solving For Least Squares
////////////////////////////////////////////////////////////////////////////////

// Note: Orthonormal sets are mutually orthogonal and each vector has unit length

(*
From "Matrix Analysis and Applied Linear Algebra"

Objective: Use B to construct an orthonormal basis O = {u1, u2 , ... , un }
           for S.
Strategy: Construct O sequentially so that Ok = {u1 , u2 , ... , uk } is an
          orthonormal basis for Sk = span {x1 , x2 , ... , xk } for k = 1, ... , n.

in our function B is called `arbitraryBasis`. Use the "modified" versions
of this function for something that is more stable.
*)
let gramSchmidtOrth m =
    let nRow = Array2D.length1 m
    let nCol = Array2D.length2 m
    let m = Array2D.copy m
    
    // outer column loop
    for outerColIndex = 0 to nCol - 1 do
        let outerCol = U.matColumn m outerColIndex
        
        // inner column loop
        for innerColIndex = 0 to outerColIndex - 1 do
            let innerCol = U.matColumn m innerColIndex
            let currDotProd = dotProd outerCol innerCol
            for row = 0 to nRow - 1 do
                outerCol.[row] <- outerCol.[row] - (innerCol.[row] * currDotProd)
        U.setMatrixColumn m outerColIndex (normalizeVector outerCol)
    m

// This is one way to QR factorize with gram schmidt othroganalization
let qrFactorizeWithGramSchmidt m =
    let q = gramSchmidtOrth m
    
    // R can also be calculated at the same time
    // as the GS orthoganalization. See wikipedia
    let r = matMult (U.transpose q) m
    
    (q, r)

// See page 314 of "Matrix Analysis and Applied Linear Algebra"
let solveLeastSquares a b =
    let (q, r) = qrFactorizeWithGramSchmidt a
    // R x = trans(Q) b
    let rightSide = U.matColumn (matMult (U.transpose q) (U.columnVector b)) 0
    backSubstituteUpper r rightSide

