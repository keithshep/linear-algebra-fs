// determines if the given number is near 0
let nearZero = 1e-8
let isNearZero x = x < nearZero && x > -(nearZero)

let reorderArray (array : 'a array) (order : int array) =
    let len = Array.length array
    assert (Array.length order = len)
    Array.init len (fun i -> array.[order.[i]])

let reorderMatrixRows matrix rowOrder =
    let rowCount = Array2D.length1 matrix
    let colCount = Array2D.length2 matrix
    assert (Array.length rowOrder = rowCount)
    
    Array2D.init rowCount colCount (fun row col -> matrix.[rowOrder.[row], col])

// applies f to each corresponding pair of elements in m1 and m2
let matZipWith (f : 'a -> 'b -> 'c) (m1 : 'a[,]) (m2 : 'b[,]) =
    let rowCount = Array2D.length1 m1
    let colCount = Array2D.length2 m1
    
    assert (Array2D.length1 m2 = rowCount && Array2D.length2 m2 = colCount)
    Array2D.init rowCount colCount (fun i j -> f m1.[i, j] m2.[i, j])

// adds two matrices
let matAdd = matZipWith (+)

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

// the identity matrix has 1's along the diagonal and 0's elsewhere
let identityMatrix n = Array2D.init n n (fun i j -> if i = j then 1.0 else 0.0)

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

let mean vec = Array.sum vec / double (Array.length vec)

let variance vec =
    let mu = mean vec
    let sumSqDiffs = ref 0.0
    for x in vec do sumSqDiffs := !sumSqDiffs + (mu - x) ** 2.0
    !sumSqDiffs

let stdDev vec = sqrt (variance vec)

// standardize the vector by mean centering and dividing by
// standard deviation
let standardize vec =
    let mu = mean vec
    let vecStdDev = stdDev vec
    [| for x in vec -> (x - mu) / vecStdDev |]

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

let isSquareMat m = Array2D.length1 m = Array2D.length2 m

let trace m =
    assert (isSquareMat m)
    let diagSum = ref m.[0, 0]
    for i = 1 to Array2D.length1 m - 1 do
        diagSum := !diagSum + m.[i, i]
    !diagSum

let transpose m =
    Array2D.init (Array2D.length2 m) (Array2D.length1 m) (fun i j -> m.[j, i])

let matInnerProd m1 m2 =
    // TODO: look further into
    //       http://en.wikipedia.org/wiki/Trace_%28linear_algebra%29#Inner_product
    //       for info about the conjugate transpose (in the case of imaginaries)
    trace (matMult (transpose m1) m2)

// calculate cos theta using the given vectors
let vecCosine v1 v2 = dotProd v1 v2 / (vecNorm v1 * vecNorm v2)

// coefficient of correlation is the same thing as `vecCosine`
let coefOfCorr = vecCosine

// two vectors are orthogonal if thier dot product is zero
let orthogonal v1 v2 = isNearZero (dotProd v1 v2)

// extract the given row number from the given matrix
let matRow m row = Array.init (Array2D.length2 m) (fun col -> m.[row, col])

// extract the given column number
let matColumn (m : 'a [,]) (col : int) = Array.init (Array2D.length1 m) (fun row -> m.[row, col])

// The fourier coefficients are the scalar values that we have to
// multiply each member of the `orthonormalSet` by in order
// to recover `vec`
let fourierCoefficients orthonormalBasis v =
    let coefForRow row = dotProd (matRow orthonormalBasis row) v
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

(*
Note: basis definition from wikipedia
In linear algebra, a basis is a set of vectors that, in a linear combination,
can represent every vector in a given vector space or free module, and such
that no element of the set can be represented as a linear combination of the
others. In other words, a basis is a linearly independent spanning set.
*)

// Note: Orthonormal sets are mutually orthogonal and each vector has unit length

let setColumn m col colVector =
    let nRow = Array2D.length1 m
    assert (nRow = Array.length colVector)
    for row = 0 to nRow - 1 do
        m.[row, col] <- colVector.[row]

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
        let outerCol = matColumn m outerColIndex
        
        // inner column loop
        for innerColIndex = 0 to outerColIndex - 1 do
            let innerCol = matColumn m innerColIndex
            let currDotProd = dotProd outerCol innerCol
            for row = 0 to nRow - 1 do
                outerCol.[row] <- outerCol.[row] - (innerCol.[row] * currDotProd)
        setColumn m outerColIndex (normalizeVector outerCol)
    m

// This is one way to QR factorize with gram schmidt othroganalization
let qrFactorizeWithGramSchmidt m =
    let q = gramSchmidtOrth m
    
    // R can also be calculated at the same time
    // as the GS orthoganalization. See wikipedia
    let r = matMult (transpose q) m
    
    (q, r)

let columnVector (v : 'a []) = Array2D.init (Array.length v) 1 (fun row _ -> v.[row])

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

// determine if the given matrix is square
let isSquare m = Array2D.length1 m = Array2D.length2 m

// use guassian elimination to zero out the coefficients under the diagonal
// Returns the resulting coefficient matrix and right-hand-side values and
// orderings. Rather than leaving zeroes in the bottom diagonal of the
// coefficient we fill it in with the multipliers that were used to zero out the
// lower diagonal
let private gaussianEliminationInternal coefMat rhsVec =
    let coefMat = Array2D.copy coefMat
    let rhsVec = Array.copy rhsVec
    let size = Array2D.length1 coefMat

    // do some sanity checking
    if size = 0 then failwith "coefMat cannot be empty"
    if not (isSquare coefMat) then
        invalidArg "coefMat" "the coefficient matrix must be a square matrix"
    if Array.length rhsVec <> size then
        let errorMsg =
            sprintf
                "There are %i rows in coefMat and %i elements in rhsVec. These \
                 values should be the same."
                size
                (Array.length rhsVec)
        failwith errorMsg

    let rowOrdering = [|0 .. size - 1|]
    for diagIndex = 0 to size - 2 do
        // perform "partial pivoting" which bubbles the max absolute coefficient
        // to the top (improves algorithm's accuracy)
        let maxIndex = ref diagIndex
        let maxAbsCoef = ref (abs coefMat.[rowOrdering.[diagIndex], diagIndex])
        for j = diagIndex + 1 to size - 1 do
            let currAbsVal = abs coefMat.[rowOrdering.[j], diagIndex]
            if currAbsVal > !maxAbsCoef then
                maxAbsCoef := currAbsVal
                maxIndex := j
        if isNearZero !maxAbsCoef then
            let errorMsg =
                sprintf
                    "the matrix is singular: could not find a non-zero \
                     coefficient under the diagonal at position %i"
                    rowOrdering.[diagIndex]
            failwith errorMsg

        // now swap the max row with the current
        if !maxIndex <> diagIndex then
            let tmp = rowOrdering.[diagIndex]
            rowOrdering.[diagIndex] <- rowOrdering.[!maxIndex]
            rowOrdering.[!maxIndex] <- tmp

        // now "zero out" the coefficients below the diagonal
        let diagCoef = coefMat.[rowOrdering.[diagIndex], diagIndex]
        for row = diagIndex + 1 to size - 1 do
            let orderedRow = rowOrdering.[row]
            let currCoef = coefMat.[orderedRow, diagIndex]
            let zeroFactor = currCoef / diagCoef
            if not (isNearZero zeroFactor) then
                for col = diagIndex + 1 to size - 1 do
                    coefMat.[orderedRow, col] <-
                        coefMat.[orderedRow, col] -
                        (zeroFactor * coefMat.[rowOrdering.[diagIndex], col])
                rhsVec.[orderedRow] <-
                    rhsVec.[orderedRow] - (zeroFactor * rhsVec.[rowOrdering.[diagIndex]])
            coefMat.[orderedRow, diagIndex] <- zeroFactor
    
    (coefMat, rhsVec, rowOrdering)

let gaussianElimination coefMat rhsVec =
    let (gaussElimCoefMat, gaussElimRHSVec, rowOrdering) = gaussianEliminationInternal coefMat rhsVec
    (reorderMatrixRows gaussElimCoefMat rowOrdering, reorderArray gaussElimRHSVec rowOrdering)

// Solve the linear equations using gaussian elimination and back
// substitution
let solveWithGaussAndBackSub coefMat rhsVec =
    let (gaussElimCoefMat, gaussElimRHSVec) = gaussianElimination coefMat rhsVec
    backSubstituteUpper gaussElimCoefMat gaussElimRHSVec

// See page 314 of "Matrix Analysis and Applied Linear Algebra"
let solveLeastSquares a b =
    let (q, r) = qrFactorizeWithGramSchmidt a
    // R x = trans(Q) b
    // two of the transposes are just for switching to/from single
    // column matrices/vectors
    let rightSide = matColumn (matMult (transpose q) (columnVector b)) 0
    backSubstituteUpper r rightSide

// EVERYTHING BELOW THIS LINE IS JUST FOR QUICK AND DIRTY TESTING

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
    printfn "%A" (matAdd (array2D [[1.0; 2.0; 3.0]]) (array2D [[0.1; 0.2; 0.3]]))
    printfn "Matrix Multiplication:"
    printfn "%A" (matMult testMat1 testMat2)
    printfn "Identity Mat 10:"
    printfn "%A" (identityMatrix 10)
    printfn "Least Squares:"
    printfn "%A" (solveLeastSquares x y)
    printfn "Gaussian/Back-substitution solution for L"
    printfn "%A" (solveWithGaussAndBackSub lCoefMatrix lRHSVector)
    printfn "Gaussian/Back-substitution solution for M"
    printfn "%A" (solveWithGaussAndBackSub mCoefMatrix mRHSVector)

