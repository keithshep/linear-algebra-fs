module Utilities

let inline transpose m =
    Array2D.init (Array2D.length2 m) (Array2D.length1 m) (fun i j -> m.[j, i])

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

let inline isSquareMat m = Array2D.length1 m = Array2D.length2 m

let inline columnVector (v : 'a []) = Array2D.init (Array.length v) 1 (fun row _ -> v.[row])

let inline setMatrixColumn m col colVector =
    let nRow = Array2D.length1 m
    assert (nRow = Array.length colVector)
    for row = 0 to nRow - 1 do
        m.[row, col] <- colVector.[row]

// extract the given row number from the given matrix
let inline matRow m row = Array.init (Array2D.length2 m) (fun col -> m.[row, col])

// extract the given column number
let inline matColumn m col = Array.init (Array2D.length1 m) (fun row -> m.[row, col])

