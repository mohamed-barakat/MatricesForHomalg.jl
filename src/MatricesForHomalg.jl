@doc """
Matrices for the homalg project
"""
module MatricesForHomalg

import Base: getindex

import AbstractAlgebra
import AbstractAlgebra: ZZ, QQ, matrix

export ZZ, QQ

TypeOfMatrixForHomalg = AbstractAlgebra.Generic.MatSpaceElem

## Constructors of homalg matrices

"""
    HomalgMatrix(L, r, c, R)

Construct a (r x c)-matrix over the ring R with L as list of entries

```jldoctest
julia> mat = HomalgMatrix([1,2,3,4,5,6], 2, 3, ZZ)
[1   2   3]
[4   5   6]
```
"""
function HomalgMatrix(entries, r, c, R)::TypeOfMatrixForHomalg
    return matrix(R, r, c, entries)
end

"""
    HomalgIdentityMatrix(r, R)

Construct a (r x r)-identity matrix over the ring R

```jldoctest
julia> mat = HomalgIdentityMatrix(3, ZZ)
[1   0   0]
[0   1   0]
[0   0   1]
```
"""
function HomalgIdentityMatrix(r, R)::TypeOfMatrixForHomalg
    return AbstractAlgebra.identity_matrix(R, r)
end

"""
    HomalgZeroMatrix(r, c, R)

Construct a (r x c)-zero matrix over the ring R

```jldoctest
julia> mat = HomalgZeroMatrix(3, 2, ZZ)
[0   0]
[0   0]
[0   0]
```
"""
function HomalgZeroMatrix(r, c, R)::TypeOfMatrixForHomalg
    return AbstractAlgebra.zero_matrix(R, r, c)
end

"""
    HomalgRowVector(L, c, R)

Construct a (1 x c)-matrix over the ring R with entries in the list L

```jldoctest
julia> mat = HomalgRowVector(1:5, 5, ZZ)
[1   2   3   4   5]
```
"""
function HomalgRowVector(entries, c, R)::TypeOfMatrixForHomalg
    return HomalgMatrix(entries, 1, c, R)
end

"""
    HomalgRowVector(L, R)

Construct a (1 x c)-matrix over the ring R with entries in the list L, where c = length(L)

```jldoctest
julia> mat = HomalgRowVector(1:5, ZZ)
[1   2   3   4   5]
```
"""
function HomalgRowVector(entries, R)::TypeOfMatrixForHomalg
    return HomalgRowVector(entries, length(entries), R)
end

"""
    HomalgColumnVector(L, r, R)

Construct a (r x 1)-matrix over the ring R with entries in the list L

```jldoctest
julia> mat = HomalgColumnVector(1:5, 5, ZZ)
[1]
[2]
[3]
[4]
[5]
```
"""
function HomalgColumnVector(entries, r, R)::TypeOfMatrixForHomalg
    return HomalgMatrix(entries, r, 1, R)
end

"""
    HomalgColumnVector(L, r, R)

Construct a (r x 1)-matrix over the ring R with entries in the list L, where r = length(L)

```jldoctest
julia> mat = HomalgColumnVector(1:5, ZZ)
[1]
[2]
[3]
[4]
[5]
```
"""
function HomalgColumnVector(entries, R)::TypeOfMatrixForHomalg
    return HomalgColumnVector(entries, length(entries), R)
end

"""
    HomalgDiagonalMatrix(L, R)

Construct a diagonal matrix over the ring R using the list L of diagonal entries

```jldoctest
julia> mat = HomalgDiagonalMatrix(1:5, ZZ)
[1   0   0   0   0]
[0   2   0   0   0]
[0   0   3   0   0]
[0   0   0   4   0]
[0   0   0   0   5]
```
"""
function HomalgDiagonalMatrix(diagonal_entries, R)::TypeOfMatrixForHomalg
    return AbstractAlgebra.block_diagonal_matrix(map(a->HomalgMatrix([a],1,1,R), diagonal_entries))
end

"""
    R * mat

Rewrite the matrix mat over the ring R (if possible)

```jldoctest
julia> mat = HomalgMatrix([1,2,3,4,5,6], 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> QQ * mat
[1//1   2//1   3//1]
[4//1   5//1   6//1]

julia> qmat = QQ * mat
[1//1   2//1   3//1]
[4//1   5//1   6//1]

julia> ZZ * qmat == mat
true
```
"""
Base.:*(R, mat) = AbstractAlgebra.change_base_ring(R, mat)

export HomalgMatrix, HomalgIdentityMatrix, HomalgZeroMatrix, HomalgRowVector, HomalgColumnVector, HomalgDiagonalMatrix

## Properties of homalg matrices

"""
    IsOne(mat)

Return true if mat is an identity matrix, otherwise false

```jldoctest
julia> mat = HomalgIdentityMatrix(3, ZZ)
[1   0   0]
[0   1   0]
[0   0   1]

julia> IsOne(mat)
true
```
"""
function IsOne(mat)::Bool
    return AbstractAlgebra.isone(mat)
end

"""
    IsZero(mat)

Return true if mat is a zero matrix, otherwise false

```jldoctest
julia> mat = HomalgZeroMatrix(3, 2, ZZ)
[0   0]
[0   0]
[0   0]

julia> IsZero(mat)
true
```
"""
function IsZero(mat)::Bool
    return AbstractAlgebra.iszero(mat)
end

"""
    IsEmptyMatrix(mat)

Return true if mat does not contain any entry, otherwise false

```jldoctest
julia> mat = HomalgMatrix([], 0, 0, ZZ)
0 by 0 empty matrix

julia> IsEmptyMatrix(mat)
true
```
"""
function IsEmptyMatrix(mat)::Bool
    return isempty(mat)
end

"""
    IsSymmetricMatrix(mat)

Return true if the matrix mat is symmetric with respect to its main diagonal, otherwise false

```jldoctest
julia> mat = HomalgMatrix([1,2,3,2,4,5,3,5,6], 3, 3, ZZ)
[1   2   3]
[2   4   5]
[3   5   6]

julia> IsSymmetricMatrix(mat)
true
```
"""
function IsSymmetricMatrix(mat)::Bool
    return AbstractAlgebra.is_symmetric(mat)
end

export IsOne, IsZero, IsEmptyMatrix, IsSymmetricMatrix

## Attributes of homalg matrices

"""
    HomalgRing(mat)

Return the ring underlying the matrix mat.

```jldoctest
julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> HomalgRing(mat)
Integers

julia> mat = HomalgMatrix(1:6, 2, 3, QQ)
[1//1   2//1   3//1]
[4//1   5//1   6//1]

julia> HomalgRing(mat)
Rationals
```
"""
function HomalgRing(mat)
    return AbstractAlgebra.base_ring(mat)
end

"""
    NumberRows(mat)

The number of rows of the matrix mat

```jldoctest
julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> NumberRows(mat)
2
```
"""
function NumberRows(mat)::Int64
    return AbstractAlgebra.nrows(mat)
end

"""
    NumberColumns(mat)

The number of columns of the matrix mat

```jldoctest
julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> NumberColumns(mat)
3
```
"""
function NumberColumns(mat)::Int64
    return AbstractAlgebra.ncols(mat)
end

"""
    TransposedMatrix(mat)

Return the transposed matrix of mat

```jldoctest
julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> TransposedMatrix(mat)
[1   4]
[2   5]
[3   6]
```
"""
function TransposedMatrix(mat)::TypeOfMatrixForHomalg
    return Base.transpose(mat)
end

"""
    ConvertMatrixToRow(mat)

Unfold the matrix M row-wise into a row.

```jldoctest
julia> mat = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> ConvertMatrixToRow(mat)
[2   3   4   5   6   7]

```
"""
function ConvertMatrixToRow(mat)::TypeOfMatrixForHomalg
    union_rows = HomalgZeroMatrix(1, 0, HomalgRing(mat))
    for i = 1:(NumberRows(mat))
        row = mat[i:i, :]
        union_rows = UnionOfColumns(HomalgRing(mat),1, [union_rows, row])
    end
    return union_rows
end

"""
    ConvertMatrixToColumn(mat)

Unfold the matrix M column-wise into a column.

```jldoctest
julia> mat = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> ConvertMatrixToColumn(mat)
[2]
[4]
[6]
[3]
[5]
[7]

```
"""
function ConvertMatrixToColumn(mat)::TypeOfMatrixForHomalg
    union_columns = HomalgZeroMatrix(0, 1, HomalgRing(mat))
    #Problem: NumberColumns liefert falschen Datentyp
    for i = 1:(NumberColumns(mat))
        column = mat[:, i:i]
        union_columns = UnionOfRows(HomalgRing(mat),1, [union_columns, column])
    end
    return union_columns
end

"""
    RowReducedEchelonForm(mat)

Return the reduced row-echelon form of the matrix mat.

```jldoctest
julia> mat = HomalgMatrix(reverse(1:9), 3, 3, ZZ)
[9   8   7]
[6   5   4]
[3   2   1]

julia> RowReducedEchelonForm(mat)
([3 0 -3; 0 1 2; 0 0 0], 2)

julia> mat = HomalgMatrix(reverse(1:9), 3, 3, QQ)
[9//1   8//1   7//1]
[6//1   5//1   4//1]
[3//1   2//1   1//1]

julia> RowReducedEchelonForm(mat)
([1 0 -1; 0 1 2; 0 0 0], 2)
"""
function RowReducedEchelonForm(mat::AbstractAlgebra.Generic.MatSpaceElem{BigInt})::Tuple{TypeOfMatrixForHomalg, Int64}
    hnf = AbstractAlgebra.hnf(mat)
    rank = AbstractAlgebra.rank(hnf)
    return hnf, rank
end

function RowReducedEchelonForm(mat::AbstractAlgebra.Generic.MatSpaceElem{Rational{BigInt}})::Tuple{TypeOfMatrixForHomalg, Int64}
    rank, rref = AbstractAlgebra.rref(mat)
    return rref, rank
end

"""
    BasisOfRows(mat)

Return the triangular form of mat

```jldoctest
julia> mat = HomalgMatrix(1:9, 3, 3, ZZ)
[1   2   3]
[4   5   6]
[7   8   9]

julia> BasisOfRows(mat)
[1   2   3]
[0   3   6]

julia> mat = HomalgMatrix(1:9, 3, 3, QQ)
[1//1   2//1   3//1]
[4//1   5//1   6//1]
[7//1   8//1   9//1]

julia> BasisOfRows(mat)
[1//1   0//1   -1//1]
[0//1   1//1    2//1]
```
"""
function BasisOfRows(mat)::TypeOfMatrixForHomalg
    rref, rank = RowReducedEchelonForm(mat)
    return rref[1:rank, :]
end

"""
    BasisOfColumns(mat)

Return the triangular form of mat

```jldoctest
julia> mat = HomalgMatrix(1:9, 3, 3, ZZ)
[1   2   3]
[4   5   6]
[7   8   9]

julia> BasisOfColumns(mat)
[1   0]
[1   3]
[1   6]

julia> mat = HomalgMatrix(1:9, 3, 3, QQ)
[1//1   2//1   3//1]
[4//1   5//1   6//1]
[7//1   8//1   9//1]

julia> BasisOfColumns(mat)
[ 1//1   0//1]
[ 0//1   1//1]
[-1//1   2//1]
```
"""
function BasisOfColumns(mat)::TypeOfMatrixForHomalg
    return TransposedMatrix(BasisOfRows(TransposedMatrix(mat)))
end

export HomalgRing, NumberRows, NumberColumns, TransposedMatrix, ConvertMatrixToRow, ConvertMatrixToColumn, RowReducedEchelonForm, BasisOfRows, BasisOfColumns

## Operations of homalg matrices

"""
    UnionOfRows(R, nr_cols, list)

Return the matrices in list stacked, where all of them have same number of columns nr_cols.

```jldoctest
julia> UnionOfRows(ZZ, 3, [])
0 by 3 empty matrix

julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> UnionOfRows(ZZ, 3, [mat, mat])
[1   2   3]
[4   5   6]
[1   2   3]
[4   5   6]
```
"""
function UnionOfRows(R, nr_cols, list)::TypeOfMatrixForHomalg
    if length(list) == 0
        return HomalgZeroMatrix(0, nr_cols, R)
    end

    return vcat(list...)
end

"""
    UnionOfColumns(R, nr_rows, list)

Return the matrices in list augmented, where all of them have same number of rows nr_rows.
.

```jldoctest
julia> UnionOfColumns(ZZ, 2, [])
2 by 0 empty matrix

julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> UnionOfColumns(ZZ, 2, [mat, mat])
[1   2   3   1   2   3]
[4   5   6   4   5   6]
```
"""
function UnionOfColumns(R, nr_rows, list)::TypeOfMatrixForHomalg
    if length(list) == 0
        return HomalgZeroMatrix(nr_rows, 0, R)
    end

    return hcat(list...)
end

"""
    CertainColumns(mat, list)

Return the matrix of which the i-th column is the k-th column of the homalg matrix M, where k=list[i].

```jldoctest
julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> CertainColumns(mat, [2, 2, 1])
[2   2   1]
[5   5   4]

julia> CertainColumns(mat, [])
2 by 0 empty matrix

```
"""
function CertainColumns(mat, list)::TypeOfMatrixForHomalg
    nr_rows = NumberRows(mat)
    union_columns = HomalgZeroMatrix(nr_rows, 0, HomalgRing(mat))
    for i = 1:(length(list))
        column = mat[:, list[i]:list[i]]
        union_columns = UnionOfColumns(HomalgRing(mat), nr_rows,[union_columns, column])
    end
    return union_columns
end

"""
    CertainRows(mat, list)

Return the matrix of which the i-th row is the k-th row of the homalg matrix M, where k=list[i].

```jldoctest
julia> mat = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> CertainRows(mat, [2, 2, 1])
[4   5]
[4   5]
[2   3]

julia> CertainRows(mat, [])
0 by 2 empty matrix

```
"""
function CertainRows(mat, list)::TypeOfMatrixForHomalg
    nr_cols = NumberColumns(mat)
    union_rows = HomalgZeroMatrix(0, nr_cols, HomalgRing(mat))
    for i = 1:(length(list))
        row = mat[list[i]:list[i], :]
        union_rows = UnionOfRows(HomalgRing(mat), nr_cols,[union_rows, row])
    end
    return union_rows
end

"""
    KroneckerMat(mat1, mat2)

Return the Kronecker (or tensor) product of the two homalg matrices mat1 and mat2.

```jldoctest
julia> mat1 = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> mat2 = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> KroneckerMat(mat1, mat2)
[ 2    3    4    6    6    9]
[ 4    5    8   10   12   15]
[ 6    7   12   14   18   21]
[ 8   12   10   15   12   18]
[16   20   20   25   24   30]
[24   28   30   35   36   42]

julia> KroneckerMat(mat2, mat1)
[ 2    4    6    3    6    9]
[ 8   10   12   12   15   18]
[ 4    8   12    5   10   15]
[16   20   24   20   25   30]
[ 6   12   18    7   14   21]
[24   30   36   28   35   42]

```
"""
function KroneckerMat(mat1, mat2)::TypeOfMatrixForHomalg
    return AbstractAlgebra.kronecker_product(mat1, mat2)
end

"""
    RightDivide(mat2, mat1)

Returns: a homalg matrix or fail

Let mat2 and mat1 be matrices having the same number of columns and defined over the same ring.
The matrix RightDivide(mat2, mat1) is a particular solution of the inhomogeneous (one sided) linear system of equations X(mat1)=(mat2) in case it is solvable.
Otherwise fail is returned. The name RightDivide suggests "X=(mat2)(mat1^(-1))".

```jldoctest
julia> mat1 = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> mat2 = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> RightDivide(mat2, mat2)
[ 1   0   0]
[ 0   1   0]
[-1   2   0]

julia> RightDivide(mat2, mat1)
"fail"
```
"""
function RightDivide(mat2, mat1)::Union{TypeOfMatrixForHomalg, String}
    try
        return AbstractAlgebra.solve_left(mat1, mat2)
    catch y
        return "fail"
    end
end

"""
    LeftDivide(mat1, mat2)

Returns: a homalg matrix or fail

Let mat1 and mat2 be matrices having the same number of rows and defined over the same ring.
The matrix LeftDivide(mat1, mat2) is a particular solution of the inhomogeneous (one sided) linear system of equations (mat1)X=(mat2) in case it is solvable.
Otherwise fail is returned. The name LeftDivide suggests "X=((mat1)^-1)(mat2)".

```jldoctest
julia> mat1 = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> mat2 = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> mat3 = HomalgMatrix([1,0,0,0,0,0], 3, 2, ZZ)
[1   0]
[0   0]
[0   0]

julia> LeftDivide(mat2, mat2)
[1   0]
[0   1]

julia> LeftDivide(mat1, mat2)
[0   -1]
[1    2]

julia> LeftDivide(mat3, mat2)
"fail"
```
"""
function LeftDivide(mat1, mat2)::Union{TypeOfMatrixForHomalg, String}
    try
        return AbstractAlgebra.solve(mat1, mat2)
    catch y
        return "fail"
    end
end

"""
    DecideZeroRows(mat1, mat2)

```jldoctest
julia> mat1 = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> mat2 = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> reduced_mat1 = DecideZeroRows(mat2, mat1)
[0   1]
[0   1]
[0   1]

julia> RightDivide(mat2, mat1)
"fail"

julia> RightDivide(mat2 - reduced_mat1, mat1)
[-1   1   0]
[-2   2   0]
[-3   3   0]

julia> DecideZeroRows(mat1, mat1)
[0   0]
[0   0]
[0   0]

julia> mat3 = HomalgMatrix([4, 6, 2, 2], 2, 2, ZZ)
[4   6]
[2   2]

julia> DecideZeroRows(mat3, mat1)
[0   0]
[0   0]

julia> DecideZeroRows(mat1, mat3)
[1   0]
[1   0]
[1   0]
```
"""
function DecideZeroRows(B,A)
    #A,B are defined over the same ring
    ring = HomalgRing(B)

    nr_rows_A = NumberRows(A)
    nr_rows_B = NumberRows(B)

    #A,B have the same number of columns
    nr_cols = NumberColumns(B)

    null_mat_a = HomalgZeroMatrix(nr_rows_A, nr_rows_B, ring)
    ident_mat_b = HomalgIdentityMatrix(nr_rows_B, ring)

    list_of_matrices = [UnionOfRows(ring, nr_rows_B, [ident_mat_b, null_mat_a]), UnionOfRows(ring, nr_cols,[B, A])]
    temp_mat = UnionOfColumns(ring, nr_rows_B + nr_rows_A,list_of_matrices)
    resulting_mat = BasisOfRows(temp_mat)

    return resulting_mat[1:nr_rows_B, nr_rows_B+1:nr_rows_B+nr_cols]
end

"""
    DecideZeroRows(mat1, mat2)

```jldoctest
julia> mat1 = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> mat3 = HomalgMatrix([3, 1, 7, 1, 11, 1], 3, 2, ZZ)
[ 3   1]
[ 7   1]
[11   1]

julia> DecideZeroColumns(mat3, mat1)
[0   0]
[0   0]
[0   0]
```
"""
function DecideZeroColumns(B, A)
    return TransposedMatrix(DecideZeroRows(TransposedMatrix(B), TransposedMatrix(A)))
end

export UnionOfRows, UnionOfColumns, KroneckerMat, CertainColumns, CertainRows, RightDivide, LeftDivide, DecideZeroRows, DecideZeroColumns

end
