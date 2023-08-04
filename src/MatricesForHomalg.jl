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
function BasisOfRows(mat::AbstractAlgebra.Generic.MatSpaceElem{BigInt})::TypeOfMatrixForHomalg
    hnf = AbstractAlgebra.hnf(mat)
    rank = AbstractAlgebra.rank(hnf)
    return hnf[1:rank, :]
end

function BasisOfRows(mat::AbstractAlgebra.Generic.MatSpaceElem{Rational{BigInt}})::TypeOfMatrixForHomalg
    rank, rref = AbstractAlgebra.rref(mat)
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

export HomalgRing, NumberRows, NumberColumns, TransposedMatrix, BasisOfRows, BasisOfColumns

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

export UnionOfRows, UnionOfColumns, KroneckerMat

end
