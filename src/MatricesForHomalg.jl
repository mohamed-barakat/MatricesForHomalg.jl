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
"""
function HomalgMatrix(entries, r, c, R)::TypeOfMatrixForHomalg
    return matrix(R, r, c, entries)
end

"""
    HomalgIdentityMatrix(r, R)

Construct a (r x r)-identity matrix over the ring R
"""
function HomalgIdentityMatrix(r, R)::TypeOfMatrixForHomalg
    return AbstractAlgebra.identity_matrix(R, r)
end

"""
    HomalgZeroMatrix(r, c, R)

Construct a (r x c)-zero matrix over the ring R
"""
function HomalgZeroMatrix(r, c, R)::TypeOfMatrixForHomalg
    return AbstractAlgebra.zero_matrix(R, r, c)
end

"""
    HomalgRowVector(L, c, R)

Construct a (1 x c)-matrix over the ring R with entries in the list L
"""
function HomalgRowVector(entries, c, R)::TypeOfMatrixForHomalg
    return HomalgMatrix(entries, 1, c, R)
end

"""
    HomalgRowVector(L, R)

Construct a (1 x c)-matrix over the ring R with entries in the list L, where c = length(L)
"""
function HomalgRowVector(entries, R)::TypeOfMatrixForHomalg
    return HomalgRowVector(entries, length(entries), R)
end

"""
    HomalgColumnVector(L, r, R)

Construct a (r x 1)-matrix over the ring R with entries in the list L
"""
function HomalgColumnVector(entries, r, R)::TypeOfMatrixForHomalg
    return HomalgMatrix(entries, r, 1, R)
end

"""
    HomalgColumnVector(L, r, R)

Construct a (r x 1)-matrix over the ring R with entries in the list L, where r = length(L)
"""
function HomalgColumnVector(entries, R)::TypeOfMatrixForHomalg
    return HomalgColumnVector(entries, length(entries), R)
end

"""
    HomalgDiagonalMatrix(L, R)

Construct a diagonal matrix over the ring R using the list L of diagonal entries
"""
function HomalgDiagonalMatrix(diagonal_entries, R)::TypeOfMatrixForHomalg
    return AbstractAlgebra.block_diagonal_matrix(map(a->HomalgMatrix([a],1,1,R), diagonal_entries))
end

"""
    R * mat

Rewrite the matrix mat over the ring R (if possible)
"""
Base.:*(R, mat) = AbstractAlgebra.change_base_ring(R, mat)

export HomalgMatrix, HomalgIdentityMatrix, HomalgZeroMatrix, HomalgRowVector, HomalgColumnVector, HomalgDiagonalMatrix

## Properties of homalg matrices

"""
    IsOne(mat)

Return true if mat is an identity matrix, otherwise false
"""
function IsOne(mat)::Bool
    return AbstractAlgebra.isone(mat)
end

"""
    IsZero(mat)

Return true if mat is a zero matrix, otherwise false
"""
function IsZero(mat)::Bool
    return AbstractAlgebra.iszero(mat)
end

export IsOne, IsZero

## Attributes of homalg matrices

"""
    NumberRows(mat)

The number of rows of the matrix mat
"""
function NumberRows(mat)::BigInt
    return AbstractAlgebra.nrows(mat)
end

"""
    NumberColumns(mat)

The number of columns of the matrix mat
"""
function NumberColumns(mat)::BigInt
    return AbstractAlgebra.ncols(mat)
end

export NumberRows, NumberColumns

end
