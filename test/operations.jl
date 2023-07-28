@testset "operations" begin
    mat = HomalgMatrix(1:6, 2, 3, ZZ)
    
    tensor = KroneckerMat(mat, mat)
    @test NumberRows(tensor) == 4
    @test NumberColumns(tensor) == 9

    empty_rows0 = UnionOfRows(ZZ, 0, [])
    empty_columns0 = UnionOfColumns(ZZ, 0, [])
    @test NumberRows(empty_rows0) == 0
    @test NumberColumns(empty_rows0) == 0
    @test NumberRows(empty_columns0) == 0
    @test NumberColumns(empty_columns0) == 0

    empty_rows = UnionOfRows(ZZ, 3, [])
    empty_columns = UnionOfColumns(ZZ, 2, [])
    @test NumberRows(empty_rows) == 0
    @test NumberColumns(empty_rows) == 3
    @test NumberRows(empty_columns) == 2
    @test NumberColumns(empty_columns) == 0

    u_rows = UnionOfRows(ZZ, 3, [mat, mat])
    u_cols = UnionOfColumns(ZZ, 2, [mat, mat])
    @test NumberRows(u_rows) == 4
    @test NumberColumns(u_rows) == 3
    @test NumberRows(u_cols) == 2
    @test NumberColumns(u_cols) == 6

end
