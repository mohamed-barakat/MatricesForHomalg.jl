@testset "operations" begin
    mat = HomalgMatrix(1:6, 2, 3, ZZ)
    
    tensor = KroneckerMat(mat, mat)
    @test NumberRows(tensor) == 4
    @test NumberColumns(tensor) == 9

end
