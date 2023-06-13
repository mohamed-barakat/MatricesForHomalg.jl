@testset "attributes" begin
    mat = HomalgMatrix(1:6, 2, 3, ZZ)

    @test NumberRows(mat) == 2
    @test NumberColumns(mat) == 3

end