@testset "Rhs" begin
    prm = AnnABase.Parameters(vₚₑ = 0)
    cell = AnnABase.Cell(prm)
    u = rand(length(cell.g))
    du = similar(u)
    @test_nowarn cell.rhs(du,u,cell,0)

    @test_throws ErrorException AnnABase.Cell(prm, mode=:err)

    cell_err = setproperties(cell,mode=:err)
    @test_throws ErrorException cell.rhs(du,u,cell_err,0)

    #alls = @benchmark AnnABase.rhs!(du,u,cell,1.0)
    #@test allocs(alls) == 0
end
