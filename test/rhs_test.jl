@testset "Rhs" begin
    prm = AnnABase.Parameters(vₚₑ = 0u"m/s")
    cell = AnnABase.Cell(prm)
    u = rand(length(cell.g))
    du = similar(u)
    @test_nowarn cell.rhs(du,u,cell,0)

    @test_throws ErrorException AnnABase.Cell(prm, mode=:err)

    #alls = @benchmark AnnABase.rhs!(du,u,cell,1.0)
    #@test allocs(alls) == 0
end
