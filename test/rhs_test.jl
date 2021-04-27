@testset "Rhs" begin
    prm = AnnA.Parameters(vₚₑ = 0u"m/s")
    cell = AnnA.Cell(prm)
    u = rand(length(cell.g))
    du = similar(u)
    @test_nowarn cell.rhs(du,u,cell,0)

    @test_throws ErrorException AnnA.Cell(prm, mode=:err)

    #alls = @benchmark AnnA.rhs!(du,u,cell,1.0)
    #@test allocs(alls) == 0


    u=zeros(10)
    bu = BigFloat.(u)
    lc = AnnA.lcache()
    @test AnnA.get_tmp(lc,u) isa Vector{Float64}
    @test AnnA.get_tmp(lc,bu) isa Vector{BigFloat}
    @test AnnA.get_tmp(lc,bu) isa Vector{BigFloat}
    @test AnnA.get_tmp(lc,u) isa Vector{Float64}
end
