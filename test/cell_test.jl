@testset "Cell" begin
    @testset "Special Parameter Sets" begin
        prm = AnnA.Parameters(vₚₑ=0u"m/s")
        @test_nowarn AnnA.Cell(prm)

        prm = AnnA.Parameters(vₚₕ=0u"m/s")
        @test_nowarn AnnA.Cell(prm)

        prm = AnnA.Parameters(τₚ=0u"m/s")
        @test_nowarn AnnA.Cell(prm)

        prm = AnnA.Parameters(Rₛₕ=Inf*u"V/A*m^2",Rₛ=Inf*u"V/A*m^2")
        @test_throws ErrorException AnnA.Cell(prm,mode=:cc)
    end
end
