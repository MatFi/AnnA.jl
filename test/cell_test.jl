@testset "Cell" begin
    @testset "Zeros in Recombination" begin
        prm = AnnABase.Parameters(vₚₑ=0u"m/s")
        @test_nowarn AnnABase.Cell(prm)

        prm = AnnABase.Parameters(vₚₕ=0u"m/s")
        @test_nowarn AnnABase.Cell(prm)

        prm = AnnABase.Parameters(τₚ=0u"m/s")
        @test_nowarn AnnABase.Cell(prm)

        #give a initialisation warning
        prm = AnnABase.Parameters(V= t -> 1)
        @test_logs (:warn,) AnnABase.Cell(prm,mode=:oc)
    end
end
