@testset "Cell" begin
    @testset "Zeros in Recombination" begin
        prm = AnnA.Parameters(vₚₑ=0u"m/s")
        @test_nowarn AnnA.Cell(prm)

        prm = AnnA.Parameters(vₚₕ=0u"m/s")
        @test_nowarn AnnA.Cell(prm)

        prm = AnnA.Parameters(τₚ=0u"m/s")
        @test_nowarn AnnA.Cell(prm)

        #give a initialisation warning
        prm = AnnA.Parameters(V= t -> 1)
        @test_logs (:warn,) AnnA.Cell(prm,mode=:oc)
    end
end
