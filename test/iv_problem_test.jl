@testset "Problems" begin

    parameters= AnnABase.Parameters(light = t->0,Rₛₕ =Inf*1u"V/A*m^2" )
    prob    = AnnABase.IVProblem(parameters,[-0.2,1.4]u"V",0.005u"V/s")
    @test prob isa AnnABase.IVProblem
    @test  AnnABase.solve!(prob)==nothing
end
