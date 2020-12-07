@testset "Problems" begin
    parameters= AnnABase.Parameters(light = t->0,Rₛₕ =Inf*1u"V/A*m^2" )
    prob = AnnABase.IVProblem(parameters,[-0.2,1.4]u"V",0.005u"V/s")
    @test prob isa AnnABase.IVProblem
    @test  AnnABase.solve(prob) isa AnnABase.ProblemSolution

    #forward direction
    prob = AnnABase.IVProblem(parameters,[-0.2,1.4]u"V",0.005u"V/s",double_sweep=false)
    @test  AnnABase.solve(prob) isa AnnABase.ProblemSolution

    #reverse direction
    prob = AnnABase.IVProblem(parameters,[1.2,0.4]u"V",0.005u"V/s",double_sweep=false)
    @test  AnnABase.solve(prob) isa AnnABase.ProblemSolution

    #test the ocvd
    prob_ocvd = AnnABase.OCVDProblem(parameters,10u"s",1e3u"s")
    @test  AnnABase.solve(prob_ocvd) isa AnnABase.OCVDSolution
end
