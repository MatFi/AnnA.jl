@testset "Problems" begin
    parameters= AnnA.Parameters(light = t->0,Rₛₕ =Inf*1u"V/A*m^2" )
    prob = AnnA.IVProblem(parameters,[-0.2,1.4]u"V",0.005u"V/s")
    @test prob isa AnnA.IVProblem

    #
    sol =  AnnA.solve(prob)
    @test sol isa AnnA.ProblemSolution
    @test isapprox(sol.df.V[1] - sol.df.V[end] .|> ustrip, 0; atol=1e-4)
    
    #forward direction
    prob = AnnA.IVProblem(parameters,[-0.2,1.4]u"V",0.005u"V/s",double_sweep=false)
    sol =  AnnA.solve(prob) 
    @test sol isa AnnA.ProblemSolution
    @test !isapprox(sol.df.V[1] - sol.df.V[end] .|> ustrip, 0; atol=1e-4)
    #reverse direction
    prob = AnnA.IVProblem(parameters,[1.2,0.4]u"V",0.005u"V/s",double_sweep=false)
    sol = AnnA.solve(prob) 
    @test sol isa AnnA.ProblemSolution
    @test !isapprox(sol.df.V[1] - sol.df.V[end] .|> ustrip, 0; atol=1e-4)

    #test the ocvd  
    prob_ocvd = AnnA.OCVDProblem(parameters,10u"s",1e3u"s")
    @test  AnnA.solve(prob_ocvd) isa AnnA.OCVDSolution
    
end