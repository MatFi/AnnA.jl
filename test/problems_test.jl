@testset "Problems" begin
    parameters= Parameters(light = t->0,Rₛₕ =Inf*1u"V/A*m^2" )
    prob = IVProblem(parameters,[-0.2,1.4]u"V",0.005u"V/s")
    @test prob isa AnnA.IVProblem

    #
    sol =  AnnA.solve(prob)
    @test sol isa AnnA.ProblemSolution
    @test isapprox(sol.df.V[1] - sol.df.V[end] .|> ustrip, 0; atol=1e-4)
    
    #forward direction
    prob = AnnA.IVProblem(parameters,[-0.2,1.4]u"V",0.005u"V/s",double_sweep=false)
    sol =  AnnA.solve(prob) 
    @test sol isa AnnA.ProblemSolution{AnnA.IV}
    @test !isapprox(sol.df.V[1] - sol.df.V[end] .|> ustrip, 0; atol=1e-4)
    #reverse direction
    prob = AnnA.IVProblem(parameters,[1.2,0.4]u"V",0.005u"V/s",double_sweep=false)
    sol = AnnA.solve(prob) 
    @test sol isa AnnA.ProblemSolution{AnnA.IV}
    @test !isapprox(sol.df.V[1] - sol.df.V[end] .|> ustrip, 0; atol=1e-4)

    #test the ocvd  
    @test OCVDProblem(parameters,10u"s",1e3u"s",alg_control=AlgControl()) isa OCVDProblem
    prob_ocvd = AnnA.OCVDProblem(parameters,10u"s",1e3u"s")
    sol = AnnA.solve(prob_ocvd)
    @test sol isa AnnA.ProblemSolution{AnnA.OCVD}

    props = [:t_on, :t_decay, :V_on, :V_decay, :sol_on, :sol_decay]
    for p in props
        @test !isnothing(getproperty(sol, p)) 
    end
    @test sol.t_on[end] ≈ 10u"s"
    @test sol.t_decay[end] ≈ 1000u"s"
    @test sol.V_on[end] > 1.2u"V"
    @test sol.V_decay[end] ≈ 0.1597653552u"V"


end