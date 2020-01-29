@testset "Problems" begin

    parameters= AnnABase.Parameters(light = t->0,Rₛₕ =Inf*1u"V/A*m^2" )
    prob    = AnnABase.IVProblem(parameters,[-0.2,1.4]u"V",0.005u"V/s")
    sols = AnnABase.solve!(prob)

end


plot(sol.prob.p.sol.V .|>ustrip,sol.prob.p.sol.J .|>ustrip .|>abs, yscale=:log10,ylims=(1e-8,1e3))
sols = AnnABase.solve!(prob)

prob.sol
