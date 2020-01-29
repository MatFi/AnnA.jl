using AnnABase,AnnAPlot, Plots, Unitful

parameters = AnnABase.Parameters(light = t->0,Rₛₕ =Inf*1u"V/A*m^2" )
prob = AnnABase.IVProblem(parameters,[-0.2,1.4]u"V",0.005u"V/s")
AnnABase.solve!(prob)
typeof(prob.sol[1])

AnnABase.calculate_currents(prob.sol[1])
