using AnnABase , Unitful #,AnnAPlot, Plots

parameters = AnnABase.Parameters(N=300,light = t->0,Rₛₕ =Inf*1u"V/A*m^2")
prob_iv = AnnABase.IVProblem(parameters,[1.2,0.3]u"V",0.005u"V/s",double_sweep=true)
sol = AnnABase.solve(prob_iv)
prob_ocvd = AnnABase.OCVDProblem(parameters,1u"s",1e3u"s")
sol = AnnABase.solve(prob_ocvd)
using BenchmarkTools
prob_ocvd.
parameters isa AnnABase.AbstractParameters
1u"s" isa Unitful.AbstractQuantity
using ProfileView
@profview AnnABase.solve(prob)
