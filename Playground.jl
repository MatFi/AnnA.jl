using AnnABase , Unitful #,AnnAPlot, Plots

parameters = AnnABase.Parameters(N=300,light = t->0,Rₛₕ =Inf*1u"V/A*m^2" )
prob = AnnABase.IVProblem(parameters,[1.2,0.3]u"V",0.005u"V/s",double_sweep=true)
sol = AnnABase.solve(prob)
using BenchmarkTools


using ProfileView
@profview AnnABase.solve(prob)
