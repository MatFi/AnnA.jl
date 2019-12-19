

using DiffEqDevTools
using OrdinaryDiffEq, DiffEqCallbacks
callback = AutoAbstol(true;init_curmax=u[1].u[1].+0.1)
callback = nothing
test_sol = TestSolution(u[1])
reltols=1 ./10 .^(1:2:7)
abstols=1 ./10 .^(1:2:7)
setups = [
        #Dict(:alg => Rodas4(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg => Rodas4P(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg =>Velds4(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg =>Trapezoid(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg =>Rosenbrock23(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg =>GRK4T(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg =>Veldd4(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),

        #Dict(:alg =>RosShamp4(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback), #GOOD
        #Dict(:alg =>ROS34PW2(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg =>ImplicitEuler(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        Dict(:alg =>ROS34PW3(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),


        Dict(:alg => Rodas5(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
    ];

wp = WorkPrecisionSet(u[2], abstols, reltols, setups; appxsol=test_sol,
                      maxiters=Int(1e3), verbose = true,numruns=2)

plot(wp)


using Plots

u.u[1][4*cell.g.N+4+2*cell.g.Nₑ+cell.g.Nₕ-10]

a= [u.u[i][4*cell.g.N+4+2*cell.g.Nₑ+cell.g.Nₕ-1] for i in 1:length(u.u)]

plotly()
plot(abs.(u.t.-1) .+1e-16 ,a,xscale=:log10)
plot!(sol.t.+1e-16,c.ndim.G.light.(sol.t),xscale=:log10)

x=collect(0:0.05:1)
plot(sol.u[1])
sol =0

a = 0

sol = nothing

a = nothing

using SharedArrays
