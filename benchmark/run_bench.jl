using DiffEqDevTools
using DiffEqBase
using OrdinaryDiffEq
using AnnA
using Unitful

# create test problem
parm = Parameters()
cell = AnnA.Cell(parm,mode=:oc)
AnnA.initial_conditions!(cell)
prob = AnnA.get_problem(cell,tstart=0u"s",tend=0.1u"s")|>deepcopy

#create test solution
sol = solve(prob,Rodas4P2(); abstol=1e-10,reltol=1e-10)
test_sol = TestSolution(sol)

callback=nothing
reltols=1 ./10 .^(1:2:9)
abstols=1 ./10 .^(1:2:9)
setups = [
        Dict(:alg => Rodas4(),:dtmin => 1e-20, :force_dtmin => false, :progress_steps => 10, :progress => true, :callback => callback),
        Dict(:alg => Rodas4P(),:dtmin => 1e-20, :force_dtmin => false, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg =>Velds4(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg =>Trapezoid(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg =>Rosenbrock23(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg =>GRK4T(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg =>Veldd4(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),

        #Dict(:alg =>RosShamp4(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback), #GOOD
        #Dict(:alg =>ROS34PW2(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
        #Dict(:alg =>ImplicitEuler(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),
  #      Dict(:alg =>ROS34PW3(),:dtmin => 1e-20, :force_dtmin => true, :progress_steps => 10, :progress => true, :callback => callback),

       
        Dict(:alg => Rodas4P2(),:dtmin => 1e-20, :force_dtmin => false, :progress_steps => 10, :progress => true, :callback => callback),
        Dict(:alg => Rodas5(),:dtmin => 1e-20, :force_dtmin => false, :progress_steps => 10, :progress => true, :callback => callback),
        Dict(:alg => ABDF2(),:dtmin => 1e-40, :force_dtmin => false, :progress_steps => 10, :progress => true, :callback => callback),
      #  Dict(:alg => Rodas3(),:dtmin => 1e-20, :force_dtmin => false, :progress_steps => 10, :progress => true, :callback => callback),
        Dict(:alg => ROS3P(),:dtmin => 1e-20, :force_dtmin => false, :progress_steps => 10, :progress => true, :callback => callback),
    ];

wp = WorkPrecisionSet(prob, abstols, reltols, setups; appxsol=test_sol,
                      maxiters=Int(1e4), dt=1e-10,verbose = true,numruns=2)
wp.wps
dwp= deepcopy(wp)
dwp.wps=wp.wps[1:4]
p= plot(wp.wps[8])
plot!.((p,),wp.wps[1:6])
p
display(p)
using Plots

solve(prob,ROS3P(); abstol=1e-5,reltol=1e-5,dmin=1e-40,force_dtmin=true,  maxiters=Int(1e3))

solve(prob,JVODE_BDF(); abstol=1e-5,reltol=1e-5,dmin=1e-40,force_dtmin=true,  maxiters=Int(1e3))

ImplicitHairerWannerExtrapolation