function solve(c::Cell;tstart=c.alg_ctl.tstart, tend=c.alg_ctl.tend,kwargs...)
    if (c.initialized==:false)
        initial_conditions!(c)
    #    if c.initialized != :Sucess
        #    @warn "Initialisation Failure"
    #    end
    end

    prob = get_problem(c,tstart=tstart,tend=tend)

    #tol_callback = AutoAbstol(false;init_curmax=c.u0 .+ 000.01)
    #c.alg_ctl.abstol = ones(length(c.u0))*c.alg_ctl.abstol
#    (sav_callback, sav_sol) = SavingCallback(c)
#    callbacks = CallbackSet(sav_callback, )
    @debug "Solve"
    τᵢ = c.parameters.τᵢ
    c.alg_ctl.dt = 1e-9*ustrip(τᵢ  |> u"s")
    c.alg_ctl.dtmin = ustrip(1e-40*τᵢ |> u"s")
    sol = solve(prob,c.alg_ctl.alg; c.alg_ctl.kwargs...,kwargs... )
    c.sol = sol
    return sol
end
function get_problem(c::Cell;tstart=0.0u"s",tend = 20.0u"s")

     ntype = c.alg_ctl.numtype 
     tspan = (convert(ntype,upreferred(tstart/c.parameters.τᵢ)), convert(ntype,upreferred(tend/c.parameters.τᵢ)))
    if !(ntype <: Float64)
        u0 = ntype.(c.u0)
        Jac = ntype.(c.Jac)
        tspan = ntype.(tspan)
        M= c.M
    else
        u0 = c.u0
        Jac =c.Jac
        M= c.M
    end
    odefun = ODEFunction(c.rhs;
        mass_matrix = M,
        jac_prototype = Jac,
        colorvec = matrix_colors(c.Jac),
    )
   

    #tol_callback = AutoAbstol(false;init_curmax=u0 .+ 0.1)

    prob = ODEProblem(odefun,u0,tspan,c)
end
