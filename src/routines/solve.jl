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

    sol = solve(prob,c.alg_ctl.alg; c.alg_ctl.kwargs...,kwargs... )
    c.sol = sol
    return sol
end
function get_problem(c::Cell;tstart=0.0u"s",tend = 20.0u"s")
    u0 = c.u0

    colorvec = matrix_colors(c.Jac)

    odefun = ODEFunction(c.rhs;
        mass_matrix = c.M,
        jac_prototype = c.Jac,
        colorvec = matrix_colors(c.Jac),
    )
    tspan = (convert(Float64,upreferred(tstart/c.parameters.τᵢ)), convert(Float64,upreferred(tend/c.parameters.τᵢ)))

    #tol_callback = AutoAbstol(false;init_curmax=u0 .+ 0.1)

    prob = ODEProblem(odefun,u0,tspan,c)
end
