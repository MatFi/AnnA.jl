function solve(c::Cell;tend=c.alg_ctl.tend)
    if (c.initialized==:false)
        initial_conditions!(c)
    #    if c.initialized != :Sucess
        #    @warn "Initialisation Failure"
    #    end
    end

    prob = get_problem(c,tend=tend)

    tol_callback = AutoAbstol(false;init_curmax=c.u0 .+ 000.00001)
    (sav_callback, sav_sol) = SavingCallback(c)
    callbacks = CallbackSet(sav_callback, )
    @info "Solve"

    sol = solve(prob,c.alg_ctl.alg; c.alg_ctl.kwargs... )
    c.sol = sav_sol
    return sol
end
function get_problem(c::Cell;tend = 20.0u"s")
    u0 = c.u0

    colorvec = matrix_colors(c.Jac)

    odefun = ODEFunction(c.rhs;
        mass_matrix = c.M,
        jac_prototype = c.Jac,
        colorvec = matrix_colors(c.Jac),
    )
    tspan = (0.0, convert(Float64,(tend/c.parameters.τᵢ)))



    tol_callback = AutoAbstol(false;init_curmax=u0 .+ 0.1)
    (sav_callback, sav_sol) = SavingCallback(c)
    callbacks = CallbackSet(sav_callback, tol_callback )
    prob = ODEProblem(odefun,u0,tspan,c)
end
