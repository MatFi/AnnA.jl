function solve(c::Cell;tend=20)
    if (c.initialized==:false)
        initial_conditions!(c)
    #    if c.initialized != :Sucess
        #    @warn "Initialisation Failure"
    #    end
    end

    prob = get_problem(c,tend=tend)

    tol_callback = AutoAbstol(false;init_curmax=c.u0 .+ 000.00001)
    (sav_callback, sav_sol) = SavingCallback(c)
    callbacks = CallbackSet(sav_callback, tol_callback )
    @info "Solve"
    #TODO: abstol needs to bedone more nice
    sol = solve(prob,c.alg_ctl.alg;
        progress_steps = 50,
        progress = true,
        callback = sav_callback,#callbacks,
        dt =  1e-15,
        dtmin = 1e-20,# ustrip(1e-20*c.parameters.τᵢ), #1e-20,
        force_dtmin = true,
        reltol = c.alg_ctl.reltol,
        abstol = 1e-12,#c.u0 .* 0, #1e-12,#c.u0 .* 0,
        maxiters= 5000
    )
    c.sol = sav_sol
    return sol
end
function get_problem(c::Cell;tend = 20.0)
    u0 = c.u0

    colorvec = matrix_colors(c.Jac)

    odefun = ODEFunction(c.rhs;
        mass_matrix = c.M,
        jac_prototype = c.Jac,
        colorvec = matrix_colors(c.Jac),
    )

    tspan = (0.0,ustrip(tend/c.parameters.τᵢ))



    tol_callback = AutoAbstol(false;init_curmax=u0 .+ 0.1)
    (sav_callback, sav_sol) = SavingCallback(c)
    callbacks = CallbackSet(sav_callback, tol_callback )
    prob = ODEProblem(odefun,u0,tspan,c)
end
