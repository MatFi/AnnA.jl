function solve(c::Cell)
    if !c.initialized
        initial_conditions!(c)
        c.initialized=true
    end
    u0 = c.u0

    colorvec = matrix_colors(c.Jac)

    odefun = ODEFunction(c.rhs;
        mass_matrix = c.M,
        jac_prototype = c.Jac,
        colorvec = matrix_colors(c.Jac),
    )

    tspan = (0.0,5000*π/ustrip(c.parameters.τᵢ))
    if c.mode==:oc
        tspan = (0.0,1e6)
    end


    tol_callback = AutoAbstol(false;init_curmax=u0 .+ 0.1)
    (sav_callback, sav_sol) = SavingCallback(c)
    callbacks = CallbackSet(sav_callback, tol_callback )
    prob = ODEProblem(odefun,u0,tspan,c)
    solve(prob,c.alg_ctl.alg;
        progress_steps = 50,
        progress = true,
        callback = callbacks,
        dt = 1e-20,
        dtmin = 1e-20,
        force_dtmin = true,
        reltol = c.alg_ctl.reltol,
        abstol = u0 .* 0
    )
    return sav_sol
end
