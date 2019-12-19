
"""
    initial_conditions(c::Cell)

Calculates the initial condition of a given `Cell` using NLsolve.
Returns the steady state solution vector to the conditions of
time `t = 0`.
"""
function initial_conditions(c::Cell)

    @info "initialisation : nlsolve on first guess in :precondition mode"
    u0 = c.u0
    # NLSolve will do the rest in no time! Actually this onnly works because
    # the algebraic varibles nowher apper as timederivative in the massmatrix
    if c.mode == :oc
        p_init = setproperties(c.parameters, V = t -> 0) #ensure that we initialize to
    else
        p_init = c.parameters
    end

    c_init = setproperties(
        c,
        parameters = p_init,
        mode = :precondition,
        Jac = get_jac_sparse_pattern(c.g; mode = :precondition),
    )
    # get the sparse colored jacobian for fast NLsolve
    function j!(jac, x, c_init)
        colors = matrix_colors(c_init.Jac)
        forwarddiff_color_jacobian!(
            jac,
            (dx, x) -> c.rhs(dx, x, c_init, 0),
            x;
            colorvec = colors,
            sparsity = c_init.Jac,
        )
        return nothing
    end

    df = OnceDifferentiable(
        (dx, x) -> c.rhs(dx, x, c_init, 0),
        (jac, x) -> j!(jac, x, c_init),
        u0,
        copy(u0),
        c_init.Jac,
    )
    u1 = nlsolve(
        df,
        u0;
        iterations = 100,
        ftol = c.alg_ctl.ss_tol,
        show_trace = haskey(ENV, "JULIA_DEBUG"),
    )
    @debug "NLsolve: " u1.f_converged u1.iterations u1.residual_norm

    # in :oc mode a second init step is needed (in case we have light)
    if c.mode == :oc
        @info "initalisatiion: stating conditions in :oc mode"
        c_init = c
        df = OnceDifferentiable(
            (dx, x) -> c.rhs(dx, x, c_init, 0),
            (jac, x) -> j!(jac, x, c_init),
            u0,
            copy(u0),
            c_init.Jac,
        )
        u1 = nlsolve(
            df,
            u1.zero;
            iterations = 100,
            ftol = c.alg_ctl.ss_tol,
            method = :trust_region,
            factor = u1.ftol,
            show_trace = haskey(ENV, "JULIA_DEBUG"),
        )
        @debug "NLsolve: " u1.f_converged u1.iterations u1.residual_norm
    end
    @info "Initialisation tolerance of $(u1.residual_norm) reached"
#=

    odefun = ODEFunction((dx,x,p,t) -> c.rhs(dx, x, p, 0);
        mass_matrix = c.M,
        jac_prototype = c.Jac,
        colorvec = matrix_colors(c.Jac),
    )

    prob = SteadyStateProblem(
        odefun,
        u1.zero,
        c;
    )

    sol = solve(prob,DynamicSS(Rodas5();abstol=1e-8,reltol=1e-6,tspan=1e3),progress = true, progress_steps=50)

    return sol.u
=#
    return u1.zero
end

"""
    initial_conditions!(c::Cell)

Inplace mutation variant of `initial_conditions(c::Cell)``
"""
function initial_conditions!(c::Cell)
    c.u0 = initial_conditions(c::Cell)
end

"""
    init_guess(g::Grid,ndim::NodimParameters)

Returns an initial guess for NLsolve root finding. The `NodimParameter` value
is needed to guarantee consitancy of the guess.
"""
function init_guess(g::Grid, ndim::NodimParameters)
    P_init = ones(size(g.x))
    phi_init = zeros(size(g.x))
    p_init = ones(size(g.x)) * ndim.kₕ .* range(0, 1, length = g.N + 1) * 1e-1
    n_init = ones(size(g.x)) * ndim.kₑ .* range(1, 0, length = g.N + 1) * 1e-1
    phiE_init = zeros(size(g.xₑ))
    nE_init = ones(size(g.xₑ)) .* range(1, n_init[1] / ndim.kₑ, length = g.Nₑ)
    phiH_init = zeros(size(g.xₕ))
    pH_init = ones(size(g.xₕ)) .* range(p_init[end] / ndim.kₕ, 1, length = g.Nₕ)
    return [
        P_init
        phi_init
        n_init
        p_init
        phiE_init
        nE_init
        phiH_init
        pH_init
    ]
end
