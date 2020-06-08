
"""
    initial_conditions(c::Cell)

Calculates the initial condition of a given `Cell` using NLsolve.
Returns the steady state solution vector to the conditions of
time `t = 0`.
"""
function initial_conditions(c::Cell)

    #@info "initialisation : nlsolve on first guess in :precondition mode"
    u0 = c.u0
    # NLSolve will do the rest in no time! Actually this onnly works because
    # the algebraic varibles nowher apper as timederivative in the massmatrix
    cc =deepcopy(c)
    #if c.mode == :oc
        p_init = setproperties(cc.parameters, V = t -> 0) #ensure that we initialize to
        c_init= setproperties(cc, parameters = p_init,
            mode = :precondition,
            Jac = get_jac_sparse_pattern(cc.g; mode = :cc))
    #else
#        c_init= deepcopy(c)
#    end
    c_init = Cell(p_init;mode = :cc,alg_ctl = c.alg_ctl)
    #u0 = nl_solve_intiter(c_init,u0;ftol=c.alg_ctl.ss_tol,factor=2).zero
#    return u0
    # in :oc mode a second init step is needed (in case we have light)
    if c.mode == :occ  #legacy
        @info "initalisatiion: stating conditions in :oc mode"
        u0 = nl_solve_intiter(c,u1.zero;ftol=c.alg_ctl.ss_tol,factor=u1.ftol).zero
    end
#    c_init=deepcopy(c)

    @debug "Init_Solve"

    odefun = ODEFunction((dx,x,p,t) -> c_init.rhs(dx, x, p, 0);
        mass_matrix = c_init.M,
        jac_prototype = c_init.Jac,
        colorvec = matrix_colors(c_init.Jac),
    )
    τᵢ = c_init.parameters.τᵢ
    prob = ODEProblem(odefun,u0,(0,1e6*ustrip(τᵢ |> u"s")),c_init)

    ss_cb = TerminateSteadyState(c.alg_ctl.ss_tol,c.alg_ctl.ss_tol)

    #abstol_cb = AutoAbstol(false;init_curmax=u0 .+ 0.001)
    cb = CallbackSet(ss_cb,)#abstol_cb )
    sol = solve(prob,Rodas5();
        progress_steps = 50,
        progress = true,
        callback = cb,
        dt =1e-9*ustrip(τᵢ  |> u"s"),
        dtmin = ustrip(1e-30*τᵢ |> u"s"), #1e-20,
        force_dtmin = false,
        reltol = c_init.alg_ctl.reltol,
        abstol = c_init.alg_ctl.abstol*ones(length(u0)),#u0 .* 0, #1e-12,#c.u0 .* 0,
        maxiters= 5000,
    )

    @debug "Initialized V_oc =" get_V(c,sol)[end] sol.t[end]*c.parameters.τᵢ
    return sol

#=
    prob = SteadyStateProblem(
        odefun,
        u0,
        c_init;
    )
    sol = solve(prob,DynamicSS(Rodas5();abstol=1e-4,reltol=1e-6,tspan= ustrip(1e4u"s"/c.parameters.τᵢ)),progress = true, progress_steps=10,force_dtmin =true ,callback = AutoAbstol(false;init_curmax=u0 .+ 0.1), dtmin=1e-10*ustrip(τᵢ),abstol=u0.*0,reltol=1e-6)
    if sol.retcode !=:Success
        #return u0
    end
    return sol

#    return u1.zero
=#
end


function nl_solve_intiter(c_init::Cell,u0;ftol=1e-6,factor=1)
    # get the sparse colored jacobian for fast NLsolve
    function j!(jac, x, c_init)
        colors = matrix_colors(c_init.Jac)
        forwarddiff_color_jacobian!(
            jac,
            (dx, x) -> c_init.rhs(dx, x, c_init, 0.0),
            x;
            colorvec = colors,
            sparsity = c_init.Jac,
        )
        return nothing
    end

    df = OnceDifferentiable(
        (dx, x) -> c_init.rhs(dx, x, c_init, 0.0),
        (jac, x) -> j!(jac, x, c_init),
        u0,
        copy(u0),
        c_init.Jac,
    )
    u1 = nlsolve(
        df,
        u0;
        iterations = 100,
        ftol = ftol,
        factor = factor,
        show_trace = haskey(ENV, "JULIA_DEBUG"),
        method = :newton,
    )
    @debug "NLsolve: " u1.f_converged u1.iterations u1.residual_norm
    return u1
end
"""
    initial_conditions!(c::Cell)

Inplace mutation variant of `initial_conditions(c::Cell)``
"""
function initial_conditions!(c::Cell)
    c.u0 = initial_conditions(c).u[end]

end

"""
    init_guess(g::Grid,ndim::NodimParameters)

Returns an initial guess for NLsolve root finding. The `NodimParameter` value
is needed to guarantee consitancy of the guess.
"""
function init_guess(g::Grid, ndim::NodimParameters)
    P_init = ones(size(g.x))
    phi_init = zeros(size(g.x))
    p_init = ones(size(g.x)) * ndim.kₕ .* range(0, 1, length = g.N + 1) #* 1e-1
    n_init = ones(size(g.x)) * ndim.kₑ .* range(1, 0, length = g.N + 1) #* 1e-1
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
