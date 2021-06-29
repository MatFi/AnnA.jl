
"""
    initial_conditions(c::Cell)

Calculates the initial condition of a given `Cell` at `t=tstart`.
Returns the steady state solution vector.
"""
function initial_conditions(c::Cell)

    u0 = c.u0
    cc =deepcopy(c)
    ts=[0.0,1e-2,1e256] 
    vs=[ustrip(cc.parameters.Vbi),
        cc.parameters.V(convert(Float64,upreferred(cc.alg_ctl.tstart/cc.parameters.τᵢ))),
        cc.parameters.V(convert(Float64,upreferred(cc.alg_ctl.tstart/cc.parameters.τᵢ))),
    ]
    v_itp = interpolate(ts, vs ,SteffenMonotonicInterpolation())
    p_init = setproperties(
        cc.parameters, 
        V = (t) -> v_itp(t),
        light=t->c.parameters.light(ustrip(upreferred(cc.alg_ctl.tstart))),
    )
    if ustrip(p_init.Rₛ)==Inf
        @warn "initialization could not be performed with Rₛ=Inf, used 0 instead"
        p_init.Rₛ=0u"V/A*m^2"
    end
    c_init = Cell(p_init;mode = :cc, alg_ctl = cc.alg_ctl, u0=Float64.(cc.u0))

    if  c.mode == :occ  #legacy
        @info "initalisatiion: stating conditions in :oc mode"
        u0 = nl_solve_intiter(c,u0;ftol=c.alg_ctl.ss_tol).zero
    end

    @debug "Init_Solve"
    τᵢ = c_init.parameters.τᵢ
    prob = get_problem(c_init,tstart=0.0u"s",tend = 1e6*τᵢ |> u"s")

    ss_cb = TerminateSteadyState(c.alg_ctl.ss_tol,c.alg_ctl.ss_tol)
    

    cb = CallbackSet(ss_cb,)#abstol_cb )
    sol = solve(prob,c_init.alg_ctl.alg;
        progress_steps = 10,
        progress = c_init.alg_ctl.progress,
        callback = cb,
        dt =1e-9*ustrip(τᵢ  |> u"s"),
        dtmin = ustrip(1e-40*τᵢ |> u"s"), #1e-20,
        force_dtmin = true,
        reltol = c_init.alg_ctl.reltol,
        abstol = c_init.alg_ctl.abstol,#*ones(length(u0)),#u0 .* 0, #1e-12,#c.u0 .* 0,
        maxiters= 5000,
        save_on=false,
        #initializealg=ShampineCollocationInit(),
        initializealg=OrdinaryDiffEq.NoInit(),
       # initializealg=BrownFullBasicInit(),
    )

        if sol.retcode ==:Success
            @warn "Initialisation did not reach steady state within ss_tol"
        elseif sol.retcode  !=:Terminated
            @warn "Initialisation failed"
            #@show "Initialized V_oc =" get_V(c,sol)[end] sol.t[end]*c.parameters.τᵢ
        end
    return sol

end


# function nl_solve_intiter(c_init::Cell,u0;ftol=1e-6,factor=1)
#     # get the sparse colored jacobian for fast NLsolve
#     function j!(jac, x, c_init)
#         colors = matrix_colors(c_init.Jac)
#         forwarddiff_color_jacobian!(
#             jac,
#             (dx, x) -> c_init.rhs(dx, x, c_init, 0.0),
#             x;
#             colorvec = colors,
#             sparsity = c_init.Jac,
#         )
#         return nothing
#     end

#     df = OnceDifferentiable(
#         (dx, x) -> c_init.rhs(dx, x, c_init, 0.0),
#         (jac, x) -> j!(jac, x, c_init),
#         u0,
#         copy(u0),
#         c_init.Jac,
#     )
#     u1 = nlsolve(
#         df,
#         u0;
#         iterations = 10000,
#         ftol = ftol,
#         factor = factor,
#         show_trace = haskey(ENV, "JULIA_DEBUG"),
#        # method = :newton,
#     )
#     @debug "NLsolve: " u1.f_converged u1.iterations u1.residual_norm
#     return u1
# end
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
function init_guess(g::Grid, ndim::NodimParameters,Vbikt)
    P_init = ones(size(g.x))
    phi_init = zeros(size(g.x))
    dn,dp = (0.0,0.0)

    if ndim.ϰ >0
        dn =  abs(ndim.ϰ )/ ndim.δ
        dp = ndim.nᵢ²*exp(Vbikt)/dn
    else
        dp = abs(ndim.ϰ)/(ndim.δ*ndim.χ)
        dn = ndim.nᵢ²*exp(Vbikt)/dp
    end

    dnend =   ndim.nᵢ²*exp(Vbikt)/ndim.kₕ
    n_init = ones(size(g.x)) .*range( ndim.kₑ,dnend, length = g.N + 1) #* 
    p_init = ndim.nᵢ²*exp(Vbikt)./n_init

    #  n_init = ones(size(g.x)) .*dn #* 1
    # p_init = ones(size(g.x)) .*dp#* 1e-1
    #  p_init = ones(size(g.x)) .*range(dp , ndim.kₕ, length = g.N + 1) #* 1e-1
    # n_init = ones(size(g.x)) .*range( ndim.kₑ,dn, length = g.N + 1) #* 1e-1
    phiE_init = zeros(size(g.xₑ))
    nE_init = ones(size(g.xₑ))# .* range(1, n_init[1] / ndim.kₑ, length = g.Nₑ)
    phiH_init = zeros(size(g.xₕ))
    pH_init = ones(size(g.xₕ)) #.* range(p_init[end] / ndim.kₕ, 1, length = g.Nₕ)
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