"""
    function calculate_currents(c::Cell, sol, dt, sol_prev)

Calculates the total current in the center of the device see: DOI: 0.1007/s10825-019-01396-2
"""
function calculate_currents(g,ndim, sol, dt, sol_prev)
    if dt == 0
        dt = Inf
    end
    N = g.Nₑ
    N = g.N
    κₙ = ndim.κₑ
#    κₙ = ndim.κₙ
    κₚ = ndim.κₚ
    dpt = ndim.dpt
    dpf = ndim.dpf
    kₑ = ndim.kₑ
    kₕ = ndim.kₕ

    ϕ = sol[2][1]
    ϕ_prev = sol_prev[2][1]
    n = sol[3][1]
    p = sol[4][1]


#    P = sol[1][1]
#    ϕ = sol[2][2]
#    ϕ_prev = sol_prev[2][2]
#    n = sol[3][2]
#    p = sol[4][1]

    dx = g.dxₑ
#    dx = g.dx
#    pos = ceil(Int, (N + 1) / 2) + 0# *0+100  # mid of grid

    pos =2
    jn = κₙ ./ dx[pos] .*
         ((n[pos+1] - n[pos]) - (n[pos+1] + n[pos]) .* (ϕ[pos+1] - ϕ[pos]) ./ 2)
#    jp = -κₚ ./ dx[pos] .*
#         (p[pos+1] - p[pos] + (p[pos+1] + p[pos]) .* (ϕ[pos+1] - ϕ[pos]) ./ 2)

    jd = dpt ./ dx[pos] .* (ϕ[pos+1] - ϕ[pos] - ϕ_prev[pos+1] + ϕ_prev[pos]) ./ dt
#    jf = -dpf ./ dx[pos] .*
#         (P[pos+1] - P[pos] + (P[pos+1] + P[pos]) .* (ϕ[pos+1] - ϕ[pos]) ./ 2)

    jₛₕ = -(sol[2][3][end] +ndim.Vbi) *ndim.σₛₕ
    return jn+ - jd + jₛₕ
    return jn + jp - jd -jf + jₛₕ
end
calculate_currents(p::Rhs, sol, dt, sol_prev) = calculate_currents(p.g,p.ndim, sol, dt, sol_prev)
function calculate_currents(sol::DiffEqBase.ODESolution)
    p = sol.prob.f.f
    J = Array{typeof(p.parameters.jay |> u"A/m^2")}(undef,length(sol.t))
    u = decompose.(sol.u,p.g)
    J[1] =calculate_currents(p, u[1], Inf, u[1]) * p.parameters.jay
    for (u,dt,u_prev,i) in zip(u[2:end],diff(sol.t),u[1:end-1],1:length(sol.t)-1)
        J[i+1]= calculate_currents(p, u, dt, u_prev) * p.parameters.jay
    end
    return J
end

get_V(sol::DiffEqBase.ODESolution,t::AbstractArray=sol.t .*sol.prob.f.f.parameters.τᵢ) = get_V.((sol,),t)
function get_V(sol::DiffEqBase.ODESolution,t)
    p = sol.prob.f.f
    if !(typeof(t) isa Unitful.AbstractQuantity)
        t = t*u"s"
    end
    # return external voltage if cc
    if p.mode == :cc
        V= p.parameters.V(t)*u"V"
        return V
    end
    u = rdim_sol(sol,t)
    V= u[2][3][end]+p.ndim.Vbi*p.parameters.VT
    return V
end

get_t(sol::DiffEqBase.ODESolution) = upreferred.(sol.t*sol.prob.f.f.parameters.τᵢ )

function (f::DiffEqBase.ODESolution)(t::Unitful.AbstractQuantity)

    p = f.prob.f.f
    sol_vec = rdim_sol(f,t)
    r=Dict{Symbol,Any}()
    r[:x]=p.g.x.*p.parameters.b
    r[:xₑ]=p.g.xₑ.*p.parameters.b
    r[:xₕ]=p.g.xₕ.*p.parameters.b
    r[:P]=sol_vec[1][1]
    r[:ϕ]=sol_vec[2][2]
    r[:ϕₑ]=sol_vec[2][1]
    r[:ϕₕ]=sol_vec[2][3]
    r[:n]=sol_vec[3][2]
    r[:nₑ]=sol_vec[3][1]
    r[:p]=sol_vec[4][1]
    r[:pₕ]=sol_vec[4][2]
    getQFL!(r,p.parameters)
    r[:V]=get_V(f,t)
    #calculate current wo dispacement
    t_ndim =Float64(t/p.parameters.τᵢ)
    ndim_sol=decompose(f(t_ndim), p.g)
    dt = minimum(abs.(f.t .- t_ndim))
    ndim_sol_prev=(f(t_ndim-abs(dt)), p.g)
    ndim_sol_pre=decompose(f(Float64(t/p.parameters.τᵢ)), p.g)
    r[:j]= -calculate_currents(p.g,p.ndim, ndim_sol, dt, ndim_sol)* p.parameters.jay
    return r
end

function getQFL!(d::Dict, p::AbstractParameters)
     d[:nₑf] = @.p.Ecₑ - p.kB * p.T * log(p.gcₑ / abs.(d[:nₑ]))
     d[:nf] = @. p.Ec - p.kB * p.T * log(p.gc / abs.(d[:n]))

     d[:pf] =@. p.Ev + p.kB * p.T * log(p.gv / abs.(d[:p]))
     d[:pₕf] = @.p.Evₕ + p.kB * p.T * log(p.gvₕ / abs.(d[:pₕ]))
     d[:QFLS] = @.d[:nf] - d[:pf]
     d[:V] = d[:ϕₕ][end]+p.Vbi/u"eV"*u"V"
    nothing
end
