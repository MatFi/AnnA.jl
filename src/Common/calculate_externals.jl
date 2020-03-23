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
    J = Array{typeof(p.parameters.jay |> u"mA/cm^2")}(undef,length(sol.t))
    u = decompose.(sol.u,p.g)
    J[1] =calculate_currents(p, u[1], Inf, u[1]) * p.parameters.jay
    for (u,dt,u_prev,i) in zip(u[2:end],diff(sol.t),u[1:end-1],1:length(sol.t)-1)
        J[i+1]= calculate_currents(p, u, dt, u_prev) * p.parameters.jay
    end
    return J
end

function get_V(sol::DiffEqBase.ODESolution)#::Array{Quantity{Float64,𝐋^2*𝐌*𝐈^-1*𝐓^-3,Unitful.FreeUnits{(Unitful.V,),𝐋^2*𝐌*𝐈^-1*𝐓^-3,nothing}},1}
    p = sol.prob.f.f
    V =Array{typeof(p.ndim.Vbi*p.parameters.VT |>u"V")}(undef,length(sol.t))
#    u = decompose.(sol.u,p.g)
    u = rdim_sol(sol)
    for (uu,i) in zip(u,eachindex(sol.t))
        V[i]= uu[2][3][end]+p.ndim.Vbi*p.parameters.VT
    end
    return V
end
function get_V(c::Cell,sol::DiffEqBase.ODESolution)
    p = sol.prob.f.f
    p = c.rhs
    V = Array{typeof(p.ndim.Vbi*p.parameters.VT |>u"V")}(undef,length(sol.t))
#    u = decompose.(sol.u,p.g)
    u = rdim_sol(c,sol)
    for (uu,i) in zip(u,eachindex(sol.t))
        V[i]= uu[2][3][end]+p.ndim.Vbi*p.parameters.VT
    end
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
    r[:V]=sol_vec[2][3][end]+p.ndim.Vbi*p.parameters.VT
    return r
end

function getQFL!(d::Dict, p::AbstractParameters)
     d[:nf] = @. p.Ec - p.kB * p.T * log(p.gc / d[:n])
     d[:nₑf] = @.p.Ecₑ - p.kB * p.T * log(p.gcₑ / d[:nₑ])
     d[:pf] =@. p.Ev + p.kB * p.T * log(p.gv / d[:p])
     d[:pₕf] = @.p.Evₕ + p.kB * p.T * log(p.gvₕ / d[:pₕ])
     d[:QFLS] = @.d[:nf] - d[:pf]
    nothing
end
