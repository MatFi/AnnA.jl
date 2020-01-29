"""
    function calculate_currents(c::Cell, sol, dt, sol_prev)

Calculates the total current in the center of the device see: DOI: 0.1007/s10825-019-01396-2
"""
function calculate_currents(c::Cell, sol, dt, sol_prev)
    if dt == 0
        dt = Inf
    end

    N = c.g.N
    κₙ = c.ndim.κₑ
    κₙ = c.ndim.κₙ
    κₚ = c.ndim.κₚ
    dpt = c.ndim.dpt
    dpf = c.ndim.dpf
    kₑ = c.ndim.kₑ
    kₕ = c.ndim.kₕ

    P = sol[1][1]
    ϕ = sol[2][2]
    ϕ_prev = sol_prev[2][2]
    n = sol[3][2]
    p = sol[4][1]


    dx = c.g.dx
    pos = ceil(Int, (N + 1) / 2) + 0# *0+100  # mid of grid

    pos =1
    jn = κₙ ./ dx[pos] .*
         ((n[pos+1] - n[pos]) - (n[pos+1] + n[pos]) .* (ϕ[pos+1] - ϕ[pos]) ./ 2)
    jp = -κₚ ./ dx[pos] .*
         (p[pos+1] - p[pos] + (p[pos+1] + p[pos]) .* (ϕ[pos+1] - ϕ[pos]) ./ 2)
    jd = dpt ./ dx[pos] .* (ϕ[pos+1] - ϕ[pos] - ϕ_prev[pos+1] + ϕ_prev[pos]) ./ dt
    jf = -dpf ./ dx[pos] .*
         (P[pos+1] - P[pos] + (P[pos+1] + P[pos]) .* (ϕ[pos+1] - ϕ[pos]) ./ 2)

    jₛₕ = -(sol[2][3][end] +c.ndim.Vbi) *c.ndim.σₛₕ
    return jn+ jp -jf - jd + jₛₕ
end

function calculate_currents(sol::DiffEqBase.ODESolution)
    p = sol.prob.p
    J = Array{typeof(p.parameters.jay |> u"mA/cm^2")}(undef,length(sol.t)-1)
    u = decompose(sol.u,p.g)
    for (u,dt,u_prev,i) in zip(u[2:end],diff(sol.t),u[1:end-1],1:length(sol.t)-1)
        J[i]= calculate_currents(p, u, dt, u_prev) * p.parameters.jay
    end
    return J
end
