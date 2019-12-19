"""
    function calculate_currents(c::Cell, sol, dt, sol_prev)

Calculates the total current in the center of the device see: DOI: 0.1007/s10825-019-01396-2
"""

function calculate_currents(c::Cell, sol, dt, sol_prev)
    if dt == 0
        dt = Inf
    end
    N = c.g.Nₑ
    κₙ = c.ndim.κₑ
    κₚ = c.ndim.κₚ
    dpt = c.ndim.dpt
    dpf = c.ndim.dpf
    kₑ = c.ndim.kₑ
    kₕ = c.ndim.kₕ

    P = sol[1][1]
    ϕ = sol[2][1]
    ϕ_prev = sol_prev[2][1]
    n = sol[3][1]
    p = sol[4][1]

    dx = c.g.dxₑ
    pos = ceil(Int, (N + 1) / 2) + 0# *0+100  # mid of grid


    jn = κₙ ./ dx[pos] .*
         ((n[pos+1] - n[pos]) - (n[pos+1] + n[pos]) .* (ϕ[pos+1] - ϕ[pos]) ./ 2)
#    jp = -κₚ ./ dx[pos] .*
#         (p[pos+1] - p[pos] + (p[pos+1] + p[pos]) .* (ϕ[pos+1] - ϕ[pos]) ./ 2)
    jd = dpt ./ dx[pos] .* (ϕ[pos+1] - ϕ[pos] - ϕ_prev[pos+1] + ϕ_prev[pos]) ./ dt
#    jf = -dpf ./ dx[pos] .*
#         (P[pos+1] - P[pos] + (P[pos+1] + P[pos]) .* (ϕ[pos+1] - ϕ[pos]) ./ 2)

    jₛₕ = -(sol[2][3][end] +c.ndim.Vbi) *c.ndim.σₛₕ
    return jn - jd  + jₛₕ
end
