"""
    rdim_sol(u,p::AbstractParameters)

a single time slice `u` of the ODESolution is redimensionalized with `p`
"""
function rdim_sol(u,p::AbstractParameters)
#    [[P], [ϕₑ,ϕ,ϕₕ], [nₑ,n], [p,pₕ]]
    P   = u[1].*p.N₀
    ϕ   = u[2].*p.VT
    n   = u[3].*p.dₑ
    p   = u[4].*p.dₕ
    return  [P, ϕ ,n, p]
end

"""
    rdim_sol(sol::DiffEqBase.ODESolution)

Redimensionalizes the complete ODESolution at all time steps.
"""
function rdim_sol(sol::DiffEqBase.ODESolution)
    p = sol.prob.f.f
    u = decompose.(sol.u,p.g)
    u_ret = rdim_sol.(u,(p.parameters,))
    return u_ret
end
"""
    rdim_sol(sol::DiffEqBase.ODESolution,t)

returns a redimensionalized solution vector at time `t`, the time is treated as
seconds if no unit is given.
"""
rdim_sol(sol::DiffEqBase.ODESolution,t::AbstractArray)= rdim_sol.((sol,),t)

function rdim_sol(sol::DiffEqBase.ODESolution,t)
    if !(t isa Unitful.AbstractQuantity)
        t = t*u"s"
    end
    p = sol.prob.f.f
    u = decompose(sol(ustrip(upreferred(t/p.parameters.τᵢ)|>ustrip)),p.g)
    u_ret = rdim_sol(u,p.parameters)
    return u_ret
end
