"""
    rdim_sol(sol,p::AbstractParameters)

a single time slice `sol` of the ODESolution is redimensionalized with `p`
"""
function rdim_sol(sol,p::AbstractParameters)
#    [[P], [ϕₑ,ϕ,ϕₕ], [nₑ,n], [p,pₕ]]
    P   = sol[1].*p.N₀
    ϕ   = sol[2].*p.VT
    n   = sol[3].*p.dₑ
    p   = sol[4].*p.dₕ

    return  [P, ϕ ,n, p]
end

rdim_sol(c::Cell,sol::Array{T,1} where T<:Real) = rdim_sol(decompose(sol,c.g), c.parameters) # legacy wrapper
#rdim_sol(sol::Vector{Float64},p::AbstractParameters) = rdim_sol(sol,p)

"""
    rdim_sol(sol::DiffEqBase.ODESolution)

Redimensionalizes the complete ODESolution at each time step.
"""
function rdim_sol(sol::DiffEqBase.ODESolution)
    p = sol.prob.f.f
    u = decompose.(sol.u,p.g)
    u_ret = rdim_sol.(u,(p.parameters,))
    return u_ret
end
function rdim_sol(c::Cell,sol::DiffEqBase.ODESolution)
    p = c.rhs
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

function rdim_sol(sol::DiffEqBase.ODESolution,t::Number)
    p = sol.prob.f.f
    u = decompose(sol(ustrip(upreferred(t/p.parameters.τᵢ)|>ustrip)),p.g)
    u_ret = rdim_sol(u,p.parameters)
    return u_ret
end
