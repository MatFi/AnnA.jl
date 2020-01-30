
function rdim_sol(sol,p::AbstractParameters)
#    [[P], [ϕₑ,ϕ,ϕₕ], [nₑ,n], [p,pₕ]]
    P   = sol[1].*p.N₀
    ϕ   = sol[2].*p.VT
    n   = sol[3].*p.dₑ
    p   = sol[4].*p.dₕ

    return  [P, ϕ ,n, p]
end

rdim_sol(c::Cell,sol) = rdim_sol(sol, c.parameters) # legacy wrapper
rdim_sol(sol::Vector{Float64},p::AbstractParameters) = rdim_sol(decompose(sol,c.g),p)
function rdim_sol(sol::DiffEqBase.ODESolution)
    p = sol.prob.f.f
    u = decompose.(sol.u,p.g)
    u_ret = rdim_sol.(u,(p.parameters,))
    return u_ret
end
