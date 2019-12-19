function rdim_sol(c::Cell,sol)
#    [[P], [ϕₑ,ϕ,ϕₕ], [nₑ,n], [p,pₕ]]
    P   = sol[1].*c.parameters.N₀
    ϕ   = sol[2].*c.parameters.VT
    n   = sol[3].*c.parameters.dₑ
    p   = sol[4].*c.parameters.dₕ

    return  [P, ϕ ,n, p]
end
