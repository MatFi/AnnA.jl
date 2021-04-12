struct ProblemSolution <: AbstractProblemSolution
    parameters
    df
    xₑ
    x
    xₕ
end
function ProblemSolution(sol::ODESolution)
    parm = sol.prob.f.f
    grid = parm.g
    ProblemSolution(parm.parameters, todf(sol), grid.xₑ, grid.x, grid.xₕ)
end



function todf(sol::ODESolution)
    rhs = sol.prob.f.f
    parm = rhs.parameters
    grid = rhs.g
    τᵢ = parm.τᵢ

    j = -calculate_currents(sol)
    V = get_V(sol)
    t = sol.t * τᵢ .|> u"s" 

    d = rdim_sol(sol)
  #  d= decompose(s_redim, grid)
    
    phi_etm = map(a -> a[2][1], d)
    n_etm =  map(a -> a[3][1], d)
    P = map(a -> a[1][1], d)
    phi = map(a -> a[2][2], d)
    n =  map(a -> a[3][2], d)
    p = map(a -> a[4][1], d)
    phi_htm = map(a -> a[2][3], d)
    p_htm = map(a -> a[4][2], d)
   
    Ecₑ = parm.Ecₑ
    Ec = parm.Ec
    Ev = parm.Ev
    Evₕ = parm.Evₕ

    Ece = map(x->(-x * u"eV/V").+Ecₑ ,phi_etm) 
    Ecp = map(x->(-x * u"eV/V") .+ Ec,phi)
    Evp = map(x->(-x * u"eV/V" ) .+ Ev, phi)
    Evh = map(x-> (-x * u"eV/V") .+ Evₕ,phi_htm)
  
    nₑf =map((e,n)-> (e .- parm.kB .* parm.T .* log.(parm.gcₑ ./ abs.(n))),Ece,n_etm)
    nf = map((e,n)-> (e .- parm.kB .* parm.T .* log.(parm.gc ./ abs.(n))),Ecp,n)

    pf = map((e,p)-> (e .+ parm.kB .* parm.T .* log.(parm.gv ./ abs.(p))),Evp,p)
    pₕf = map((e,p) -> (e .+ parm.kB .* parm.T .* log.(parm.gvₕ ./ abs.(p))),Evh,p_htm)
    

    light = parm.light.(ustrip(sol.t * τᵢ .|> u"s"))

    df = DataFrame(
        t=t,
        V=V,
        j=j,
        light=light,
        nₑ=n_etm,
        n=n,
        p=p,
        pₕ=p_htm,
        ϕₑ=phi_etm,
        ϕ=phi,
        ϕₕ=phi_htm,
        P=P,
        Ecₑ=Ece,
        Ec=Ecp,
        Ev=Evp,
        Evₕ=Evh,
        nₑf=nₑf,
        nf=nf,
        pf=pf,
        pₕf=pₕf,
    )
end



