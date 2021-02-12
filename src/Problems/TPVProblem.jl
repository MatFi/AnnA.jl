struct TPVProblem{P,C}
    parameters::P
    alg_control::C
end


function TPVProblem(
    parameters,
    I₀,
    ΔI,
    tend;
    alg_control = AnnABase.AlgControl(
        dtmin = 1e-22,
        dt = 1e-20,
        reltol = 1e-8,
        abstol = 1e-8,
        tend = 1*u"s",
        ss_tol=  1e-11
        ),
    Δt=80e-12,
    Fₚₕ=1.4e21u"m^-2*s^-1",
)

    alg_control =setproperties(alg_control,tend=tend*1u"s")
    p=pulse(; Δh=ΔI, h₀=I₀, w=1e-12, Δt=Δt, tₑ=Δt+2e-12)
    parameters= setproperties(parameters,light=(t)->p(t),Fₚₕ=Fₚₕ)

    TPVProblem(parameters,alg_control)
end

function solve(p::TPVProblem, kwargs...)
    c=Cell(p.parameters;mode = :oc,alg_ctl = p.alg_control)
    solve(c,kwargs...)
end