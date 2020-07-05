struct JscVocProblem{
    V<:Number,
    VR<:Number,
    P<:Parameters,
    C<:AlgControl,
} <: AbstractProblem
    parameters::P
    max_intensity::V
    max_t::VR
    alg_control::C
end

"""

"""
function JscVocProblem(
    parm::AbstractParameters,
    max_intensity::Number,
    max_t::Number;
    alg_control = missing,
)

    if alg_control isa Missing

        alg_control = AlgControl(
            dtmin = 1e-20,
            dt = 1e-15,
            reltol = 1e-5,
            abstol = 1e-8,
            tend = max_t,
            force_dtmin=false,
        )
    else
        #enforce correct tend
        alc_control = setproperty!(
            alg_control,
            :tend,copy(max_t)
        )
    end
    max_t=ustrip(upreferred(copy(max_t)))
    max_intensity=ustrip(upreferred(copy(max_intensity)))
#   lightt = (t) -> max_intensity*exp(-ustrip(max_t )/(t+eps(Float64)) +1)
#    lightt = (t) -> max_intensity*(exp(-ustrip(t/max_t ))-1)/(exp(1)-1)
#    lightt = (t) -> max_intensity*((exp(8*log(max_intensity)*t/ustrip(max_t) )-1)/(exp(8*log(max_intensity))-1))^2

#   lightt = (t) -> max_intensity*(ustrip(t/max_t))^4

   i0=1e-16
   k=log((max_intensity+i0)/i0)
   lightt = (t) -> i0*(exp(k*(t/max_t))-1)

    parm = setproperties(parm, light = lightt, V=t->0)
    prob = JscVocProblem(parm, max_intensity, max_t, alg_control)

end

struct JscVocSolution <: AbstractProblemSolution
    sol_oc::Union{DiffEqBase.ODESolution}
    sol_cc::Union{DiffEqBase.ODESolution}
    prob::JscVocProblem
    Jsc::AbstractArray
    Voc::AbstractArray
    Illumination::AbstractArray
end

JscVocSolution(sol_oc,sol_cc, p::JscVocProblem) = JscVocSolution(
    sol_oc,
    sol_cc,
    p,
    #get_V(sol_cc;t=sol_oc.t*p.parameters.τᵢ)./1e-12u"V/A*m^2",
    calculate_currents(sol_cc),
    get_V(sol_oc,t=sol_cc.t*p.parameters.τᵢ),
    p.parameters.light.(ustrip(sol_cc.t*p.parameters.τᵢ))
)



function solve(p::JscVocProblem; alg_control = p.alg_control,kwargs...)
    tend = p.max_t
    #init

    τᵢ = p.parameters.τᵢ
    if alg_control isa Missing
        alg_control = AlgControl(
            dtmin = 1e-22*ustrip(τᵢ  |> u"s"),
            dt = 1e-8*ustrip(τᵢ  |> u"s"),
            reltol = 1e-8,
            abstol = 1e-8,
            force_dtmin=false,
            tend = tend,
        )
    end
    tstops=10.0 .^(-13:0.3:log10(ustrip(tend/p.parameters.τᵢ)))
    init_oc = Cell(
        p.parameters,
        mode = :oc,
        alg_ctl = alg_control
    )
    sol_oc = solve(init_oc,tstops=tstops)

    init_cc = Cell(
    #    p.parameters,
        setproperties(p.parameters, Rₛₕ= Inf*u"V/A*m^2"),
        mode = :cc,
        alg_ctl = alg_control
    )
    init_cc.initialized=:true
    init_cc.u0=init_oc.u0
    sol_cc = solve(init_cc,tstops=sol_oc.t )
    return JscVocSolution(sol_oc,sol_cc, p)
end
