struct JscVocProblem{V<:Number,VR<:Number,P<:Parameters,C<:AlgControl} <: AbstractProblem
    parameters::P
    max_intensity::V
    max_t::VR
    alg_control::C
end


"""
    function JscVocProblem(
        parm::AbstractParameters,
        max_intensity::Number,
        max_t::Number;
        alg_control = missing,
    )

Creates an `JscVocProblem`. Where the backgrond illumination is increased to `max_intensity*Fₚₕ` logerythmical within the time `max_t`.
"""
function JscVocProblem( 
    parm::AbstractParameters,
    max_intensity::Number,
    max_t::Number;
    alg_control = missing,
)

    τᵢ = parm.τᵢ
    if alg_control isa Missing
        alg_control = AlgControl(
            dtmin = 1e-22 * ustrip(τᵢ |> u"s"),
            dt = 1e-8 * ustrip(τᵢ |> u"s"),
            reltol = 1e-6,
            abstol = 1e-6,
            force_dtmin = false,
            tend = max_t,
        )

    else
        #enforce correct tend
        alc_control = setproperty!(alg_control, :tend, copy(max_t))
    end
    max_t = ustrip(upreferred(copy(max_t)))
    max_intensity = ustrip(upreferred(copy(max_intensity)))

    i0 = 1e-16
    k = log((max_intensity + i0) / i0)
    lightt = (t) -> i0 * (exp(k * (t / max_t)) - 1)

    parm = setproperties(parm, light = lightt, V = t -> 0)
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

JscVocSolution(sol_oc, sol_cc, p::JscVocProblem) = JscVocSolution(
    sol_oc,
    sol_cc,
    p,
    #get_V(sol_cc;t=sol_oc.t*p.parameters.τᵢ)./1e-12u"V/A*m^2",
    calculate_currents(sol_cc),
    get_V(sol_oc, t = sol_cc.t * p.parameters.τᵢ),
    p.parameters.light.(ustrip(sol_cc.t * p.parameters.τᵢ)),
)



function solve(p::JscVocProblem; alg_control = p.alg_control, kwargs...)
    tend = p.max_t
    #init

    tstops = 10.0 .^ (-13:0.3:log10(ustrip(tend / p.parameters.τᵢ)))
    init_oc = Cell(p.parameters, mode = :oc, alg_ctl = alg_control)
    sol_oc = solve(init_oc, tstops = tstops)

    init_cc = Cell(
        #    p.parameters,
        setproperties(p.parameters, Rₛₕ = Inf * u"V/A*m^2"),
        mode = :cc,
        alg_ctl = alg_control,
    )
    init_cc.initialized = :true
    init_cc.u0 = init_oc.u0
    sol_cc = solve(init_cc, tstops = sol_oc.t)
    return JscVocSolution(sol_oc, sol_cc, p)
end
