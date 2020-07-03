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
    max_t::Unitful.AbstractQuantity;
    alg_control = missing,
)
    if alg_control isa Missing

        alg_control = AlgControl(
            dtmin = 1e-20,
            dt = 1e-8,
            reltol = 1e-5,
            abstol = 1e-8,
            tend = max_intensity / abs(rate),
        )
    else
        #enforce correct tend
        alc_control = setproperty!(
            alg_control,
            :tend,max_t
        )
    end

    lightt = t -> max_intensity*exp(-ustrip(max_t |> u"s")/(t) +1)



    parm = setproperties(parm, light = lightt, V=t->0)
    prob = JscVocProblem(parm, max_intensity, max_t, alg_control)

end

struct JscVocSolution <: AbstractProblemSolution
    sol::Union{DiffEqBase.ODESolution,Nothing}
    prob::IVProblem
    Jsc::Union{AbstractArray,Nothing}
    Voc::Union{AbstractArray,Nothing}
end

JscVocSolution(sol_oc,sol_cc, p::JscVocProblem) = JscVocSolution(
    sol,
    p,
    calculate_currents(sol_oc),
    get_V(sol_cc),
)



function solve(p::JscVocProblem, alg_control = p.alg_control, args...)
    tend = p.max_t
    #init

    τᵢ = p.parameters.τᵢ
    if alg_control isa Missing
        alg_control = AlgControl(
            dtmin = 1e-22*ustrip(τᵢ  |> u"s"),
            dt = 1e-8*ustrip(τᵢ  |> u"s"),
            reltol = 1e-6,
            abstol = 1e-6,
            force_dtmin=false,
            tend = tend,
        )
    end

    init_c = Cell(
        p.parameters,
        mode = :oc,
        alg_ctl = alg_control
    )
    sol_oc = solve(init_c)

    init_c = Cell(
        p.parameters,
        mode = :cc,
        alg_ctl = alg_control
    )
    sol_cc = solve(init_c)
    return IVSolution(sol_oc,sol_cc, p)
end
