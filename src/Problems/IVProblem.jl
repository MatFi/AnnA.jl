mutable struct IVProblem{B<:Bool,V<:AbstractArray,VR<:Number,P<:Parameters,C<:AlgControl}
    parameters::P
    voltage_range::V
    sweep_rate::VR
    double_sweep::B
    alg_control::C
end

"""
    IVProblem(
        parm::Parameters,
        range::AbstractArray,
        rate::Unitful.AbstractQuantity;
        double_sweep = true,
        alg_control = AlgControl(dtmin = 1e-20,
            dt = 1e-6,
            reltol = 1e-4,
            abstol = 1e-12,
            tend = (maximum(range) - minimum(range)) / abs(rate)
        ),
    )

Creates an `IVProblem`. The voltage parameter `V` defined in `parm` gets overwritten by the linear voltage sweep defined as `V = t-> first(range) + rate*t`. `AlgControl.tend` is forced to be `(maximum(range) - minimum(range)) / abs(rate)`, an overwrite can be done on the finalized object using Setfield:
    `p = Setfield.setproperties(p::IVProblem, alg_control=AlgControl(...))`
"""
function IVProblem(
    parm::AbstractParameters,
    range::AbstractArray,
    rate::Unitful.AbstractQuantity;
    double_sweep = true,
    alg_control = missing,
)

    if alg_control isa Missing
        alg_control = AlgControl(
            dtmin = 1e-20,
            dt = 1e-6,
            reltol = 1e-4,
            abstol = 1e-12,
            tend = (maximum(range) - minimum(range)) / abs(rate),
        )
    else
        #enforce correct tend
        alc_control = setproperty!(
            alg_control,
            tend = (maximum(range) - minimum(range)) / abs(rate),
        )
    end

    Vt = t -> ustrip(uconvert(u"V", first(range))) +
              sign(-ustrip(first(range) + last(range))) *
              abs(ustrip(uconvert(u"V/s", rate))) * ustrip(t)

    parm = setproperties(parm, V = Vt)
    prob = IVProblem(parm, range, rate, double_sweep, alg_control)

end

struct IVSolution <: AbstractProblemSolution
    sol_fwd::Union{DiffEqBase.ODESolution,Nothing}
    sol_rwd::Union{DiffEqBase.ODESolution,Nothing}
    prob::IVProblem
    I_fwd::Union{AbstractArray,Nothing}
    I_rwd::Union{AbstractArray,Nothing}
    V_fwd::Union{AbstractArray,Nothing}
    V_rwd::Union{AbstractArray,Nothing}
end

IVSolution(fwd, rwd, p::IVProblem) = IVSolution(
    fwd,
    rwd,
    p,
    calculate_currents(fwd),
    calculate_currents(rwd),
    get_V(fwd),
    get_V(rwd),
)

IVSolution(fwd, rwd::Nothing, p::IVProblem) =
    IVSolution(fwd, nothing, p, calculate_currents(fwd), nothing, get_V(fwd), nothing)

IVSolution(fwd::Nothing, rwd, p::IVProblem) =
    IVSolution(nothing, rwd, p, nothing, calculate_currents(rwd), nothing, get_V(rwd))


function solve(p::IVProblem, args...)
    tend = (maximum(p.voltage_range) - minimum(p.voltage_range)) / abs(p.sweep_rate)
    #init
    init_c = Cell(
        p.parameters,
        mode = :cc,
        alg_ctl = AlgControl(
            dtmin = 1e-20,
            dt = 1e-6,
            reltol = 1e-4,
            abstol = 1e-12,
            tend = tend,
        ),
    )
    s1 = solve(init_c)

    if p.double_sweep == true
        p2 = setproperties(
            p.parameters,
            V = t -> p.parameters.V(ustrip(tend |> u"s") - ustrip(t)),
        )
        init_c = Cell(
            p2;
            u0 = s1.u[end],
            mode = :cc,
            alg_ctl = AlgControl(
                dtmin = 1e-20,
                dt = 1e6,
                reltol = 1e-6,
                abstol = 1e-8,
                tend = tend,
            ),
        )
        s2 = solve(init_c)#.u[end]
        (fwd, rwd) = p.parameters.V(0) < p2.V(0) ? (s1, s2) : (s2, s1)

        return IVSolution(fwd, rwd, p)
    end

    (fwd, rwd) = p.parameters.V(0) < p.parameters.V(tend) ? (s1, nothing) : (nothing, s1)
    return IVSolution(fwd, rwd, p)
end
