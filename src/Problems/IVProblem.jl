mutable struct IVProblem{B<:Bool,V<:AbstractArray,VR<:Number,P<:Parameters}
    parameters::P
    voltage_range::V
    sweep_rate::VR
    double_sweep::B
    sol
end

"""
    IVProblem(
        parm::Parameters,
        range::AbstractArray,
        rate::Unitful.AbstractQuantity;
        double_sweep = true,
    )

Creates an `IVProblem`. The voltage parameter `V` defined in `parm` gets overwritten by the linear voltage sweep defined as `V = t-> first(range) + rate*t`

"""
function IVProblem(
    parm::Parameters,
    range::AbstractArray,
    rate::Unitful.AbstractQuantity;
    double_sweep = true,
)

    Vt = t -> ustrip(uconvert(u"V", first(range))) + ustrip(uconvert(u"V/s", rate)) * t

    parm = setproperties(parm, V = Vt)
    prob = IVProblem(parm, range, rate, double_sweep, double_sweep ? [] : nothing)

end

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
    s1 = solve(init_c)#.u[end]
    if p.double_sweep==true
        p2 = setproperties(p.parameters, V = t -> p.parameters.V(ustrip(tend |> u"s") - t))
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
        s2 = solve(init_c)
        s1 = (s1,s2)
    end
    return s1
end

function solve!(p::IVProblem, args...)
    p.sol = solve(p,args...)
    return nothing
end
