struct OCVDProblem{P<:AbstractParameters,T<:Number,C<:AlgControl}  <: AbstractProblem
    parameters::P
    on_time::T
    decay_time::T
    alg_control::C
end

function OCVDProblem(
    parm::AbstractParameters,
    on_time::Unitful.AbstractQuantity,
    decay_time::Unitful.AbstractQuantity;
    alg_control = missing,
)
    on_time = (on_time |> u"s")
    if alg_control isa Missing
        alg_control = AlgControl(
            dtmin = ustrip(1e-20 * parm.τᵢ),
            dt = ustrip(1e-12 * parm.τᵢ),
            reltol = 1e-4,
            abstol = 1e-12,
            tend = on_time ,
        )
    else
        #enforce correct tend
        alc_control = setproperty!(alg_control, tend = ustrip(on_time / parm.τᵢ))
    end

    lght = pulse(
        tₑ = Float64(ustrip(on_time)) + 2e-15,
        w = Float64(ustrip(on_time)),
        Δh = 1.0,
        Δt = 1e-15,
    )
    parm = setproperties(parm; light =lght, V= t->0)

    prob = OCVDProblem(parm, promote(on_time, decay_time)..., alg_control)

end

struct OCVDSolution{S<:DiffEqBase.ODESolution,P<:OCVDProblem,A,B} <: AbstractProblemSolution
    sol_on::S
    sol_decay::S
    prob::P
    t_on::Array{A,1}
    t_decay::Array{A,1}
    V_on::Array{B,1}
    V_decay::Array{B,1}
end

OCVDSolution(on, decay, p::OCVDProblem) = OCVDSolution(
    on,
    decay,
    p,
    get_t(on) ,
    get_t(decay),
    get_V(on),
    get_V(decay),
)

function solve(p::OCVDProblem, args...)

    tend = p.on_time
    parms= p.parameters

    init_c = Cell(parms;mode = :oc,alg_ctl = p.alg_control)


    s1 = solve(init_c)
    t0=ustrip(upreferred(copy(p.on_time)))
    p2 = setproperties(
        p.parameters,
        light  = t -> p.parameters.light(t+t0),
    )
    alg_ctl =  setproperties(
            p.alg_control,
            tend= p.decay_time,
        )
    c = Cell(p2;u0 = s1.u[end],mode = :oc,alg_ctl =alg_ctl)

    s2 = solve(c)#.u[end]


    return OCVDSolution(s1, s2, p)

end
