struct OCVD end
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
    τᵢ = parm.τᵢ
    on_time = (on_time |> u"s")
    if alg_control isa Missing
        alg_control = AlgControl(
            dt =1e-9*ustrip(τᵢ  |> u"s"),
            dtmin = ustrip(1e-40*τᵢ |> u"s"), #1e-20,
            reltol = 1e-5,
            abstol = 1e-8,
            force_dtmin=false,
            ss_tol=1e-5,
            tend = on_time ,
        )
    else
        #enforce correct tend
        alc_control = setproperty!(alg_control, :tend,  decay_time )
    end
    
    lght = pulse(
        tₑ = 0.,#Float64(ustrip(on_time)) + 2e-15,
        w = Float64(ustrip(on_time))-1e-12,
        Δh = 1.0,
        Δt = 1e-12,
    )

    parm = setproperties(parm; light =lght, V= (t)->0, Rₛ= Inf*u"V/A*m^2")

    prob = OCVDProblem(parm, promote(on_time, decay_time)..., alg_control)

end

# backwards compatibility and conveniant dispatches
getproperty(p::ProblemSolution{OCVD},n::Symbol) = getproperty(p, Val{n}())
getproperty(p::ProblemSolution{OCVD},::Val{S}) where {S} = getfield(p, S)
getproperty(p::ProblemSolution{OCVD},::Val{:t_on}) = begin r=filter(row -> row.t <= 0u"s", p.df); r.t .-=r.t[1]; r.t end
getproperty(p::ProblemSolution{OCVD},::Val{:t_decay}) = filter(row -> row.t >= 0u"s", p.df).t
getproperty(p::ProblemSolution{OCVD},::Val{:V_on}) = begin r=filter(row -> row.t <= 0u"s", p.df); r.t .-=r.t[1]; r.V end
getproperty(p::ProblemSolution{OCVD},::Val{:V_decay}) = filter(row -> row.t >= 0u"s", p.df).V
getproperty(p::ProblemSolution{OCVD},::Val{:sol_on}) =begin r=filter(row -> row.t <= 0u"s", p.df); r.t .-=r.t[1]; r end
getproperty(p::ProblemSolution{OCVD},::Val{:sol_decay}) = filter(row -> row.t >= 0u"s", p.df)


function solve(pp::OCVDProblem, args...)
    p=deepcopy(pp)
    tend = p.on_time
    parms= p.parameters
    alg_ctl =  setproperties(
            p.alg_control,
            tend= p.decay_time,
            tstart= -p.on_time
        )
    c = Cell(parms;mode = :oc,alg_ctl =alg_ctl)
    sol = solve(c)#.u[end]
    return ProblemSolution(sol, OCVD)
end
