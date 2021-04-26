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
            dtmin = ustrip((1e-30 * parm.τᵢ )|>u"s"),
            dt = ustrip((1e-10 * parm.τᵢ)|>u"s"),
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

    parm = setproperties(parm; light =lght, V= (t)->0)

    prob = OCVDProblem(parm, promote(on_time, decay_time)..., alg_control)

end

struct OCVDSolution{S,P} <: AbstractProblemSolution
    df::S
    #sol_decay::S
    prob::P
    #t_on::Array{A,1}
    #t_decay::Array{A,1}
    #V_on::Array{B,1}
    #V_decay::Array{B,1}
end


# backwards compatibility

getproperty(p::OCVDSolution,n::Symbol) = getproperty(p::OCVDSolution, Val{n}())
getproperty(p::OCVDSolution,::Val{S}) where {S} = getfield(p, S)
getproperty(p::OCVDSolution,::Val{:t_on}) = begin r=filter(row -> row.t <= 0u"s", p.df); r.t .-=r.t[1]; r.t end
getproperty(p::OCVDSolution,::Val{:t_decay}) = filter(row -> row.t >= 0u"s", p.df).t
getproperty(p::OCVDSolution,::Val{:V_on}) = begin r=filter(row -> row.t <= 0u"s", p.df); r.t .-=r.t[1]; r.V end
getproperty(p::OCVDSolution,::Val{:V_decay}) = filter(row -> row.t >= 0u"s", p.df).V
getproperty(p::OCVDSolution,::Val{:sol_on}) =begin r=filter(row -> row.t <= 0u"s", p.df); r.t .-=r.t[1]; r end
getproperty(p::OCVDSolution,::Val{:sol_decay}) = filter(row -> row.t >= 0u"s", p.df)


# OCVDSolution(on, decay, p::OCVDProblem) = OCVDSolution(
#     on,
#     decay,
#     p,
#     get_t(on) ,
#     get_t(decay),
#     get_V(on),
#     get_V(decay),
# )

function solve(pp::OCVDProblem, args...)
    p=deepcopy(pp)
    tend = p.on_time
    parms= p.parameters

    #  init_c = Cell(parms;mode = :oc,alg_ctl = p.alg_control)

    # @debug p.on_time
    #s1 = solve(init_c,tend=p.on_time)

    #@debug (get_V(s1))[1] (get_V(s1))[end]

    #t0=ustrip(upreferred(copy(p.on_time)))
    #p2 = setproperties(
    #    p.parameters,
    #    light  = t -> p.parameters.light(t+t0),
    #)

    alg_ctl =  setproperties(
            p.alg_control,
            tend= p.decay_time,
            tstart= -p.on_time
        )
    c = Cell(parms;mode = :oc,alg_ctl =alg_ctl)
    sol = solve(c)#.u[end]
    return OCVDSolution(todf(sol), p)

end
