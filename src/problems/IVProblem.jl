struct IVProblem{
    B<:Bool,
    V<:Union{AbstractArray,Tuple},
    VR<:Number,
    P<:Parameters,
    C<:AlgControl,
} <: AbstractProblem
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
    range::Union{AbstractArray,Tuple},
    rate::Unitful.AbstractQuantity;
    double_sweep = true,
    alg_control = missing,
)
    range = ustrip.(uconvert.((u"V",),range))
    rate =ustrip(uconvert(u"V/s", rate))
    tend = (maximum(range) - minimum(range)) / abs(rate) * (1+double_sweep)
    if alg_control isa Missing

        alg_control = AlgControl(
            dtmin = 1e-20,
            dt = 1e-8,
            reltol = 1e-5,
            abstol = 1e-8,
            tend = tend*u"s",
        )
    else
        #enforce correct tend
        alc_control = setproperty!(
            alg_control,
            :tend, tend*u"s"
        )
    end
    rev = diff(range)[1] <0 
    amplitude = abs(first(range) - last(range))
    freq= π/2/tend * (1+double_sweep) *(1-2*rev)
    Vt= t -> amplitude*trianglewave(freq*ustrip(t)+(rev * π ))+minimum(range)

 
   # Vt = t -> ustrip(uconvert(u"V", first(range))) +
    #          sign(ustrip(-first(range) + last(range))) *
    #          abs(ustrip(uconvert(u"V/s", rate))) * ustrip(t)

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


function solve(p::IVProblem, alg_control = p.alg_control, args...)

    #init
    τᵢ = p.parameters.τᵢ


    init_c = Cell(
        p.parameters,
        mode = :cc,
        alg_ctl = alg_control
    )
    s1 = solve(init_c)

  #  s=ProblemSolution(s1)
  #  s.df.V[1]=u"V"
  #  d = diff(s.df.V)    
    return ProblemSolution(s1)
end

function steadyStateI(V::Number,s::IVSolution;ss_tol=1e-6)
    τᵢ = s.prob.parameters.τᵢ

    #get time for V, assuming backwards and forward will end up in same state
    dyn_t =  s.sol_rwd.t #|>reverse
    dyn_v = (s.V_rwd .|> upreferred .|> ustrip)# |>reverse
    p = sortperm(dyn_v)
    dyn_t= dyn_t[p]
    dyn_v=dyn_v[p]

    itp =  LinearInterpolation(dyn_v,dyn_t)
    t_ss=itp(V |> upreferred |> ustrip)
    prob = deepcopy(s.sol_rwd.prob)
    prob_fun = deepcopy(prob.f)

#    pf =(du, u, p, t) ->
#    prob_fun = setproperties(prob.f;f=pf)
#    prob = setproperties(s.sol_rwd.prob; f=prob_fun)

    function j!(jac, x,f)
        #colors = matrix_colors(c_init.Jac)
        forwarddiff_color_jacobian!(
            jac,
            (dx, x) -> f.f(dx, x, 0, t_ss),
            x;
            colorvec = f.colorvec,
            sparsity = f.jac_prototype,
        )
        return nothing
    end

    #solfe the SS
    df = OnceDifferentiable(
        (dx, x) -> prob_fun.f(dx, x,0, t_ss),
        (jac, x) -> j!(jac, x,prob_fun),
        s.sol_rwd(t_ss),
        copy(s.sol_rwd(t_ss)),
        prob_fun.jac_prototype,
    )
    u1 = nlsolve(
        df,
        s.sol_rwd(t_ss);
        iterations = 2,
        ftol = ss_tol,
        #factor = 1e1,
        show_trace =  true, #haskey(ENV, "JULIA_DEBUG"),
        method = :trust_region#:newton,

    )

    calculate_currents(prob_fun.f, decompose(u1.zero,prob_fun.f.g), 0, decompose(u1.zero,prob_fun.f.g)) * prob_fun.f.parameters.jay

end

"""
    calc_ideality(s::IVSolution;at_V=:all)

returns the ideality of an IV solution at the unitful voltage `at_V`. If no number is provided a array is returned.
"""
function calc_ideality(s::IVSolution;at_V=:all)
    lnI= log.(s.I_rwd.|>upreferred.|>ustrip .|>abs)
    Vi = s.V_rwd

    id = diff(lnI)./diff(Vi) .* s.sol_rwd.prob.f.f.parameters.VT
    id = 1 ./ id
    if at_V isa Number
        idx = findmin(abs.((Vi[1:end-1].+Vi[2:end])./2 .- at_V))[2]
        return id[idx]
    end
    id
end

#= Have to do this with interpolations
"""
    (j::IVSolution)(v::Unitful.Voltage,dir=:reverse)

returns the current at a particular voltage, evaluated from the sweep direction dir (:reverse or :forward). Uses linear interpolation of datapoints

"""
function (j::IVSolution)(v::Unitful.Voltage,dir=:reverse)
    i=j.I_rwd
    u=j.V_rwd
    if dir == :forward
        i=j.I_fwd
        u=j.V_fwd
    end
    idx_low = argmin(abs.(u.-v))
    #we clip the endpoints for now,
    #TODO: implement extraposation / not clipping ends
    ((idx_low == 1 )||( idx_low == length(i))) && return i[idx_low]
    idx_high = abs(u[idx_low-1] - v) < abs(u[idx_low+1]-v) ? idx_low+1 :  idx_low - 1

    r = (v-u[idx_low])/(u[idx_high]-u[idx_low])
    j_low=i[idx_low]
    j_high = i[idx_high]

    return -r*(j_high-j_low)+j_low
end
=#
