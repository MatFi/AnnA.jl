struct IV end
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
        range::Union{AbstractArray,Tuple},
        rate;
        double_sweep = true,
        alg_control = AlgControl(dtmin = 1e-20,
            dt = 1e-6,
            reltol = 1e-4,
            abstol = 1e-12,
            tend = abs(sum(diff(range))) / abs(rate) * (1+double_sweep)
        ),  
    )

Creates an `IVProblem`. The voltage parameter `V` defined in `parm` will be overwritten by the linear voltage sweep defined as `V = t-> first(range) + rate*t`.
`AlgControl.tend` is forced to be `abs(sum(diff(range))) / abs(rate) * (1+double_sweep)`, an overwrite can be done on the finalized object using Setfield:
`p = Setfield.setproperties(p::IVProblem, alg_control=AlgControl(...))`. 

!!! note "Units"
	
    If no units provided for `range` and `rate`, `V` and `V/s` is assumed. 
        
"""
function IVProblem(
    parm::AbstractParameters,
    range::Union{AbstractArray,Tuple},
    rate;
    double_sweep = true,
    alg_control = missing,
)
    range = ustrip.(upreferred.(range))
    rate = ustrip(upreferred(rate))
    tend = abs(sum(diff(range))) / abs(rate) * (1+double_sweep)
    if alg_control isa Missing

        alg_control = AlgControl(
            dtmin = ustrip((1e-30 * parm.τᵢ )|>u"s"),
            dt = ustrip((1e-10 * parm.τᵢ)|>u"s"),
            reltol = 1e-5,
            abstol = 1e-8,
            force_dtmin=false,
            ss_tol=1e-5,
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
    freq= π/2/tend * (1+double_sweep)
    if rev
        Vt= t -> amplitude*(trianglewave(-freq*ustrip(t)))+minimum(range)
    else
        Vt= t -> amplitude*(trianglewave(freq*ustrip(t)))+minimum(range)
    end
   # Vt = t -> ustrip(uconvert(u"V", first(range))) +
    #          sign(ustrip(-first(range) + last(range))) *
    #          abs(ustrip(uconvert(u"V/s", rate))) * ustrip(t)

    parm = setproperties(parm, V = Vt)
    prob = IVProblem(parm, range, rate, double_sweep, alg_control)

end

function solve(p::IVProblem, alg_control = p.alg_control, args...)

    #init
    τᵢ = p.parameters.τᵢ
    init_c = Cell(
        p.parameters,
        mode = :cc,
        alg_ctl = alg_control
    )
    s1 = solve(init_c)
    s=ProblemSolution(s1,IV)
    d = diff(s.df.V).>=0u"V"    
    s.df.fwd=vcat(d[1],d)
    return s
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
