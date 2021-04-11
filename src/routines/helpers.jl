"""
    (LinearAlgebra.Tridiagonal(N::Integer, a::T, b::T, c::T) where T <: Union{Number, AbstractArray{L, 1}}) where L <: Number

A simple wrapper for creating tridiagonal matrices a bit more convineantly

# Arguments:
- `N`: Size of the NxN tridiagonal
- `a`: row a
- `b`: row b
- `c`: row c

# Example
```jldoctest
julia> AnnA.Tridiagonal(4,1,-2,-1)
4×4 LinearAlgebra.Tridiagonal{Int64, Vector{Int64}}:
 -2  -1   ⋅   ⋅
  1  -2  -1   ⋅
  ⋅   1  -2  -1
  ⋅   ⋅   1  -2
```
"""
function LinearAlgebra.Tridiagonal(N::Integer,a::T,b::T,c::T) where T<:Union{Number,AbstractArray{L,1}} where L<:Number
    da  = a*ones(T,N-1)
    db  = b*ones(T,N)
    dc  = c*ones(T,N-1)

    LinearAlgebra.Tridiagonal(da,db,dc)
end


"""
    struct Pulse{A, AA, T, TT, I} <: Function

Holds all information about the pulse formation.

# Arguments:
- `Δh`: Pulse Amplitude
- `h₀`: Baseline
- `w`: Pulse width
- `Δt`: rise and fall time
- `tₑ`: end time (end of fall)
- `ts`: array of timepoints
- `hs`: array of amplitude points
- `itp`: interploating function
"""
struct Pulse{A,AA,T,TT,I} <: Function
    Δh::A   # pulse hight
    h₀::A   # base line
    w::T    # pulse width
    Δt::T   # rise/fall time
    tₑ::T   # endtime (end of fall)
    ts::TT  # time point array
    hs::AA  # ampl point array
    itp::I  # interpolation function
end
"""
    (p::Pulse)(t)

calls the interpolating function from `Pulse`
"""
function (p::Pulse)(t)
    p.itp(t)
end

"""
    pulse(; Δh=1.0, h₀=0.0, w=1.0 - 2.0e-12, Δt=1.0e-12, tₑ=1.0)

creates a pulse. Monotonic interpolation is used, so the pulse shape is
differentiable, allowing for autodiff.

# Arguments:
- `Δh`: Pulse Amplitude
- `h₀`: Baseline
- `w`: Pulse width
- `Δt`: rise and fall time
- `tₑ`: end time (end of fall)
"""
function pulse(;Δh = 1.0, h₀ = 0.0, w = (1.0-2e-12), Δt = 1.0e-12, tₑ = 1.0)
    lim=1.7976931348623157e+308
    ts=[-lim,0, Δt,Δt+w,2*Δt+w,lim].-(1*Δt+w-tₑ)
    hs=[h₀ ,h₀,h₀+Δh,h₀+Δh,h₀,h₀]
    itp = interpolate(ts, hs ,SteffenMonotonicInterpolation())
    Pulse(Δh, h₀, w, Δt, tₑ, ts, hs, itp)
end

"""
    unpac_struct(s...)

Returns a `OrderedDict` containing all `field => value` pairs. Allows for multiple structs as input. The result is then a single Dict. The letter arguments will overwrite previous one is existing already
"""
function unpac_struct(s...)

    p=OrderedDict{Symbol,Any}()
    for arg in s
        prm_names = propertynames(arg,true)
        for n in prm_names
            #@eval $n = $p.$(n)
            p[n] = getproperty(arg,n)
        end
    end
    return p
end

"""
    trianglewave(φ)

2π-periodic triangle wave of φ with amplitude 1.
"""
function trianglewave(x)
    xm = mod2pi(x + π/2)
    xm < π ? 2xm/π - 1 : -2xm/π + 3
end