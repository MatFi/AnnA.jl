"""
    struct NodimParameters{T, F <: Function, Fl <: Function, Fr <: Function, Fg <: Function, Fv <: Function}

holds all nondimensionalized parameters
# Arguments:
- `λ`: DESCRIPTION
- `λ²`: DESCRIPTION
- `δ`: DESCRIPTION
- `χ`: DESCRIPTION
- `σ`: DESCRIPTION
- `κₚ`: DESCRIPTION
- `κₙ`: DESCRIPTION
- `αb`: DESCRIPTION
- `dpt`: DESCRIPTION
- `dpf`: DESCRIPTION
- `Vbi`: Vbi optained from bolzmann statistics
- `wₑ`: DESCRIPTION
- `wₕ`: DESCRIPTION
- `κₑ`: DESCRIPTION
- `κₕ`: DESCRIPTION
- `rₑ`: DESCRIPTION
- `rₕ`: DESCRIPTION
- `λₑ`: DESCRIPTION
- `λₑ²`: DESCRIPTION
- `λₕ²`: DESCRIPTION
- `λₕ`: DESCRIPTION

- `kₑ`: proportionality factor between ETL / PVSK electron densities
- `kₕ`: proportionality factor between PVSK / HTL hole densities

- `R`: Bulk recombination function
- `Rₗ`: ETL/PVSK surface recombnation fnction
- `Rᵣ`: PVSK/HTL surface recombination function
- `G`: Carrier generation function
- `V`: Voltage/kt
"""
struct NodimParameters{T,F<:Function,Fl<:Function,Fr<:Function, Fg<:Function, Fv<:Function, σ, Tf1,Tf2,Tf3,Tf4,Tf5,Tf6,Tf7,Tf8,Tf9}

    λ::Tf1     #Debye length
    λ²::Tf2     #Debye length square
    δ::T      #ratio of typical electron/ion densities
    χ::T      #ratio of typical hole/electron densities
    ϰ::T  # ratio of perovskite doping and ionic concentration
    σ::T      #ratio of carrier and ionic timescales
    κₚ::T     #hole current Parameter
    κₙ::T     #electron current parameter
    αb::T     #Beer-Lambert Law parameter
    dpt::Tf3    # Dislacement current density factor
    dpf::T    # ionic flux displacement factor

    Vbi::T     #Built in potential
    """transport layer params"""
    wₑ::T     # relative width of ETL
    wₕ::T     # relative width of HTL
    κₑ::T     # ETL current parameter
    κₕ::T     # HTL current parameter
    rₑ::Tf4     # relative ETL permitivity
    rₕ::Tf5     # relaitve HTL permitivity
    λₑ::Tf6     # relative ETL Debye length
    λₑ²::Tf7    # relative ETL Debye length square
    λₕ²::Tf8    # relaitve HTL Debye length square
    λₕ::Tf9     # relaitve HTL Debye Length
    nᵢ²::T     #nodim intrinsc conc
    σₛₕ::T # nundim Shunt conductivity
    σₛ::σ  # nundim series conductivity
    """Interface parameters"""
    kₑ::T     # ratio between electron densities ETL/perovskite interface
    kₕ::T     # ratio between hole densities perovkite/HTL interface

    """Recombination parameters"""
    R::F    # Bulk Recombination rate
    Rₗ::Fl  # ETL/perovskite interface recombination rate
    Rᵣ::Fr   # perovskite/HTL interface recombination rate

    G::Fg   # Generation rate function
    V::Fv   # Voltage function

end

"""
    NonDimensionalize(p::AbstractParameters)

Performs the nondimensionalisation of the input parameters and returns a `NodimParameters` struct
"""
function NonDimensionalize(parm::AbstractParameters)
    prm = deepcopy(parm)
    #Unpack parameters
    p=unpac_struct(prm)

    #pack filednames to dict
    d=OrderedDict{Symbol,Any}()
    for key in fieldnames(NodimParameters)
        d[key]=0.0
    end
    #recalc


    LD= t-> sqrt(p[:kB]*p[:T]*p[:εₚ](t)/(p[:q]^2*p[:N₀]))
    d[:λ]   = t-> LD(t)/p[:b]
    d[:λ²]  = t->d[:λ](t)^2
    #d[:nᵢ²] = nᵢ^2/(dₑ*dₕ)
    d[:δ]   = p[:dₑ]/p[:N₀]
    d[:χ]   = p[:dₕ]/p[:dₑ]
    d[:ϰ]   = p[:dₚ]/p[:N₀] 
    d[:σ]   = p[:dₑ]/(p[:G₀]*p[:τᵢ])
    d[:κₚ]  = p[:Dₚ]*p[:dₕ]/(p[:G₀]*p[:b]^2)
    d[:κₙ]  = p[:Dₙ]*p[:dₑ]/(p[:G₀]*p[:b]^2)
    d[:αb]  = p[:α]*p[:b]
    d[:dpt] = t-> p[:εₚ](t)*uconvert(Unitful.NoUnits,p[:VT]/(p[:q]*p[:G₀]*p[:b]^2*p[:τᵢ]))
    d[:dpf] = p[:Dᵢ]*p[:N₀]/(p[:G₀]*p[:b]^2)
    d[:Vbi] = p[:Vbi]/p[:VT]/p[:q]
    d[:wₑ]  = p[:bₑ]/p[:b]
    d[:wₕ]  = p[:bₕ]/p[:b]
    d[:κₑ]  = p[:Dₑ]*d[:κₙ]/p[:Dₙ]
    d[:κₕ]  = p[:Dₕ]*d[:κₚ]/p[:Dₚ]
    d[:rₑ]  = t->p[:εₑ]/p[:εₚ](t)
    d[:rₕ]  = t->p[:εₕ]/p[:εₚ](t)
    d[:λₑ²] = t-> d[:rₑ](t)*p[:N₀]/p[:dₑ]*d[:λ²](t)
    d[:λₑ]  = t-> sqrt(d[:λₑ²](t))
    d[:λₕ²] = t->d[:rₕ](t)*p[:N₀]/p[:dₕ]*d[:λ²](t)
    d[:λₕ]  = t->sqrt(d[:λₕ²](t))

    d[:kₑ]  = p[:gc]/p[:gcₑ]*exp((p[:Ecₑ]-p[:Ec])/(p[:kB]*p[:T]))
    d[:kₕ]  = p[:gv]/p[:gvₕ]*exp((p[:Ev]-p[:Evₕ])/(p[:kB]*p[:T]))

    d[:σₛₕ] = 1/p[:Rₛₕ]/p[:jay]*p[:VT]

    ti=ustrip(p[:τᵢ]|> u"s")
 
    if p[:Rₛ] isa Function
        d[:σₛ] = t-> 1/p[:Rₛ](t*ti)/p[:jay]*p[:VT]
    else
        d[:σₛ] = t-> 1/p[:Rₛ]/p[:jay]*p[:VT]
    end
    d[:σₛ] = σ(p[:Rₛ],parm)
    d[:nᵢ²] = uconvert(Unitful.NoUnits,p[:nᵢ]^2/(p[:dₑ]*p[:dₕ]))
    #Test for nodimensionalty
    for i in eachindex(d)
        if typeof(d[i]) <: Unitful.Quantity
            try
                d[i]= uconvert(Unitful.NoUnits,d[i])
            catch e
                error("Nondimensionalisation failed. The Perameter $i has a dimension left")
            end
        end
    end

    nᵢ² = d[:nᵢ²] 
    """Bulk Recombination"""
    if !iszero(p[:τₚ])# && τₙ>0u"s"

        kk  = uconvert(Unitful.NoUnits,p[:k₂]*p[:dₑ]*p[:dₕ]/p[:G₀])
        γ   = uconvert(Unitful.NoUnits,p[:dₕ]/(p[:τₚ]*p[:G₀]))
        τᵣ  = uconvert(Unitful.NoUnits,p[:τₙ]*p[:dₕ]/(p[:τₚ]*p[:dₑ]))
        rtrap = uconvert(Unitful.NoUnits,(p[:τₙ]+p[:τₚ])*p[:nᵢ]/(p[:τₚ]*p[:dₑ]))

        d[:R]= Rec_function(promote(copy(nᵢ²),copy(kk),copy(γ),copy(τᵣ),copy(rtrap))...)
    #    d[:R]= R!(R,n,p) -> @. R = (n*p-nᵢ²)*(kk+γ/(n+τᵣ*p+rtrap))
    else
        kk  = uconvert(Unitful.NoUnits,p[:k₂]*p[:dₑ]*p[:dₕ]/(p[:G₀]))
        γ   = 0
        τᵣ  = 1
        rtrap = 1
        d[:R]= Rec_function(promote(copy(nᵢ²),copy(kk),copy(γ),copy(τᵣ),copy(rtrap))...)
    end

    """Interface ETM recombination"""
    if !iszero(p[:vₙₑ]) #&& vₙₑ>0u"m/s
        kk  = uconvert(Unitful.NoUnits,p[:k₂ₑ]*p[:dₑ]*p[:dₕ]/(p[:b]*p[:G₀]))
        γ   = uconvert(Unitful.NoUnits,p[:dₕ]*p[:vₚₑ]/(p[:b]*p[:G₀]))
        τᵣ  = uconvert(Unitful.NoUnits,p[:dₕ]*p[:vₚₑ]/(p[:vₙₑ]*p[:dₑ]))
        rtrap = uconvert(Unitful.NoUnits,(1/d[:kₑ]+p[:vₚₑ]/p[:vₙₑ])*p[:nᵢ]/p[:dₑ])
        d[:Rₗ]= Rec_function(promote(copy(nᵢ²),copy(kk),copy(γ),copy(τᵣ),copy(rtrap))...)
    else
        kk  = uconvert(Unitful.NoUnits,p[:k₂ₑ]*p[:dₑ]*p[:dₕ]/(p[:b]*p[:G₀]))
        γ   = 0
        τᵣ  = 1
        rtrap = 1
        d[:Rₗ]= Rec_function(promote(copy(nᵢ²),copy(kk),copy(γ),copy(τᵣ),copy(rtrap))...)
    end

    """Interface HTM recombination"""
    if !iszero(p[:vₚₕ]) # && vₙₕ>0u"m/s"
        kk  = uconvert(Unitful.NoUnits,p[:k₂ₕ]*p[:dₑ]*p[:dₕ]/(p[:b]*p[:G₀]))
        γ   = uconvert(Unitful.NoUnits,p[:dₑ]*p[:vₙₕ]/(p[:b]*p[:G₀]))
        τᵣ  = uconvert(Unitful.NoUnits,p[:dₑ]*p[:vₙₕ]/(p[:vₚₕ]*p[:dₕ]))
        rtrap = uconvert(Unitful.NoUnits,(1/d[:kₕ]+p[:vₙₕ]/p[:vₚₕ])*p[:nᵢ]/p[:dₕ])
        d[:Rᵣ]= Rec_function(promote(copy(nᵢ²),copy(kk),copy(γ),copy(τᵣ),copy(rtrap))...)
    else
        kk  = uconvert(Unitful.NoUnits,p[:k₂ₕ]*p[:dₑ]*p[:dₕ]/(p[:b]*p[:G₀]))
        γ   = 0
        τᵣ  = 1
        rtrap = 1
        d[:Rᵣ]= Rec_function(promote(copy(nᵢ²),copy(kk),copy(γ),copy(τᵣ),copy(rtrap))...)
    end

#    d[:G] = Gen_function(uconvert(Unitful.NoUnits,p[:α]*p[:b]), Float64(p[:dir]), t -> p[:light](t*ustrip(p[:τᵢ]|> u"s")), uconvert(Unitful.NoUnits,p[:τᵢ]/1u"s"))
    
    d[:G] = Gen_function(uconvert(Unitful.NoUnits,p[:α]*p[:b]), Float64(p[:dir]), t -> prm.light(t*ti), uconvert(Unitful.NoUnits,p[:τᵢ]/1u"s"))
    
    d[:V] = Pot_function(uconvert(Unitful.NoUnits,p[:VT]/1u"V"),p[:V],uconvert(Unitful.NoUnits,p[:τᵢ]/1u"s"))

    #when the ions are frozen no corrent is calculated afterwards
    d[:dpf] = prm.freeze_ions ? Float64(0) : d[:dpf]
    NodimParameters(collect(values(d))...)
end

"""
    struct Rec_function{T} <: Function

after constructing this type using `r = Rec_function(nᵢ², kk, γ, τᵣ, rtrap)` givs  access to the recombination function `r(R, n, p)`

#Arguments:
- `nᵢ²`: nondimensionalized intrinsic carrier square
- `kk`: nondimensionalized bimolecular recombination rate
- `γ`: nondimensionalized rate constant for hole-dominated recombination
- `τᵣ`: ratio of SRH carrier lifetime
- `rtrap`: nondimensionalized k₂ (deep trap) parameter
"""
struct Rec_function{T} <: Function
    nᵢ²::T
    kk::T
    γ::T
    τᵣ::T
    rtrap::T
end

"""
    (r!::Rec_function)(R::AbstractArray, n::AbstractArray, p::AbstractArray)

Call the inplace recombination function. This is where the bulk recombination
is calculated

#Arguments:
- `R`: Buffer Array
- `n`: Electron Array
- `p`: Hole Array
"""
function (r!::Rec_function)(R::AbstractArray, n::AbstractArray, p::AbstractArray)
        @. R = (n*p-r!.nᵢ²)*(r!.kk+r!.γ/(n+r!.τᵣ*p+r!.rtrap))
end

"""
    (r::Rec_function)(n::T, p::T) where T

Out of place definition of the recombination function. Calculates the
recombination for a given n and p. This is intended to calculate the sruface
crecombinations
"""
function (r::Rec_function)(n::Number, p::Number)
        return (n*p-r.nᵢ²)*(r.kk+r.γ/(n+r.τᵣ*p+r.rtrap))
end

"""
    struct Gen_function{T, F, S} <: Function

after constructing this type using `g! = Gen_function(αb, dir, light, tion)` givs
access to the generation function `g!(R, n, p)`

#Arguments:
- `αb`: Nondim absorption coefficient
- `dir`: Light direction +1/-1
- `light`: scalar tempoal function of generation
- `tion`: ionic timescale, needed for redimensionalize the time
"""
struct Gen_function{T,F,S} <: Function
    αb::T
    dir::T
    light::F
    tion::S
end

"""
    (g!::Gen_function)(G, x, t)

Inplace generation function

#Arguments:
- `G`: Buffer Array
- `x`: Nodimensionalized positions Array
- `t`: Nondimensionalized time
"""
function (g!::Gen_function)(G,x,t)
    #put the light function in rhs solves a forward diff problem
    @.  G   =  g!.αb/(1-exp(-g!.αb))*exp(-g!.αb*(g!.dir*x + (1-g!.dir)/2));
#    @.  G   =  (g!.light(t*g!.tion))*g!.αb/(1-exp(-g!.αb))*exp(-g!.αb*(g!.dir*x + (1-g!.dir)/2));
end


"""
    struct Pot_function{T, F} <: Function

after constructing this type using `v =  Pot_function(VT, V, tion)` givs  access to the generation function `v(t)`

#Arguments:
- `VT`: thermal voltage
- `V`: temporal potential function
- `tion`: ionic timescale, needed for redimensionalize the time
"""
struct Pot_function{T,F} <: Function
    VT::T
    V::F
    tion::T
end

"""
    (v::Pot_function)(t)

Returns the nondimensionalized potential at time `t`
"""
function (v::Pot_function)(t)
    v.V(t*v.tion)/v.VT
end

struct TVariable end
struct TConstant end
struct σ{FTag,T}
    s::T
end
"""
    σ(R,p::AbstractParameters)

Wrapper for nondimensionalization of timedependant conductivity


# Example
```jldoctest
julia> p = Parameters(Rₛ = 1.0e-20 * u"V/A*mm^2",);

julia> s = AnnA.σ(p.Rₛ,p)
AnnA.σ{AnnA.TConstant, Float64}(1.15893217665509e22)

julia> s.(0:4)
5-element Vector{Float64}:
 1.15893217665509e22
 1.15893217665509e22
 1.15893217665509e22
 1.15893217665509e22
 1.15893217665509e22

julia> pt = Parameters(Rₛ = t-> cos(t)^2*1.0e-20 * u"V/A*mm^2",);

julia> st = AnnA.σ(pt.Rₛ,pt)
AnnA.σ{AnnA.TVariable, Function}(AnnA.var"#22#23"{var"#1#2", Parameters, Float64}(var"#1#2"(), Parameters(400, 8.854187817e-12 F m^-1, 1.6021766209e-19 C, 8.61733035e-5 eV K^-1, 4.0e-7 m, 24.1, -3.7 eV, -5.3 eV, 0 m^-3, 0.00017 m^2 s^-1, 0.00017 m^2 s^-1, 0.2, 0.2, 1.6e24 m^-3, 6.5e-8 m^2 s^-1, 0.58 eV, false, 300 K, 1.3e7 m^-1, 1.4e21 m^-2 s^-1, 1, AnnA.Pulse{Float64, Vector{Float64}, Float64, Vector{Float64}, Interpolations.MonotonicInterpolation{Float64, Float64, Float64, Float64, Float64, Interpolations.SteffenMonotonicInterpolation, Vector{Float64}, Vector{Float64}}}(1.0, 0.0, 2.0, 1.0e-12, 1.0, [-1.0e308, -1.000000000001, -1.0, 1.0, 1.000000000001, 1.0e308], [0.0, 0.0, 1.0, 1.0, 0.0, 0.0], 6-element interpolate(::Vector{Float64}, ::Vector{Float64}, Interpolations.SteffenMonotonicInterpolation()) with element type Float64:
  0.0
  0.0
  1.0000000000000002
  1.0
 -2.220446049250313e-16
  0.0), AnnA.var"#3#5"(), 1.0e6 m^2 V A^-1, var"#1#2"(), 3.0e-7 s, 3.0e-7 s, 3.22e-17 m^3 s^-1, 0 m^4 s^-1, 0 m^4 s^-1, 0 m s^-1, 0 m s^-1, 0 m s^-1, 0 m s^-1, 1.0e18 cm^-3, 1.5, -4.0 eV, 1.0e-7 m, 3, 1.0e-7 m^2 s^-1, 1.0e18 cm^-3, 12, -5 eV, 1.0e-7 m, 3, 1.0e-7 m^2 s^-1), 158.17425072938966))

julia> st.(0:4)
5-element Vector{Float64}:
 1.15893217665509e22
 5.5156335151157665e22
 3.4478995483303833e22
 1.182707158333871e22
 1.0789052098079947e23

```
"""
function σ(R,p::AbstractParameters)
    if R isa Function
        ti=ustrip(p.τᵢ|> u"s")
        σ{TVariable,Function}(t-> uconvert(Unitful.NoUnits,1/R(t*ti)/p.jay*p.VT))
    else
        s=uconvert(Unitful.NoUnits,1/R/p.jay*p.VT)
        σ{TConstant,typeof(s)}(s)
    end
end
function (s::σ{TVariable,T})(t) where T
    s.s(t)
end
function (s::σ{TConstant,T})(t) where T
    s.s
end


