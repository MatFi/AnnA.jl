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
struct NodimParameters{T,F<:Function,Fl<:Function,Fr<:Function, Fg<:Function, Fv<:Function}

    λ::T     #Debye length
    λ²::T     #Debye length square
    δ::T      #ratio of typical electron/ion densities
    χ::T      #ratio of typical hole/electron densities
    ϰ::T  # ratio of perovskite doping and ionic concentration
    σ::T      #ratio of carrier and ionic timescales
    κₚ::T     #hole current Parameter
    κₙ::T     #electron current parameter
    αb::T     #Beer-Lambert Law parameter
    dpt::T    # Dislacement current density factor
    dpf::T    # ionic flux displacement factor

    Vbi::T     #Built in potential
    """transport layer params"""
    wₑ::T     # relative width of ETL
    wₕ::T     # relative width of HTL
    κₑ::T     # ETL current parameter
    κₕ::T     # HTL current parameter
    rₑ::T     # relative ETL permitivity
    rₕ::T     # relaitve HTL permitivity
    λₑ::T     # relative ETL Debye length
    λₑ²::T    # relative ETL Debye length square
    λₕ²::T    # relaitve HTL Debye length square
    λₕ::T     # relaitve HTL Debye Length
    nᵢ²::T     #nodim intrinsc conc
    σₛₕ::T # nundim Shunt conductivity

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
function NonDimensionalize(prm::AbstractParameters)

    #Unpack parameters
    p=unpac_struct(prm)

    #pack filednames to dict
    d=OrderedDict{Symbol,Any}()
    for key in fieldnames(NodimParameters)
        d[key]=0.0
    end
    #recalc


    LD= sqrt(p[:kB]*p[:T]*p[:εₚ]/(p[:q]^2*p[:N₀]))
    d[:λ]   = LD/p[:b]
    d[:λ²]  = d[:λ]^2
    #d[:nᵢ²] = nᵢ^2/(dₑ*dₕ)
    d[:δ]   = p[:dₑ]/p[:N₀]
    d[:χ]   = p[:dₕ]/p[:dₑ]
    d[:ϰ]   = p[:dₚ]/p[:N₀] 
    d[:σ]   = p[:dₑ]/(p[:G₀]*p[:τᵢ])
    d[:κₚ]  = p[:Dₚ]*p[:dₕ]/(p[:G₀]*p[:b]^2)
    d[:κₙ]  = p[:Dₙ]*p[:dₑ]/(p[:G₀]*p[:b]^2)
    d[:αb]  = p[:α]*p[:b]
    d[:dpt] = p[:εₚ]*p[:VT]/(p[:q]*p[:G₀]*p[:b]^2*p[:τᵢ])
    d[:dpf] = p[:Dᵢ]*p[:N₀]/(p[:G₀]*p[:b]^2)
    d[:Vbi] = p[:Vbi]/p[:VT]/p[:q]
    d[:wₑ]  = p[:bₑ]/p[:b]
    d[:wₕ]  = p[:bₕ]/p[:b]
    d[:κₑ]  = p[:Dₑ]*d[:κₙ]/p[:Dₙ]
    d[:κₕ]  = p[:Dₕ]*d[:κₚ]/p[:Dₚ]
    d[:rₑ]  = p[:εₑ]/p[:εₚ]
    d[:rₕ]  = p[:εₕ]/p[:εₚ]
    d[:λₑ²] = d[:rₑ]*p[:N₀]/p[:dₑ]*d[:λ²]
    d[:λₑ]  = sqrt(d[:λₑ²])
    d[:λₕ²] = d[:rₕ]*p[:N₀]/p[:dₕ]*d[:λ²]
    d[:λₕ]  = sqrt(d[:λₕ²])

    d[:kₑ]  = p[:gc]/p[:gcₑ]*exp((p[:Ecₑ]-p[:Ec])/(p[:kB]*p[:T]))
    d[:kₕ]  = p[:gv]/p[:gvₕ]*exp((p[:Ev]-p[:Evₕ])/(p[:kB]*p[:T]))

    d[:σₛₕ] = 1/p[:Rₛₕ]/p[:jay]*p[:VT]
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

    d[:G] = Gen_function(uconvert(Unitful.NoUnits,p[:α]*p[:b]), Float64(p[:dir]), t -> p[:light](t*ustrip(p[:τᵢ]|> u"s")), uconvert(Unitful.NoUnits,p[:τᵢ]/1u"s"))
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
