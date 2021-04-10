abstract type AbstractParameters end
const vPar = Union{Number,Tuple,AbstractArray,Function,Array{Function,1},T} where T

Base.@kwdef mutable struct Parameters <: AbstractParameters
    """Model Parameters"""
    N::Integer  =   400     # Subintervals in perovskite layer, resulting in N+1 Grid points
    
    # Physical parameters
    ε₀::vPar = 8.854187817e-12u"F/m"# Permitivity of free space
    q::vPar  = 1.6021766209e-19u"C" # Elemntary charge of a proton
    kB::vPar = 8.61733035e-5u"eV/K" # Bolzmann konstant
   
    
    # Perovskite parameters
    b::vPar  = 400e-9u"m"           # Perovskite Layer thickness
    ε::vPar = 24.1                  # Perovskite permitivity
    Ec::vPar = -3.7u"eV"            # Perovskite Conduction band energy
    Ev::vPar = -5.3u"eV"            # Perovskite Valence band energy
    dₚ   = 0u"m^-3"                 # Peroviskite doping concentration
    Dₙ::vPar = 1.7e-4u"m^2/s"       # Perovskite electron diffusion coefficient
    Dₚ::vPar = 1.7e-4u"m^2/s"       # Perovskite hole diffusion coefficient
    mₑ::vPar = 0.2                  # Perovskite electron mass
    mₕ::vPar = 0.2                  # Perovskite hole mass
    
    # Ion Parameters
    N₀::vPar = 1.6e24u"m^-3"        # Typical density of ion vacancys
    Dᵢ₀::vPar = 6.5e-8u"m^2/s"      # Diffusion constant
    Eᵢₐ::vPar = 0.58 * u"eV"        # Ativation energy of vacancy diffusion
    freeze_ions = false
    
    # Environment Parameters
    T::vPar  = 300u"K"              # Temperature
    α::vPar  = 1.3e7u"1/m"          # Perovskite absorption koefficient
    Fₚₕ::vPar = 1.4e21u"m^-2*s^-1"  # Photonflux (bandgap dependent)
    dir::Integer =   1              # Shine light trough  1 -> ETL, -1 -> HTL
    light::vPar = pulse(tₑ=1.0)
    V::vPar = t -> 0                 # Voltage protocoll
    Rₛₕ::vPar = 1e3u"V/A*m^2"       # Shunt resistance

    # Recombination Parameters
    τₙ::vPar = 3e-7u"s"             # electron pseudo lifetime
    τₚ::vPar = 3e-7u"s"             # hole pseudo lifetime
    k₂::vPar = 3.22e-17u"m^3/s"     # second order rate constant
   
    # Interface Recombination
    k₂ₑ::vPar = 0u"m^4/s"           # ETL/perovskite bimolecular recombination rate
    k₂ₕ::vPar = 0u"m^4/s"           # perovskite/HTL bimolecular recombination rate
    vₙₑ::vPar = 0u"m/s"             # electron recombination velocity for SHR/ETL
    vₚₑ::vPar = 0u"m/s"             # hole recombination velocity for SHR/ETL
    vₙₕ::vPar = 0u"m/s"             # electron recombination velocity for SHR/HTL
    vₚₕ::vPar = 0u"m/s"             # hole recombination velocity for SHR/HTL
   
    # ELT Parameters
    dₑ::vPar = 1e18u"cm^-3"          # ETL effective doping density
    mcₑ::vPar = 1.5                  # ETL electron mass
    Ecₑ::vPar = -4.0 * u"eV"           # ETL conduction band energy
    bₑ::vPar = 100e-9u"m"           # ETL width
    εₑᵣ::Real = 3                   # ETL permitivity
    Dₑ::vPar = 1e-7u"m^2/s"         # ETL electron diffusion coeficcient
   
    # HTL Parameters
    dₕ::vPar = 1e18u"cm^-3"          # HTL effective doping density
    mvₕ::vPar = 12                   # HTL hole mass
    Evₕ::vPar = -5 * u"eV"             # HTL valence band energy
    bₕ::vPar = 100e-9u"m"           # HTL width
    εₕᵣ::Real = 3                   # HTL permitivity
    Dₕ::vPar = 1e-7u"m^2/s"         # HTL electron diffusion coeficcient
end

propertynames(p::AbstractParameters,private::Bool=true) = begin
    pubs = fieldnames(typeof(p))
    privs = (:gc, :gv, :gcₑ, :gvₕ, :εₚ, :εₕ, :εₑ, :Eg, :Dᵢ, :nᵢ, :Efₑ, :Efₕ, :Vbi, :VT, :τᵢ, :G₀, :jay, :n0, :p0)
    private && return (pubs..., privs...)
    (pubs...,)
end

setproperty!(p::AbstractParameters,n::Symbol,x) = setproperty!(p::AbstractParameters, Val{n}(), x)
setproperty!(p::AbstractParameters,::Val{S},x) where {S} = setfield!(p, S, x)
setproperty!(p::AbstractParameters, ::Val{:nᵢ},x) = begin
    @warn "nᵢ is transformed to gc and gv via gv=gc = sqrt(nᵢ²exp(Eg/kT))"
    g = sqrt(x^2 * exp(uconvert(Unitful.NoUnits, p.Eg / (p.kB * p.T))));
    setfield!(p, :gc, g);
    setfield!(p, :gv, g);
end
setproperty!(p::AbstractParameters, ::Val{:gc},x) = setfield!(p, :mₑ, (x / 2)^(2 / 3) * (4.13566769692e-15u"eV*s")^2 / (2π * 9.10938356e-31u"kg" * p.kB * p.T) |> upreferred)
setproperty!(p::AbstractParameters, ::Val{:gv},x) = setfield!(p, :mₕ, (x / 2)^(2 / 3) * (4.13566769692e-15u"eV*s")^2 / (2π * 9.10938356e-31u"kg" * p.kB * p.T) |> upreferred)
setproperty!(p::AbstractParameters, ::Val{:gcₑ},x) = setfield!(p, :mcₑ, (x / 2)^(2 / 3) * (4.13566769692e-15u"eV*s")^2 / (2π * 9.10938356e-31u"kg" * p.kB * p.T) |> upreferred)
setproperty!(p::AbstractParameters, ::Val{:gvₕ},x) = setfield!(p, :mvₕ, (x / 2)^(2 / 3) * (4.13566769692e-15u"eV*s")^2 / (2π * 9.10938356e-31u"kg" * p.kB * p.T) |> upreferred)

getproperty(p::AbstractParameters,n::Symbol) = getproperty(p::AbstractParameters, Val{n}())
getproperty(p::AbstractParameters,::Val{S}) where {S} = getfield(p, S)
getproperty(p::AbstractParameters,::Val{:εₚ}) = p.ε * p.ε₀
getproperty(p::AbstractParameters,::Val{:εₕ}) = p.εₕᵣ * p.ε₀
getproperty(p::AbstractParameters,::Val{:εₑ}) = p.εₑᵣ * p.ε₀
getproperty(p::AbstractParameters,::Val{:Eg}) = p.Ec - p.Ev |> u"eV"
getproperty(p::AbstractParameters,::Val{:Dᵢ}) = p.Dᵢ₀ * exp(-p.Eᵢₐ / (p.kB * p.T)) |> u"m^2/s"
getproperty(p::AbstractParameters,::Val{:nᵢ}) = sqrt(p.gc * p.gv * exp(-p.Eg / (p.kB * p.T)))
getproperty(p::AbstractParameters,::Val{:Efₑ}) = p.Ecₑ - p.kB * p.T * log(p.gcₑ / p.dₑ) |> u"eV"
getproperty(p::AbstractParameters,::Val{:Efₕ}) = p.Evₕ + p.kB * p.T * log(p.gvₕ / p.dₕ) |> u"eV"
getproperty(p::AbstractParameters,::Val{:Vbi}) = (p.Efₑ - p.Efₕ) |> u"eV"
getproperty(p::AbstractParameters,::Val{:VT}) = p.kB * p.T / p.q  |> u"V"
getproperty(p::AbstractParameters,::Val{:τᵢ}) = p.b / p.Dᵢ * sqrt(p.kB * p.T * p.εₚ / p.N₀ / p.q^2) |> u"s"
getproperty(p::AbstractParameters,::Val{:G₀}) = p.Fₚₕ / p.b * (1 - exp(-p.α * p.b)) |> upreferred
getproperty(p::AbstractParameters,::Val{:jay}) = p.q * p.G₀ * p.b |> u"A/m^2"
getproperty(p::AbstractParameters,::Val{:n0}) = p.dₑ * p.gc / p.gcₑ * exp((p.Ecₑ - p.Ec) / (p.kB * p.T))  |> u"m^-3"
getproperty(p::AbstractParameters,::Val{:p0}) = p.dₕ * p.gv / p.gvₕ * exp((p.Ev - p.Evₕ) / (p.kB * p.T)) |> u"m^-3"
getproperty(p::AbstractParameters,::Val{:gc}) = 2 * (2π * p.mₑ * 9.10938356e-31u"kg" * p.kB * p.T / (4.13566769692e-15u"eV*s")^2)^(3 / 2) |> u"m^-3"
getproperty(p::AbstractParameters,::Val{:gv}) = 2 * (2π * p.mₕ * 9.10938356e-31u"kg" * p.kB * p.T / (4.13566769692e-15u"eV*s")^2)^(3 / 2) |> u"m^-3"
getproperty(p::AbstractParameters,::Val{:gcₑ}) = 2 * (2π * p.mcₑ * 9.10938356e-31u"kg" * p.kB * p.T / (4.13566769692e-15u"eV*s")^2)^(3 / 2) |> u"m^-3"
getproperty(p::AbstractParameters,::Val{:gvₕ}) = 2 * (2π * p.mvₕ * 9.10938356e-31u"kg" * p.kB * p.T / (4.13566769692e-15u"eV*s")^2)^(3 / 2) |> u"m^-3"

Base.@kwdef mutable struct AlgControl
# TODO: make this just a Dict
    # Logging::LogLevel   = Warn    # Sets the log level:  Debug ⊂ Info ⊂ Warn ⊂ Error

    ss_tol::Float64     = 1e-6      # ftol for steady state detection
    reltol::Float64     = 1e-6      # Relative tolerance for the solvers (auto_abstol =  u*reltol)
    abstol              = 1e-5
    autodiff::Bool      = true      # Enable autodifferentiation (highly suggested)
    alg     = Rodas4P(autodiff=autodiff) # Integrator Algorithm. Optional suggestions are: Rodas5(), ROS34PW2(), ROS3P(), QBDF(), QNDF()
    dtmin   = 1e-15
    dt      = 1e-14
    force_dtmin = true
    progress = false
    progress_steps = 50
    maxiters = 1000
    tstart = 0.0u"s"
    tend = 1e5u"s"
    initializealg = OrdinaryDiffEq.NoInit()

    kwargs = (; :reltol => reltol,
        :abstol => abstol,
        :dtmin => dtmin,
        :dt => dt,
        :force_dtmin => force_dtmin,
        :progress_steps => progress_steps,
        :progress => progress,
        :maxiters => maxiters,
        :initializealg => initializealg,)

end
getproperty(p::AlgControl,n::Symbol) = getproperty(p::AlgControl, Val{n}())
getproperty(p::AlgControl,::Val{S}) where {S} = getfield(p, S)
getproperty(p::AlgControl,::Val{:kwargs}) = (; :reltol => p.reltol,
    :abstol => p.abstol,
    :dtmin => p.dtmin,
    :dt => p.dt,
    :force_dtmin => p.force_dtmin,
    :progress_steps => p.progress_steps,
    :progress => p.progress,
    :maxiters => p.maxiters,
    :initializealg => p.initializealg,)
