abstract type AbstractParameters end

"""
    Parameters(;       
        N  =   400,                # Subintervals in perovskite layer, 
                                   # resulting in N+1 Grid points

        # Physical parameters
        ε₀ = 8.854187817e-12u"F/m",# Permitivity of free space
        q  = 1.6021766209e-19u"C", # Elemntary charge of a proton
        kB = 8.61733035e-5u"eV/K", # Bolzmann konstant
    
        # Perovskite parameters
        b  = 400e-9u"m",           # Perovskite layer thickness
        ε = 24.1,                  # Perovskite permitivity
        Ec = -3.7u"eV",            # Perovskite conduction band energy
        Ev = -5.3u"eV",            # Perovskite valence band energy
        dₚ = 0u"m^-3",             # Perovskite doping concentration
        Dₙ = 1.7e-4u"m^2/s",       # Perovskite electron diffusion coefficient
        Dₚ = 1.7e-4u"m^2/s",       # Perovskite hole diffusion coefficient
        mₑ = 0.2,                  # Perovskite effective electron mass
        mₕ = 0.2,                  # Perovskite effective hole mass
        
        # Ion Parameters
        N₀ = 1.6e24u"m^-3",        # Typical density of ion vacancys
        Dᵢ₀ = 6.5e-8u"m^2/s",      # Diffusion constant
        Eᵢₐ = 0.58 * u"eV",        # Ativation energy of vacancy diffusion
        freeze_ions = false,
        
        # Environment Parameters
        T  = 300u"K",              # Temperature
        α  = 1.3e7u"1/m",          # Perovskite absorption koefficient
        Fₚₕ = 1.4e21u"m^-2*s^-1",  # 1 Sun absorbed photonflux 
        dir = 1,                   # Light trough  1 -> ETL, -1 -> HTL
        light = pulse(tₑ=1.0,w=2.),     # Light(t) function 
        V = t -> 0,                # Voltage(t) function
        Rₛₕ = 1e6u"V/A*m^2",       # Shunt resistance
        Rₛ = 1e-4u"V/A*m^2",        # Series resistance
        # Recombination Parameters
        τₙ = 3e-7u"s",             # electron pseudo lifetime
        τₚ = 3e-7u"s",             # hole pseudo lifetime
        k₂ = 3.22e-17u"m^3/s",     # second order rate constant
    
        # Interface Recombination
        k₂ₑ = 0u"m^4/s",           # ETL/perovskite bimolecular recombination rate
        k₂ₕ = 0u"m^4/s",           # perovskite/HTL bimolecular recombination rate
        vₙₑ = 0u"m/s",             # electron recombination velocity for SHR/ETL
        vₚₑ = 0u"m/s",             # hole recombination velocity for SHR/ETL
        vₙₕ = 0u"m/s",             # electron recombination velocity for SHR/HTL
        vₚₕ = 0u"m/s",             # hole recombination velocity for SHR/HTL
    
        # ELT Parameters
        dₑ = 1e18u"cm^-3",         # ETL effective doping density
        mcₑ = 1.5,                 # ETL effective electron mass
        Ecₑ = -4.0 * u"eV",        # ETL conduction band energy
        bₑ = 100e-9u"m",           # ETL width
        εₑᵣ = 3,                   # ETL permitivity
        Dₑ = 1e-7u"m^2/s",         # ETL electron diffusion coeficcient
    
        # HTL Parameters
        dₕ = 1e18u"cm^-3",         # HTL effective doping density
        mvₕ = 12,                  # HTL hole mass
        Evₕ = -5 * u"eV",          # HTL valence band energy
        bₕ = 100e-9u"m",           # HTL width
        εₕᵣ = 3,                   # HTL permitivity
        Dₕ = 1e-7u"m^2/s",         # HTL electron diffusion coeficcient
    )

Constructs the parameters object from the defaults.

# Examples
```@example; output = false
def_parm = Parameters()         # Use the default parameters
mod_parm = Parameters(          # Use modified default parameters 
    b  = 432u"nm",              # Perovskite Layer thickness
    ε = 42,                     # Perovskite permitivity
)
```
"""
Base.@kwdef mutable struct Parameters <: AbstractParameters
    N::Integer  =   400       # Subintervals in perovskite layer, resulting in N+1 Grid points
    
    # Physical parameters
    ε₀ = 8.854187817e-12u"F/m"# Permitivity of free space
    q  = 1.6021766209e-19u"C" # Elemntary charge of a proton
    kB = 8.61733035e-5u"eV/K" # Bolzmann konstant
   
    
    # Perovskite parameters
    b  = 400e-9u"m"           # Perovskite Layer thickness
    ε = 24.1                  # Perovskite permitivity
    Ec = -3.7u"eV"            # Perovskite Conduction band energy
    Ev = -5.3u"eV"            # Perovskite Valence band energy
    dₚ   = 0u"m^-3"           # Peroviskite doping concentration
    Dₙ = 1.7e-4u"m^2/s"       # Perovskite electron diffusion coefficient
    Dₚ = 1.7e-4u"m^2/s"       # Perovskite hole diffusion coefficient
    mₑ = 0.2                  # Perovskite electron mass
    mₕ = 0.2                  # Perovskite hole mass
    
    # Ion Parameters
    N₀ = 1.6e24u"m^-3"        # Typical density of ion vacancys
    Dᵢ₀ = 6.5e-8u"m^2/s"      # Diffusion constant
    Eᵢₐ = 0.58 * u"eV"        # Ativation energy of vacancy diffusion
    freeze_ions = false
    
    # Environment Parameters
    T  = 300u"K"              # Temperature
    α  = 1.3e7u"1/m"          # Perovskite absorption coefficient
    Fₚₕ = 1.4e21u"m^-2*s^-1"  # 1 Sun photonflux  
    dir::Integer =   1        # Light trough  1 -> ETL, -1 -> HTL
    light::Function = pulse(tₑ=1.0,w=2.)     # Light(t) function  [Sun]
    V::Function = t -> 0               # Voltage(t) function [V]
    Rₛₕ = 1e6u"V/A*m^2"       # Shunt resistance
    Rₛ = 1e-4u"V/A*m^2"        # Series resistance

    # Recombination Parameters
    τₙ = 3e-7u"s"             # electron pseudo lifetime
    τₚ = 3e-7u"s"             # hole pseudo lifetime
    k₂ = 3.22e-17u"m^3/s"     # second order rate constant
   
    # Interface Recombination
    k₂ₑ = 0u"m^4/s"           # ETL/perovskite bimolecular recombination rate
    k₂ₕ = 0u"m^4/s"           # perovskite/HTL bimolecular recombination rate
    vₙₑ = 0u"m/s"             # electron recombination velocity for SHR/ETL
    vₚₑ = 0u"m/s"             # hole recombination velocity for SHR/ETL
    vₙₕ = 0u"m/s"             # electron recombination velocity for SHR/HTL
    vₚₕ = 0u"m/s"             # hole recombination velocity for SHR/HTL
   
    # ELT Parameters
    dₑ = 1e18u"cm^-3"         # ETL effective doping density
    mcₑ = 1.5                 # ETL electron mass
    Ecₑ = -4.0 * u"eV"        # ETL conduction band energy
    bₑ = 100e-9u"m"           # ETL width
    εₑᵣ::Real = 3             # ETL permitivity
    Dₑ = 1e-7u"m^2/s"         # ETL electron diffusion coefficient
   
    # HTL Parameters
    dₕ = 1e18u"cm^-3"         # HTL effective doping density
    mvₕ = 12                  # HTL hole mass
    Evₕ = -5 * u"eV"          # HTL valence band energy
    bₕ = 100e-9u"m"           # HTL width
    εₕᵣ::Real = 3             # HTL permitivity
    Dₕ = 1e-7u"m^2/s"         # HTL electron diffusion coefficient
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
    setproperty!(p, :gc, g);
    setproperty!(p, :gv, g);
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
getproperty(p::AbstractParameters,::Val{:Nc}) = p.gc
getproperty(p::AbstractParameters,::Val{:Nv}) = p.gv
getproperty(p::AbstractParameters,::Val{:Ncₑ}) = p.gcₑ
getproperty(p::AbstractParameters,::Val{:Nvₕ}) = p.gvₕ

Base.@kwdef mutable struct AlgControl
# TODO: make this just a Dict
    # Logging::LogLevel   = Warn    # Sets the log level:  Debug ⊂ Info ⊂ Warn ⊂ Error

    ss_tol::Float64     = 1e-6      # ftol for steady state detection
    reltol::Float64     = 1e-6      # Relative tolerance for the solvers (auto_abstol =  u*reltol)
    abstol              = 1e-6
    autodiff::Bool      = true      # Enable autodifferentiation (highly suggested)
    alg     = Rodas4P2() # Integrator Algorithm. Optional suggestions are: Rodas5(), ROS34PW2(), ROS3P(), QBDF(), QNDF()
    dtmin   = 1e-15
    dt      = 1e-14
    force_dtmin = true
    progress = false
    progress_steps = 50
    maxiters = 5000
    tstart = 0.0u"s"
    tend = 1e5u"s"
    initializealg = OrdinaryDiffEq.NoInit()
    numtype = Float64
    linsolve=OrdinaryDiffEq.DefaultLinSolve()
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
