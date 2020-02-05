abstract type AbstractParameters end

Base.@kwdef mutable struct Parameters <: AbstractParameters

    """Model Parameters"""

    N::Integer  =   400     # Subintervals in perovskite layer, resulting in N+1 Grid points
    rtol::Float64 = 1e-6    # Relative tolerance for the solver (abtol is automatic)

    #Physical parameters
    ε₀::Number = 8.854187817e-12u"F/m"  # Permitivity of free space
    q::Number  = 1.6021766209e-19u"C"   # Elemntary charge of a proton
    kB::Number = 8.61733035e-5u"eV/K"  # Bolzmann konstant
    T::Number  = 300u"K"                # Temperature

    # Perovskite parameters
    b::Number  = 400e-9u"m"           # Perovskite Layer thickness
    ε::Number = 24.1
    #εₚ::Number = ε*ε₀            # Perovskite permitivity
    Ec::Number = -3.7u"eV"        # Perovskite Conduction band energy
    Ev::Number = -5.3u"eV"        # Perovskite Valence band energy
    #Eg::Number = Ec-Ev              # Perovskite Bandgap
    Dₙ::Number = 1.7e-4u"m^2/s"     # Perovskite electron diffusion coefficient
    Dₚ::Number = 1.7e-4u"m^2/s"     # Perovskite hole diffusion coefficient
    gc::Number = 8.1e24u"m^-3"      # Perovskite conduction band DOS
    gv::Number = 5.8e24u"m^-3"      # Perovskite valence band DOS
    # Ion Parameters
    N₀::Number = 1.6e25u"m^-3"      # Typical density of ion vacancys
    Dᵢ₀::Number= 6.5e-8u"m^2/s"     # Diffusion constant
    Eᵢₐ::Number= 0.58*q*u"V"        # Ativation energy of vacancy diffusion
#    Dᵢ::Number = Dᵢ₀*exp(-Eᵢₐ/(kB*T)) #ext #Diffusion coefficient of ions

    # Light Parameters
    α::Number  = 1.3e7u"1/m"        # Perovskite absorption koefficient
    Fₚₕ::Number= 1.4e21u"m^-2*s^-1" # Photonflux (bandgap dependent)
    dir::Integer =   1              # Shine light trough  1 -> ETL, -1 -> HTL
    light::Function = pulse(tₑ=1.0)

    # Recombination Parameters
    #nᵢ::Number = sqrt(gc*gv)*exp(-Eg/(2*kB*T)) #ext  #Intrinsic carrier density
    τₙ::Number = 3e-9u"s"           #electron pseudo lifetime
    τₚ::Number = 3e-7u"s"           # hole pseudo lifetime
    k₂::Number = 0u"m^3/s"          # second order rate constant
    # Interface Recombination
    k₂ₑ::Number = 0u"m^4/s"         # ETL/perovskite bimolecular recombination rate
    k₂ₕ::Number = 0u"m^4/s"         # perovskite/HTL bimolecular recombination rate
    vₙₑ::Number = 1e5u"m/s"         # electron recombination velocity for SHR/ETL
    vₚₑ::Number = 10u"m/s"          # hole recombination velocity for SHR/ETL
    vₙₕ::Number = 0.1u"m/s"         # electron recombination velocity for SHR/HTL
    vₚₕ::Number = 1e5u"m/s"         # hole recombination velocity for SHR/HTL

    # ELT Parameters
    dₑ::Number = 1e24u"m^-3"        # ETL effective doping density
    gcₑ::Number= 5e25u"m^-3"        # ETL conduction band DOS
    Ecₑ::Number= -4.0*u"eV"          # ETL conduction band energy
    bₑ::Number = 100e-9u"m"         # ETL width
    εₑᵣ::Real = 3
    #εₑ::Number = εₑᵣ*ε₀               # ETL permitivity
    Dₑ::Number = 1e-5u"m^2/s"       # ETL electron diffusion coeficcient
    #Efₑ::Number = Ecₑ-kB*T*log(gcₑ /dₑ) #ext #HTL Fermi level on edge
    # HTL Parameters
    dₕ::Number = 1e24u"m^-3"        # HTL effective doping density
    gvₕ::Number= 5e25u"m^-3"        # HTL valence band DOS
    Evₕ::Number= -5*u"eV"          # HTL valence band energy
    bₕ::Number = 100e-9u"m"         # HTL width
    εₕᵣ::Real = 3
#    εₕ::Number = εₕᵣ*ε₀               # HTL permitivity
    Dₕ::Number = 1e-5u"m^2/s"       # HTL electron diffusion coeficcient
    #Efₕ::Number = Evₕ+kB*T*log(gvₕ /dₕ) #ext #HTL Fermi level on edge

    V::Function = t-> 0#uconvert(Unitful.NoUnits,(Efₑ-Efₕ)/1u"V")      # Voltage protocoll
    Rₛₕ::Number = 1e3u"V/A*m^2"     # shunt resistance

end

propertynames(p::AbstractParameters,private=true) = begin
    pubs = fieldnames(typeof(p))
    privs =(:εₚ,:εₕ,:εₑ,:Eg,:Dᵢ,:nᵢ,:Efₑ,:Efₕ,:Vbi,:VT,:τᵢ,:G₀,:jay,:n0,:p0)
    private && return (pubs...,privs...)
    (pubs...,)
end

setproperty(p::AbstractParameters,n::Symbol) = setproperty(p::AbstractParameters,Val{n}())
setproperty(p::AbstractParameters,::Val{S}) where {S} = setfield!(p,S)
setproperty(p::AbstractParameters, ::Val{:nᵢ}) = begin
    @warn "nᵢ is transformed to gc nd gv via gv= gc =2*kB*T*ln(nᵢ)/Eg"
    g = 2*log(nᵢ)/p.Eg*p.kB*p.T
    setfield!(p,:gc,g)
    setfield!(p,:gv,g)
end

getproperty(p::AbstractParameters,n::Symbol) = getproperty(p::AbstractParameters,Val{n}())
getproperty(p::AbstractParameters,::Val{S}) where {S} = getfield(p,S)
getproperty(p::AbstractParameters,::Val{:εₚ}) = p.ε*p.ε₀
getproperty(p::AbstractParameters,::Val{:εₕ}) = p.εₕᵣ*p.ε₀
getproperty(p::AbstractParameters,::Val{:εₑ}) = p.εₑᵣ*p.ε₀
getproperty(p::AbstractParameters,::Val{:Eg}) = p.Ec-p.Ev
getproperty(p::AbstractParameters,::Val{:Dᵢ}) = p.Dᵢ₀*exp(-p.Eᵢₐ/(p.kB*p.T)) |> u"m^2/s"
getproperty(p::AbstractParameters,::Val{:nᵢ}) = sqrt(p.gc*p.gv)*exp(-p.Eg/(2*p.kB*p.T))
getproperty(p::AbstractParameters,::Val{:Efₑ}) = p.Ecₑ-p.kB*p.T*log(p.gcₑ /p.dₑ) |> u"eV"
getproperty(p::AbstractParameters,::Val{:Efₕ}) = p.Evₕ+p.kB*p.T*log(p.gvₕ /p.dₕ) |> u"eV"
getproperty(p::AbstractParameters,::Val{:Vbi}) = (p.Efₑ-p.Efₕ) |> u"eV"
getproperty(p::AbstractParameters,::Val{:VT}) = p.kB*p.T/p.q  |> u"V"
getproperty(p::AbstractParameters,::Val{:τᵢ}) = p.b/p.Dᵢ*sqrt(p.kB*p.T*p.εₚ/p.N₀/p.q^2) |> u"s"
getproperty(p::AbstractParameters,::Val{:G₀}) = p.Fₚₕ/p.b*(1-exp(-p.α*p.b)) |> upreferred
getproperty(p::AbstractParameters,::Val{:jay}) = p.q*p.G₀*p.b |> u"A/m^2"
getproperty(p::AbstractParameters,::Val{:n0}) = p.dₑ* p.gc/p.gcₑ*exp((p.Ecₑ-p.Ec)/(p.kB*p.T))  |> u"m^-3"
getproperty(p::AbstractParameters,::Val{:p0}) = p.dₕ* p.gv/p.gvₕ*exp((p.Ev-p.Evₕ)/(p.kB*p.T)) |> u"m^-3"


Base.@kwdef struct AlgControl

    #Logging::LogLevel   = Warn    # Sets the log level:  Debug ⊂ Info ⊂ Warn ⊂ Error

    ss_tol::Float64     = 1e-6      # ftol for steady state detection
    reltol::Float64     = 1e-6      # Relative tolerance for the solvers (auto_abstol =  u*reltol)
    abstol              = 1e-9
    autodiff::Bool      = true      # Enable autodifferentiation (higly suggested)
    alg     = Rodas4P(autodiff=autodiff) # Integrator Algorithm. Optinal suggestions are: Rodas5(), ROS34PW2(), ROS3P(), QBDF(), QNDF()
    dtmin   = 1e-15
    dt      = 1e-14
    force_dtmin = true
    progress = true
    progress_steps = 50
    maxiters = 5000
    tend = 1e5u"s"

    kwargs = (; :reltol => reltol,
        :abstol => abstol,
        :dtmin => dtmin,
        :dt => dt,
        :force_dtmin => force_dtmin,
        :progress_steps => progress_steps,
        :progress => progress,
        :maxiters => maxiters)

end
