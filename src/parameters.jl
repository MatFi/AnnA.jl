abstract type AbstractParameters end

Base.@kwdef struct Parameters <: AbstractParameters

    """Model Parameters"""

    N::Integer  =   400     # Subintervals in perovskite layer, resulting in N+1 Grid points
    rtol::Float64 = 1e-6    # Relative tolerance for the solver (abtol is automatic)

    #Physical parameters
    ε₀::Number = 8.854187817e-12u"F/m"  # Permitivity of free space
    q::Number  = 1.6021766209e-19u"C"   # Elemntary charge of a proton
    kB::Number = q*8.61733035e-5u"V/K"  # Bolzmann konstant
    T::Number  = 300u"K"                # Temperature

    # Perovskite parameters
    b::Number  = 400e-9u"m"           # Perovskite Layer thickness
    εₚ::Number = 24.1*ε₀            # Perovskite permitivity
    Ec::Number = -3.7*q*u"V"        # Perovskite Conduction band energy
    Ev::Number = -5.3*q*u"V"        # Perovskite Valence band energy
    Eg::Number = Ec-Ev              # Perovskite Bandgap
    Dₙ::Number = 1.7e-4u"m^2/s"     # Perovskite electron diffusion coefficient
    Dₚ::Number = 1.7e-4u"m^2/s"     # Perovskite hole diffusion coefficient
    gc::Number = 8.1e24u"m^-3"      # Perovskite conduction band DOS
    gv::Number = 5.8e24u"m^-3"      # Perovskite valence band DOS
    # Ion Parameters
    N₀::Number = 1.6e25u"m^-3"      # Typical density of ion vacancys
    Dᵢ₀::Number= 6.5e-8u"m^2/s"     # Diffusion constant
    Eᵢₐ::Number= 0.58*q*u"V"        # Ativation energy of vacancy diffusion
    Dᵢ::Number = Dᵢ₀*exp(-Eᵢₐ/(kB*T)) #Diffusion coefficient of ions

    # Light Parameters
    α::Number  = 1.3e7u"1/m"        # Perovskite absorption koefficient
    Fₚₕ::Number= 1.4e21u"m^-2*s^-1" # Photonflux (bandgap dependent)
    dir::Integer =   1              # Shine light trough  1 -> ETL, -1 -> HTL
    light::Function = pulse(tₑ=1.0)

    # Recombination Parameters
    nᵢ::Number = sqrt(gc*gv)*exp(-Eg/(2*kB*T))  #Intrinsic carrier density
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
    Ecₑ::Number= -4*q*u"V"          # ETL conduction band energy
    bₑ::Number = 100e-9u"m"         # ETL width
    εₑ::Number = 3*ε₀               # ETL permitivity
    Dₑ::Number = 1e-5u"m^2/s"       # ETL electron diffusion coeficcient
    Efₑ::Number = Ecₑ-kB*T*log(gcₑ /dₑ) #HTL Fermi level on edge
    # HTL Parameters
    dₕ::Number = 1e24u"m^-3"        # HTL effective doping density
    gvₕ::Number= 5e25u"m^-3"        # HTL valence band DOS
    Evₕ::Number= -5*q*u"V"          # HTL valence band energy
    bₕ::Number = 100e-9u"m"         # HTL width
    εₕ::Number = 3*ε₀               # HTL permitivity
    Dₕ::Number = 1e-5u"m^2/s"       # HTL electron diffusion coeficcient
    Efₕ::Number = Evₕ+kB*T*log(gvₕ /dₕ) #HTL Fermi level on edge

    V::Function = t-> 0#uconvert(Unitful.NoUnits,(Efₑ-Efₕ)/1u"V")      # Voltage protocoll
    Rₛₕ::Number = 1e3u"V/A*m^2"     # shunt resistance
    """Internal Parameter"""
    """Dont Touch"""
    Vbi::Number = Efₑ-Efₕ       #Build In Potential
    τᵢ::Number = b/Dᵢ*sqrt(kB*T*εₚ/N₀/q^2) #Characteristic time scale of ion motion
    G₀::Number = Fₚₕ/b*(1-exp(-α*b))        #typical rate of photogeneration
#    G::Function = (x,t) -> @. light(t)*α*b/(1-exp(-α*b))*exp(-α*b*(dir*x + (1-dir)/2));
    VT::Number = kB*T/q      #Thermal voltage
    jay::Number = q*G₀*b  # rescaling factor for current density to be in mAcm^-2

    n0::Number = dₑ* gc/gcₑ*exp((Ecₑ-Ec)/(kB*T)) #typical electron density in perovskite
    p0::Number = dₕ* gv/gvₕ*exp((Ev-Evₕ)/(kB*T)) #typical hole density in perovskite

end


Base.@kwdef struct AlgControl

    #Logging::LogLevel   = Warn    # Sets the log level:  Debug ⊂ Info ⊂ Warn ⊂ Error

    ss_tol::Float64     = 1e-6      # ftol for steady state detection
    reltol::Float64     = 1e-6      # Relative tolerance for the solvers (auto_abstol =  u*reltol)
    autodiff::Bool      = true      # Enable autodifferentiation (higly suggested)
    alg     = Rodas4P(autodiff=autodiff) # Integrator Algorithm. Optinal suggestions are: Rodas5(), ROS34PW2(), ROS3P(), QBDF(), QNDF()
end
