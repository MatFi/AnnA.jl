(
    N  =   400,       # Subintervals in perovskite layer, resulting in N+1 Grid points
    
    # Physical parameters
    ε₀ = 8.854187817e-12u"F/m",# Permitivity of free space
    q  = 1.6021766209e-19u"C", # Elemntary charge of a proton
    kB = 8.61733035e-5u"eV/K", # Bolzmann konstant

    # Perovskite parameters
    b  = 400e-9u"m",           # Perovskite Layer thickness
    ε = 24.1,                  # Perovskite permitivity
    Ec = -3.7u"eV",            # Perovskite Conduction band energy
    Ev = -5.3u"eV",            # Perovskite Valence band energy
    dₚ = 0u"m^-3",             # Peroviskite doping concentration
    Dₙ = 1.7e-4u"m^2/s",       # Perovskite electron diffusion coefficient
    Dₚ = 1.7e-4u"m^2/s",       # Perovskite hole diffusion coefficient
    mₑ = 0.2,                  # Perovskite electron mass
    mₕ = 0.2,                  # Perovskite hole mass
    
    # Ion Parameters
    N₀ = 1.6e24u"m^-3",        # Typical density of ion vacancys
    Dᵢ₀ = 6.5e-8u"m^2/s",      # Diffusion constant
    Eᵢₐ = 0.58 * u"eV",        # Ativation energy of vacancy diffusion
    freeze_ions = false,
    
    # Environment Parameters
    T  = 300u"K",              # Temperature
    α  = 1.3e7u"1/m",          # Perovskite absorption coefficient
    Fₚₕ = 1.4e21u"m^-2*s^-1",  # 1 Sun photonflux  
    dir = 1,                   # Light trough  1 -> ETL, -1 -> HTL
    light = pulse(tₑ=1.0, w=2.),     # Light(t) function  [Sun]
    V = t -> 0,                # Voltage(t) function [V]
    Rₛₕ = 1e6u"V/A*m^2",       # Shunt resistance

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
    mcₑ = 1.5,                 # ETL electron mass
    Ecₑ = -4.0 * u"eV",        # ETL conduction band energy
    bₑ = 100e-9u"m",           # ETL width
    εₑᵣ = 3,                   # ETL permitivity
    Dₑ = 1e-7u"m^2/s",         # ETL electron diffusion coefficient
   
    # HTL Parameters
    dₕ = 1e18u"cm^-3",         # HTL effective doping density
    mvₕ = 12,                  # HTL hole mass
    Evₕ = -5 * u"eV",          # HTL valence band energy
    bₕ = 100e-9u"m",           # HTL width
    εₕᵣ = 3,                   # HTL permitivity
    Dₕ = 1e-7u"m^2/s",         # HTL electron diffusion coefficient
)