```@meta
CurrentModule = AnnA
```

# Simulation Parameters
The interface is designed for a convenient access to the core functionality. One of the most important aspects the user schould be able to control are the input Parameters. This can be archieved via the [`Parameters()`](@ref) constructor.
## Transient Parameters

The fields `light` and `V` , correspond to the illumination and voltage transients during the simulation. They must be both of type `Function` and take one argument of type `Real` which corresponds to the time in seconds. 

```@example interface; output = false
using AnnA # hide
Parameters(V = t -> 0.42)   # Constant voltage of 420 mV
Parameters(V = t -> 0.5*t)  # Voltage ramp with a slope of 0.5V/s

function l(t)
    1/2*(sin(t)+1)
end

Parameters(light = l)       # Sinusiodal light excitation with 1 Sun amplitude
``` 

The solver generaly likes continous functions, but may also work for a broad range of discontinous `light`. 

### Ionic Motion

The mobile ions are controlled by the concentration `N₀`, diffusion constant `Dᵢ₀` with the activation energy
`Eᵢₐ`. 

!!! note To disable ion movement, use `freeze_ions = true`. The setting `N₀=0`, `Dᵢ₀ = 0` or `Eᵢₐ=Inf` results in zero divisions during non-dimensioning and lead to an error.

## Implicit Parameters

Values for other important device parameters are implicitly provided from the input parameters in  [`Parameters()`](@ref). 
### Effective Density Of States

### Built In Potential

### Fermi Level

### Ionic Timescale


