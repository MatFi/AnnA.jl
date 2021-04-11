```@meta
CurrentModule = AnnA
```

# Simulation Parameters
The interface is designed for a convenient access to the core functionality. One of the most important aspects the user schould be able to control are the input Parameters. This can be archieved via the [`Parameters()`](@ref) constructor.

## Transient Parameters

the fields `light` and `V` , correspond to the illumination and voltage transients during the simulation. They must be both of type `Function` and take one argument of type `Real` which corresponds to the time in seconds. 

```@example interface; output = false
using AnnA # hide
Parameters(V = t -> 0.5*t) # Voltage ramp with a slope of 0.5V/s

function l(t)
    1/2*(sin(t)+1)
end

Parameters(light = l)      # Sinusiodal light excitation with 1 Sun amplitude
``` 

The solver generaly likes continous functions, but may also work for a broad range of discontinous `light`. 

### Ionic Motion

## Implicit Parameters

### Density Of States

### Built In Potential

### Fermi Level

### Ionic Timescale


