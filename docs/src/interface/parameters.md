```@meta
CurrentModule = AnnA
```

# Simulation Parameters
The interface is designed for a convenient access to the core functionality. One of the most important aspects the user schould be able to control are the input Parameters. This can be archieved via the `Parameters()` constructor.

```@docs
Parameters
```

## Transient Parameters

The fields `light` and `V` , correspond to the illumination and voltage transients during the simulation. They must be both of type `Function` and take one argument of type `Real` which corresponds to the time in seconds. 

```@example interface; output = false
using AnnA # hide
using Unitful # hide 
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

!!! note "Disable ionic motion"
	
	To disable ion movement, use `freeze_ions = true`. The setting `N₀=0`, `Dᵢ₀ = 0` or `Eᵢₐ=Inf` results in zero divisions during non-dimensioning and lead to an error.

## Implicit Parameters

Values for other important device parameters are implicitly provided from the input in [`Parameters()`](@ref). 
### Effective Density Of States
This conduction / valence band DOS ``~	N_{c/v} `` is calculated from the *free electron gas* approximation, electron/hole effective mass ``m_{e/h}``, the Bolzmann constant ``k_B``, temperature ``T``, Electron rest mass ``m_0`` and Planck constant ``h`` via:
```math
N_{c/v} = 2\left(\frac{2\pi m_{e/h}  m_0  k_B T}{h^2}   \right)^{3/2}
```	

```@repl interface
parm = Parameters(mₑ=0.2, mₕ=0.2, T=300u"K");	
parm.Nc
parm.Nv
``` 

### Intrinsic Carrier Density
``n_i`` is defined by the Temperature, ``N_{c/v}`` and the Bandgap ``E_g`` by the mass action law:
```math
n_i^2 = N_{c}N_{v} \cdot e^\frac{- E_g }{k_B T}
```

```@repl interface	
parm.nᵢ
```
### Fermi Level
The Fermilevel ``Ef_{e/h}`` in the electron / hole trasportlayer can be expressed by thair doping density  ``d_{e/h}``, effective DOS and bandenergy ``E_{ce} / E_{vh} ``,  

```math
\begin{aligned}
Ef_{e}  &= E_{ce} - k_BT \log(N_{ce}/d_e) \\
Ef_h    &= E_{vh} + k_BT \log(N_{vh}/d_h)
\end{aligned}
```

```@repl interface	
parm.Efₑ
parm.Efₕ
```

### Built In Potential
The difference in fermi levels between ETM and HTM define the built in potential ``V_{Bi}``

```@repl interface	
parm.Efₑ - parm.Efₕ
parm.Vbi
```


