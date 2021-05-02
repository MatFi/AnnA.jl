# Working Principle
The implementation is based on Courtiers finite element scheme. The original Matlab implementation [IonMonger](https://github.com/PerovskiteSCModelling/IonMonger) and internally uses the same numerical finite element scheme.

# Physical Concepts
## Fundamental Conservation Laws
In a homogeneous thin-film solar cell where the device area is magnitudes larger than the thickness, it is sufficient to consider only the out of plane dimension. The electrons ``n``, holes ``p``, and ionic species ``P`` obey the conservation laws listed below: 
```math
    \begin{aligned}
        \frac{\partial n}{\partial t} &= \frac{1}{q}\frac{\partial j_n}{\partial x}  +G - R, 
        &j_n = qD^-\left(\frac{\partial n}{\partial x} - \frac{n~q}{k_BT}\frac{\partial \phi}{\partial x}\right) \\
        \frac{\partial p}{\partial t} &= \frac{-1}{q}\frac{\partial j_p}{\partial x} +G - R, 
        &j_p = -qD^+\left(\frac{\partial p}{\partial x} + \frac{p~q}{k_BT}\frac{\partial \phi}{\partial x}\right) \\
        \frac{\partial P}{\partial t} &= \frac{-1}{q}\frac{\partial j_{ion}}{\partial x}, 
        &j_{ion} = -qD_{ion}\left(\frac{\partial P}{\partial x} + \frac{P~q}{k_BT}\frac{\partial \phi}{\partial x}\right)   
    \end{aligned}
``` 
The respective current densities ``j_n``, ``j_p`` and ``j_{ion}`` are composed of diffusion and field driven components, whereby the electric field is provided by the Poisson equation: 
```math
    0 =\frac{\partial^2\phi}{\partial x^2} -\frac{q(N_{ion}-P+n-p-d_n+d_p)}{\varepsilon_{P/E/H}}
```
### Transport layers
For the transport layers simplified equations were applied. These specifically consisted of the assumptions that the HTM / ETM only conducts holes / electrons ``(D^- / D^+ = 0 )`` and zero charge generation ``G`` and recombination ``R`` takes place, as well as there are no mobile ions present ``(P=N_{ion}=0)``.
In case of the active layer we assumed an absence of doping, so that ``d_n = d_p = 0``. ``D^-`` and ``D^+`` are the corresponding diffusion coefficients.
For neutrality reasons the Poisson equation must also contain the immobile counterpart of the mobile ionic species, which is assumed to be uniformly distributed in the active layer.
### Generation profile
For charge carrier generation ``G(x)`` a simple Lambert Beer absorption of the incident photon flux ``F-{ph}`` is assumed
```math
    G(x) = F_{ph} e^{-\alpha x}.
```
Its important to note, that this assumption may be incorrect in case of thin devices where interference can happen.
### Recombination mechanisms
The bulk recombination behavior in the simulations is composed of Shokley-Read-Hall, and bimolecular recombination 
```math
    R\bigl(n,p\bigr) = \frac{(n\cdot p - n_i)}{\tau_p \cdot n+\tau_n \cdot p}+k_{rad}\cdot(n\cdot p - n_i)
```
For surface recombination a similar expression is assumed. Here, the charge carriers from the transport layers are assumed to recombine together with their opposite ones from the active layer. We also allow the surface recombinationion velocities ``\nu_{p_E}, \nu_{n_E}`` at the ETM/active layer interface to differ from ``\nu_{p_H}, \nu_{n_H}`` at the active layer/HTM interface.
```math
    R_{E/H}\bigl(n_{l/r},p_{r/l}\bigr) = \frac{(n_{l/r}\cdot p_{r/l} - n_i)}{\nu_{p_{E/H}} \cdot n_{l/r}+\nu_{n_{E/H}} \cdot p_{r/l}}
```
where ``n_{l/r}, p_{l/r}``, indicates the carrier densities left/right of the interface if the layer sequence ETM-Intrinsic-HTM is present.

        
### Density of states and doping
### Intrinsic carrier density and temperature

# Couriers Finite Element Scheme
The interested reader should consult Courtiers original publications [^1]
here we will provide only a few important fundamentals necessary to understand a few important limitations from fundamental physical perspective as well as in some numerical aspects.
## Spatial *tanh* grid
In order to coupe with the spatial stiffness the PDE is descretized by a non uniform grid which is more dense within the ionic accumulation layers at the interface. This helps to overcome numerical inaccuracies the finite difference gradients during the integration of the PDE.  
## Non dimensionalization
# Implementation
## OrdinaryDiffEq.jl Solvers
## Sparsity
## Automatic Differentiation


[^1]:

    ```
    Courtier, N.E., Cave, J.M., Walker, A.B. et al.
    IonMonger: a free and fast planar perovskite solar cell simulator with coupled ion vacancy and charge carrier dynamics.
    J Comput Electron 18, 1435–1449 (2019). https://doi.org/10.1007/s10825-019-01396-2
    ```




