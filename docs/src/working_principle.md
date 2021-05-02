# Working Principle
The implementation is based on Courtiers finite element scheme. The original Matlab implementation [IonMonger](https://github.com/PerovskiteSCModelling/IonMonger) and internally uses the same numerical finite element scheme.

# Physical Concepts
## Fundamental Conservation Laws
In a homogeneous thin-film solar cell where the device area is magnitudes larger than the thickness, it is sufficient to consider only the out of plane dimension. The electrons ``n``, holes ``p``, and ionic species ``P`` obey the conservation laws listed below: 
```math
    \begin{aligned}
        \frac{\partial n}{\partial t} &= \frac{1}{q}\frac{\partial j_n}{\partial x}  +G - R, 
        &j_n = qD^-\left(\frac{\partial n}{\partial x} - \frac{n~q}{k_BT}\frac{\partial \phi}{\partial x}\right) \label{eq:cnt_n}\\
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
## Algebraic Partial Differential Equations
## Recombination mechanisms
## Density of states and doping
## Intrinsic carrier density and temperature

# Couriers Finite Element Scheme
The interested reader should consult Courtiers original publications [^ionmonger]
here we will provide only a few important fundamentals necessary to understand a few important limitations from fundamental physical perspective as well as in some numerical aspects.
## Spatial *tanh* grid
In order to coupe with the spatial stiffness the PDE is descretized by a non uniform grid which is more dense within the ionic accumulation layers at the interface. This helps to overcome numerical inaccuracies the finite difference gradients during the integration of the PDE.  
## Non dimensionalization
# Implementation
## OrdinaryDiffEq.jl Solvers
## Sparsity
## Automatic Differentiation


[^ionmonger]:

    ```
    Courtier, N.E., Cave, J.M., Walker, A.B. et al. IonMonger: a free and fast planar perovskite solar cell simulator with coupled ion vacancy and charge carrier dynamics. J Comput Electron 18, 1435–1449 (2019). https://doi.org/10.1007/s10825-019-01396-2
    ```




