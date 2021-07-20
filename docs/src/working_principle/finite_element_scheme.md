# Courtiers Finite Element Scheme
The interested reader should consult Courtiers original publications [^1]
here we will provide only a few important fundamentals necessary to understand a few important limitations from fundamental physical perspective as well as in some numerical aspects.
## Spatial *tanh* grid
In order to coupe with the spatial stiffness the PDE is discretized by a non uniform grid which is more dense within the ionic accumulation layers at the interface. This helps to overcome numerical inaccuracies the finite difference gradients during the integration of the PDE.  
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




