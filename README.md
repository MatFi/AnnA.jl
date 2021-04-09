# AnnA.jl
&nbsp; &nbsp; &nbsp; &nbsp; _a shameless rewrite of Courtier's [IonMonger](https://github.com/PerovskiteSCModelling/IonMonger) in Julia_
___

## Heavy development
we trying to fill the docs: [documentation (placeholder)](https://matfi.gitlab.io/AnnABase.jl)

## Purpose and Aim of this package
This package implements all the basic functionality to conduct time resolved drift-diffusion simulations of perovksite solar cells and related systems in one spacial dimension. It applies Courtier's finite element scheme and solves the resulting equations using the awesome [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) package by [Chris Rackauckas](https://github.com/JuliaDiffEq/DifferentialEquations.jl/commits?author=ChrisRackauckas).

## Usage in short

 First load the package and define the parameter set for the simulation run.
 The Defaults, are stored in /src/Parameters.jl and will be overwritten by giving the corresponding values as arguments to `Parameters()`
```julia
using AnnA, Unitful

parameters = AnnABase.Parameters(light = t->0,Rₛₕ =Inf*1u"V/A*m^2" )
```
The whole Package uses Unitful.jl to handle the units nicely. This also means that all parameters must carry a unit (subtype of `Unitful.AbstractQuantity`). The Unit can easily set using the `Unitful.@u_str` macro.
For the simplest experimental protocols an wrapper is implemented to create them conveniently, e.g. a IV-Sweep as `IVProblem(p::AbstractParameters, range::AbstractArray, sweep_rate::AbstractQuantity)`:
```julia
prob = AnnABase.IVProblem(parameters,[-0.2,1.4]u"V",0.005u"V/s")
AnnABase.solve!(prob)
```
