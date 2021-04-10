# AnnA.jl
&nbsp; &nbsp; &nbsp; &nbsp; _a shameless rewrite of Courtier's [IonMonger](https://github.com/PerovskiteSCModelling/IonMonger) in Julia_
___

| **Documentation**          | **Build Status**        |
|:--------------------------:|:-----------------------:|
| [![][docs-dev-img]][docs-dev-url] | [![][ci-img]][i-url] [![][codecov-img]][codecov-url]|

## Purpose and Aim of this package
This package implements basic functionality to conduct time resolved drift-diffusion simulations of ionic solar cells like perovskites and related systems in one spacial dimension. It applies Courtier's finite element scheme and solves the resulting equations using the awesome [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) package by [Chris Rackauckas](https://github.com/JuliaDiffEq/DifferentialEquations.jl/commits?author=ChrisRackauckas).

## Usage in short
 First load the package and define the parameter set for the simulation run.
 The defaults, are stored in /src/Parameters.jl and will be overwritten by giving the corresponding values as arguments to the `Parameters()` constructor
```julia
using AnnA, Unitful

parameters = AnnABase.Parameters(light = t->0,Rₛₕ =Inf*1u"V/A*m^2" )
```
The whole package uses Unitful.jl to handle the units nicely. This also means that all parameters must carry a unit (subtype of `Unitful.AbstractQuantity`). The unit can easily set by the `Unitful.@u_str` macro.
For the simplest experimental protocols an wrapper is implemented to create them conveniently, e.g. a IV-Sweep as `IVProblem(p::AbstractParameters, range::AbstractArray, sweep_rate::AbstractQuantity)`:
```julia
prob = AnnABase.IVProblem(parameters,[-0.2,1.4]u"V",0.005u"V/s")
AnnABase.solve!(prob)
```



[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://matfi.github.io/AnnA.jl/dev/

[codecov-img]: https://codecov.io/gh/MatFi/AnnA.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/MatFi/AnnA.jl/
[ci-img]: https://github.com/MatFi/AnnA.jl/actions/workflows/CI.yml/badge.svg?branch=master
[ci-url]: https://github.com/MatFi/AnnA.jl/actions/workflows/CI.yml

[issues-url]: https://github.com/MatFi/AnnA.jl/issues