# Usage
## Solution handling
Handling the solution calculated by the solver is a multidimensional problem. Besides time and space, the problem contains also different variables have also different variables `ϕ(t, x)`, `I(t, x)`, `n(t, x)` and `p(t, x)`. Since we use monotonic interpolation internally, every `t` and `x` can be selected. But mostly the solution as a function of `x` at a distinct time `t`, or vice versa is wanted.

since not all parameters exist for all `x` but for all times `t`. the time dependence is put to the solution object `sol` and the spacial dependence to the variables:

```julia
# get parameter n at position xₚ as function of time
sol.n(xₚ) -> Array of Tuples (tᵢ, n(tᵢ,xₚ))

# get parameter n at time tₚ as function of x
sol(tₚ).n  -> Array of Tuples (xᵢ,n(tₚ,xᵢ)

# get parameter n at position xₚ and time tₚ
sol(tₚ).n(xₚ)  ->  n(tₚ,xₚ)



sol(t).n.x[1:end]
sol(t).n.n[1:end]

```
