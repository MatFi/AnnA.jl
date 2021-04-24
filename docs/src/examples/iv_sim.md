# IV Simulations
The simulation of an IV curve can be done by the `IVProblem` constructor


```@docs
IVProblem
```
Here is a simple example on how to us
```@example iv
using AnnA
using Unitful, UnitfulRecipes
parm = Parameters(light = t -> 1.0,   
    vₙₕ= 10u"m/s" ,                 # electron surface recombination vel. at HTM
    vₚₕ= 0.01u"m/s" ,               # hole surface recombination vel. at HTM
    N=500,                          # grid size
    N₀=1e18u"cm^-3"                 # ionic concentration
)
prob = IVProblem(parm, [-0.5,1.7]u"V", 0.2u"V/s")
sol  = solve(prob)
```
The `ProblemSolution` object contains also grid and spatial information. All timesteps are stored in a `DataFrame` an can be acessed via the `df` field of `sol`.
If we just want to plot the IV characteristics we can do: 

```@example iv
using Plots
using UnitfulRecipes  # To interface Unitful with Plots

sol=sol.df
sol.j = sol.j |>u"mA/cm^2" # scale to common units
plt = plot(sol.V[sol.fwd], sol.j[sol.fwd],label="Forward", ylims=(-25,40),xlims=(-0.5,1.3),legend=:topleft);
plot!(plt,sol.V[.!sol.fwd], sol.j[.!sol.fwd],label="Backward");
plt
```

# Jsc vs. Voc Curve
The wrapper implements ``j_{sc}(V_{oc})`` simulation by two consecutive simulation runs, where the illumination is increased exponentially over time. This reflects the aspect of a slow ``V_{oc}`` built up under low illumination intensities.


# Open Circuit Voltage Decay (OCVD)

To simulate a open circuit voltage decay a `OCVDProblem` is implemented:

```@example iv
prob_ocvd = OCVDProblem(
    parm,       # input parameter set
    50u"s",     # illumination time
    1e5u"s",    # time the decay will be simulated to
) 

sol = solve(prob_ocvd)  
plot(sol.t_decay,sol.V_decay,xscale=:log10)
```