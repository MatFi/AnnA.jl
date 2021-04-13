# IV Simulations
The simulation of an IV curve can be done by the `IVProblem` constructor


```@docs
IVProblem
```
Here is a simple example on how to us
```@example iv
using AnnA
using Unitful
parm = Parameters(light = t -> 1.0,   
    vₙₕ= 10u"m/s" ,                 # electron surface recombination vel. at HTM
    vₚₕ= 0.01u"m/s" ,                   # hole surface recombination vel. at HTM
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

sol=sol.df
sol.j =sol.j .|> u"mA/cm^2"     # scale j to conveniant units
sol= sol  .|>ustrip             # strip the units, so we can plot the result with Plots

plt = plot(sol.V[sol.fwd], sol.j[sol.fwd] , ylabel="j (mA/cm²)",xlabel="Voltage (V)",label="Forward", ylims=(-25,40),xlims=(-0.5,1.3),legend=:topleft);
plot!(plt,sol.V[.!sol.fwd], sol.j[.!sol.fwd],label="Backward");
plt
```

# Jsc vs. Voc Curve
The wrapper implements ``j_{sc}(V_{oc})`` simulation by two consecutive simulation runs, where the illumination is increased exponentially over time. This reflects the aspect of a slow ``V_{oc}`` built up under low illumination intensities.

