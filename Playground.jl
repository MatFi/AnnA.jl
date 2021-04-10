using AnnA
using Plots
using Unitful
parm = Parameters(light = t -> 1.0,
    vₙₕ= 10u"m/s" ,        
    vₚₕ= 0.01u"m/s" ,
    N=500,
    N₀=1e18u"cm^-3"
)
prob = IVProblem(parm, [-0.5,1.7]u"V", 0.2u"V/s")
sol  = solve(prob)

sol=sol.df
sol.j =sol.j .|> u"mA/cm^2"
sol= sol  .|>ustrip

plot(sol.V[sol.fwd], sol.j[sol.fwd] .|>abs, yscale=:log10, ylabel="j (A/m²)",xlabel="Voltage (V)",legend=false)
plot!(sol.V[.!sol.fwd], sol.j[.!sol.fwd] .|>abs )
plot(sol.V[sol.fwd], sol.j[sol.fwd] , ylabel="j (mA/cm²)",xlabel="Voltage (V)",label="Forward", ylims=(-25,40),xlims=(-0.5,1.3),legend=:topleft)
plot!(sol.V[.!sol.fwd], sol.j[.!sol.fwd],label="Backward" )


    using  Unitful #,AnnAPlot, Plots
include("./src/AnnABase.jl")
parameters = AnnABase.Parameters(
    N=300,mₑ =0.2 ,mₕ =0.2, light = t->0,Rₛₕ =Inf*1u"V/A*m^2",
    dₚ=5e19u"m^-3", 
    dₕ=1e23u"m^-3",
    dₑ=1e23u"m^-3",
    vₙₑ= 0*1e-4u"m/s" ,        # electron recombination velocity for SHR/ETL
    vₚₑ= 0*1000u"m/s"  ,        # hole recombination velocity for SHR/ETL
    vₙₕ= 1000u"m/s" ,        # electron recombination velocity for SHR/HTL
    vₚₕ = 1e-4u"m/s" , 
   # Fₚₕ= 1.4e25u"m^-2*s^-1",
    τₙ= 3e-7u"s" ,          #electron pseudo lifetime
    τₚ = 3e-7u"s" , 
#    bₑ = 200u"nm",
 #   bₕ = 500u"nm",
    N₀=1e22u"m^-3   " ,
    Evₕ= -5. *u"eV",  
    Ecₑ= -4. *u"eV",  
    mcₑ = 1.5, 
    mvₕ = 1.5 ,

    freeze_ions=false   )

prob_ocvd = AnnABase.OCVDProblem(
    parameters,
    1u"s",
    1e3u"s",
    alg_control = AnnABase.AlgControl(ss_tol = 1e-6, reltol = 1e-7, abstol = 1e-7,maxiters=10000,dtmin=1e-40, progress=true)
)
sol = AnnABase.solve(prob_ocvd)

prob_iv = AnnABase.IVProblem(
    parameters,
    [-0.3,1.2]u"V",
    0.00000002u"V/s",
    double_sweep=true,
    alg_control  = AnnABase.AlgControl(ss_tol = 1e-6, reltol = 1e-7, abstol = 1e-7,maxiters=10000,dtmin=1e-40, progress=true)
)
sol = AnnABase.solve(prob_iv)

@pgf TikzPicture(
       Axis({ymode="log"},
            PlotInc({ only_marks },
                Table(; x = sol.df.V, y = abs.(sol.df.j))
            ),

        )
    )

sol_decay=AnnABase.ProblemSolution(sol.sol_decay)
using PGFPlotsX
PGFPlotsX.print_tex(io::IO, data::Unitful.AbstractQuantity) = PGFPlotsX.print_tex(io, ustrip(upreferred(data)))
@pgf TikzPicture(
       Axis({ymode="log"},
            PlotInc({ only_marks },
                Table(; x = sol_decay.x, y = sol_decay.df.n[1] )
            ),
            PlotInc(
                Table(; x = sol_decay.x, y = sol_decay.df.p[1])
            )
        )
    )

@pgf TikzPicture(
       Axis({xmode="log"},
            PlotInc({ only_marks },
                Table(; x = sol_decay.df.t, y = sol_decay.df.V)
            ),

        )
    )

sol_decay=AnnABase.ProblemSolution(sol.sol_decay)
@pgf TikzPicture(
        Axis(
             PlotInc({ only_marks },
                 Table(; x = sol_decay.xₑ, y = sol_decay.df.ϕₑ[end] )
             ),
             PlotInc(
                 Table(; x = sol_decay.x, y = sol_decay.df.ϕ[end])
             ),
             PlotInc(
                Table(; x = sol_decay.xₕ, y = sol_decay.df.ϕₕ[end])
            )
         )
     )

sol_decay.df.V[1] 
sol_decay.df.p[end]

ENV["JULIA_DEBUG"]="all"