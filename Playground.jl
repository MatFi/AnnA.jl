using  Unitful #,AnnAPlot, Plots
include("./src/AnnABase.jl")
parameters = AnnABase.Parameters(
    N=300,mₑ =1 ,mₕ =1, light = t->0,Rₛₕ =Inf*1u"V/A*m^2",
    dₚ=1e21u"m^-3",
    vₙₑ= 0u"m/s" ,        # electron recombination velocity for SHR/ETL
    vₚₑ= 0u"m/s"  ,        # hole recombination velocity for SHR/ETL
    vₙₕ= 0u"m/s" ,        # electron recombination velocity for SHR/HTL
    vₚₕ = 0u"m/s" , 
    N₀=1e22u"m^-3   "    )

prob_ocvd = AnnABase.OCVDProblem(
    parameters,
    1u"s",
    1e3u"s",
    alg_control = AnnABase.AlgControl(ss_tol = 1e-6, reltol = 1e-6, abstol = 1e-6,maxiters=10000,dtmin=1e-40, progress=true)
)
sol = AnnABase.solve(prob_ocvd)

prob_iv = AnnABase.IVProblem(
    parameters,
    [-0.3,1.2]u"V",
    0.005u"V/s",
    double_sweep=true,
    alg_control = AnnABase.AlgControl(ss_tol = 1e-6, reltol = 1e-6, abstol = 1e-6,maxiters=10000,dtmin=1e-40, progress=true)
)
sol = AnnABase.solve(prob_iv)
sol_decay=AnnABase.ProblemSolution(sol.sol_decay)
using PGFPlotsX
PGFPlotsX.print_tex(io::IO, data::Unitful.AbstractQuantity) = PGFPlotsX.print_tex(io, ustrip(upreferred(data)))
@pgf TikzPicture(
       Axis({ymode="log"},
            PlotInc({ only_marks },
                Table(; x = sol_decay.x, y = sol_decay.df.n[end] )
            ),
            PlotInc(
                Table(; x = sol_decay.x, y = sol_decay.df.p[end])
            )
        )
    )

sol_decay.df.n[end] 
sol_decay.df.p[end]

ENV["JULIA_DEBUG"]="all"