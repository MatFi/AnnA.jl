#using Test, AnnA
#using Logging
@testset "Solve" begin

prm = AnnA.Parameters(
    V = t -> 0,
    light= t->0,Ecₑ=-3.9u"eV",
    N=400, 
    Fₚₕ = 1.4e21u"m^-2*s^-1" ,
    Rₛₕ = 1e6*1u"V/A*m^2",
    bₕ = 100e-9u"m",           # HTL width
    bₑ = 100e-9u"m",           # HTL width 
    freeze_ions=false,
    )

    tend=1e6u"s"
    cell = AnnA.Cell(
        prm;
        mode = :oc,
        alg_ctl = AnnA.AlgControl(ss_tol = 1e-5,tend=tend, reltol = 1e-8,abstol=1e-8),
    )

    sol = AnnA.solve(cell)
    @test isapprox(AnnA.integrate(cell.g.x, cell.g(cell.u0)[1][1]) , 1,atol = 1e-14)

    println("$(tend) long term :oc voltagedrift: ", AnnA.rdim_sol(sol)[end][2][3][end]+cell.ndim.Vbi * cell.parameters.VT)
    @test isapprox(
        ustrip(AnnA.rdim_sol(sol)[end][2][3][end]),
        ustrip(-cell.ndim.Vbi * cell.parameters.VT),atol = 3e-3
    )

end
