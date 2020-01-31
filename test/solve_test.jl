#using Test, AnnABase
#using Logging
@testset "Solve" begin
    prm = AnnABase.Parameters(V = t -> 0)
    #=
    cell = AnnABase.Cell(
        prm;
        mode = :cc,
        alg_ctl = AnnABase.AlgControl(ss_tol = 1e-8),
    )
    AnnABase.initial_conditions!(cell)
    @test AnnABase.integrate(cell.g.x, cell.g(cell.u0)[1][1]) ≈ 1
    =#
    cell = AnnABase.Cell(
        prm;
        mode = :oc,
        alg_ctl = AnnABase.AlgControl(ss_tol = 1e-8, reltol = 1e-6),
    )
    sol = AnnABase.solve(cell)
    @test isapprox(AnnABase.integrate(cell.g.x, cell.g(cell.u0)[1][1]) , 1,atol = 1e-4)


    @test isapprox(
        ustrip(AnnABase.rdim_sol(sol)[2][3][end]),
        ustrip(-cell.ndim.Vbi * cell.parameters.VT),atol = 1e-4
    )

end
