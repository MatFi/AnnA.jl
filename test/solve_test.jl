#using Test, AnnA
#using Logging
@testset "Solve" begin
    prm = AnnA.Parameters(V = t -> 0)
    #=
    cell = AnnA.Cell(
        prm;
        mode = :cc,
        alg_ctl = AnnA.AlgControl(ss_tol = 1e-8),
    )
    AnnA.initial_conditions!(cell)
    @test AnnA.integrate(cell.g.x, cell.g(cell.u0)[1][1]) ≈ 1
    =#
    cell = AnnA.Cell(
        prm;
        mode = :oc,
        alg_ctl = AnnA.AlgControl(ss_tol = 1e-8, reltol = 1e-6),
    )
    sol = AnnA.solve(cell)
    @test isapprox(AnnA.integrate(cell.g.x, cell.g(cell.u0)[1][1]) , 1,atol = 1e-4)


    @test isapprox(
        ustrip(AnnA.rdim_sol(sol)[end][2][3][end]),
        ustrip(-cell.ndim.Vbi * cell.parameters.VT),atol = 1e-4
    )

end
