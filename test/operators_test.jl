@testset "Operators" begin
    cell_a  = AnnA.Cell(AnnA.Parameters();op_flavor=:matrix_free)
    cell_b  = AnnA.Cell(AnnA.Parameters();op_flavor=:dense)
    cell_c  = AnnA.Cell(AnnA.Parameters();op_flavor=:sparse)


    ϕ   = rand(cell_a.g.N+1)
    mE_free = deepcopy(DiffEqBase.get_tmp(cell_a.rhs.mE,ϕ))

    mE_dense = deepcopy(mE_free)
    mE_sparse = deepcopy(mE_free)
    mul!(mE_free,cell_a.o.𝔇,ϕ)
    mul!(mE_dense,cell_b.o.𝔇,ϕ)
    mul!(mE_sparse,cell_c.o.𝔇,ϕ)
    @test mE_dense ≈ mE_free
    @test mE_dense ≈ mE_sparse

    mul!(mE_free,cell_a.o.𝕴,ϕ)
    mul!(mE_dense,cell_b.o.𝕴,ϕ)
    mul!(mE_sparse,cell_c.o.𝕴,ϕ)
    @test mE_dense ≈ mE_free
    @test mE_dense ≈ mE_sparse

    #ϕ   = rand(length(cell_a.g.x))
    cd_free = DiffEqBase.get_tmp(cell_a.rhs.cd,ϕ)
    cd_dense =DiffEqBase.get_tmp(cell_b.rhs.cd,ϕ)
    cd_sparse =DiffEqBase.get_tmp(cell_c.rhs.cd,ϕ)
    mul!(cd_free,cell_a.o.𝔏,ϕ)
    mul!(cd_dense,cell_b.o.𝔏,ϕ)
    mul!(cd_sparse,cell_c.o.𝔏,ϕ)
    @test cd_dense ≈ cd_free
    @test cd_dense ≈ cd_sparse

    Test.@test_throws ErrorException AnnA.Cell(
                                            AnnA.Parameters();op_flavor=:err)
    Test.@test_throws ErrorException mul!(mE_free,cell_a.o.𝕴, ϕ[2:end] )
    Test.@test_throws ErrorException mul!(mE_free,cell_a.o.𝔇, ϕ[2:end] )
    Test.@test_throws ErrorException mul!(cd_free,cell_a.o.𝔏, ϕ[2:end] )

end
