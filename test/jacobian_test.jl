#using ForwardDiff #moved to runtests.jl
#using Setfield   #moved to runtests.jl

##Test if the jacobian result is reasonable
@testset "Jacobian" begin
    prm = AnnA.Parameters()
    cell=AnnA.Cell(prm; mode=:cc)

    N=AnnA.length(cell.g)
    x=rand(N)
    dx=zeros(N)
    f!(du,u) = cell.rhs(du,u,cell,0)
    J =ForwardDiff.jacobian(f!,dx,x)
    Jsp_rnd= [!iszero(J[i,j]) for i in 1:N, j in 1:N ]
    J  = AnnA.get_jac_sparse_pattern(cell.g; mode = cell.mode)
    Jsp= [!iszero(J[i,j]) for i in 1:N, j in 1:N ]
    @test Jsp == Jsp_rnd

    cell=AnnA.Cell(prm; mode=:oc)
    N=AnnA.length(cell.g)
    x=rand(N)
    dx=zeros(N)
    f!(du,u) = cell.rhs(du,u,cell,0)
    J =ForwardDiff.jacobian(f!,dx,x)
    Jsp_rnd= [!iszero(J[i,j]) for i in 1:N, j in 1:N ]
    J  = AnnA.get_jac_sparse_pattern(cell.g; mode = cell.mode)
    Jsp= [!iszero(J[i,j]) for i in 1:N, j in 1:N ]
    @test Jsp == Jsp_rnd

    cell=AnnA.Cell(prm; mode=:precondition)
    N=AnnA.length(cell.g)
    x=rand(N)
    dx=zeros(N)
    f!(du,u) = cell.rhs(du,u,cell,0)
    J =ForwardDiff.jacobian(f!,dx,x)
    Jsp_rnd= [!iszero(J[i,j]) for i in 1:N, j in 1:N ]
    J  = AnnA.get_jac_sparse_pattern(cell.g; mode = cell.mode)
    Jsp= [!iszero(J[i,j]) for i in 1:N, j in 1:N ]
    @test Jsp == Jsp_rnd

    @test_throws ErrorException AnnA.get_jac_sparse_pattern(cell.g; mode = :err)
end
