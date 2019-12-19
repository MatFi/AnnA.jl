#using ForwardDiff #moved to runtests.jl
#using Setfield   #moved to runtests.jl

##Test if the jacobian result is reasonable
@testset "Jacobian" begin
    prm = AnnABase.Parameters()
    cell=AnnABase.Cell(prm; mode=:cc)

    N=AnnABase.length(cell.g)
    x=rand(N)
    dx=zeros(N)
    f!(du,u) = cell.rhs(du,u,cell,0)
    J =ForwardDiff.jacobian(f!,dx,x)
    Jsp_rnd= [!iszero(J[i,j]) for i in 1:N, j in 1:N ]
    J  = AnnABase.get_jac_sparse_pattern(cell.g; mode = cell.mode)
    Jsp= [!iszero(J[i,j]) for i in 1:N, j in 1:N ]
    @test Jsp == Jsp_rnd

    cell=AnnABase.Cell(prm; mode=:oc)
    N=AnnABase.length(cell.g)
    x=rand(N)
    dx=zeros(N)
    f!(du,u) = cell.rhs(du,u,cell,0)
    J =ForwardDiff.jacobian(f!,dx,x)
    Jsp_rnd= [!iszero(J[i,j]) for i in 1:N, j in 1:N ]
    J  = AnnABase.get_jac_sparse_pattern(cell.g; mode = cell.mode)
    Jsp= [!iszero(J[i,j]) for i in 1:N, j in 1:N ]
    @test Jsp == Jsp_rnd

    cell=AnnABase.Cell(prm; mode=:precondition)
    N=AnnABase.length(cell.g)
    x=rand(N)
    dx=zeros(N)
    f!(du,u) = cell.rhs(du,u,cell,0)
    J =ForwardDiff.jacobian(f!,dx,x)
    Jsp_rnd= [!iszero(J[i,j]) for i in 1:N, j in 1:N ]
    J  = AnnABase.get_jac_sparse_pattern(cell.g; mode = cell.mode)
    Jsp= [!iszero(J[i,j]) for i in 1:N, j in 1:N ]
    @test Jsp == Jsp_rnd

    @test_throws ErrorException AnnABase.get_jac_sparse_pattern(cell.g; mode = :err)
end
