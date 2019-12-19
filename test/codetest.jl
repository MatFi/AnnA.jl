using AnnABase
using Test
using BenchmarkTools
using LinearAlgebra
using Setfield
N=100

a= AnnABase.Parameters()
#a=setproperties(a,N=10)
b=AnnABase.NonDimensionalize(a)
@test_broken typeof(b) == AnnABase.NodimParameters
#test if structurs sharing elements (problematic in case of unpacking)
inarray = (a,b) ->any([in(i,b) for i in a])
@test !inarray(fieldnames(typeof(a)),fieldnames(typeof(b)))

c=AnnABase.Grid(a)
@test typeof(c) <: AnnABase.Grid
@test !inarray(fieldnames(typeof(b)),fieldnames(typeof(c)))
d=AnnABase.Operators(c;flavor=:dense)
    @test typeof(d) <: AnnABase.AbstractOperators

M=AnnABase.mass_matrix(c,b)

cell= AnnABase.Cell(a)

using RecursiveArrayTools
function bench()

    a= AnnABase.Parameters()
    cell= AnnABase.Cell(a; op_flavor=:sparse)
    u =ones(length(Matrix(cell.M)[:,1]))
    du =similar(u)
    #u= AnnABase.initial_conditions(cell).zero
    du = similar(u)
#    @code_warntype     AnnABase.rhs!(du,u,cell,0)
    @btime AnnABase.rhs!($du,$u,$cell,$0)
end
bench()

using Setfield
a= AnnABase.Parameters()

cell= AnnABase.Cell(a,mode=:precondition);
a=rand(40);
b=collect(range(0,1,length=40));
c=rand(40);
@btime cell.ndim.G($a,$b,$0)
@code_warntype cell.ndim.G(a,b,10)
@code_warntype cell.ndim.R(c,a,b)
@btime cell.ndim.R($c,$a,$b)
test = (R,n,p) ->  R = (n*p-1)*(1+1/(n+1*p+1))
@btime test.($c,$a,$b)
u =ones(length(Matrix(cell.M)[:,1]))
using NLSolversBase
using NLsolve
u(ones(2023),ones(2023))
nlsolve(u,ones(2023) )
du =similar(u)
u= AnnABase.initial_conditions(cell)
du = similar(u)
@btime AnnABase.rhs!($du,$u,$0,$cell)
AnnABase.rhs!(du,u,0,cell)
u
Threads.nthreads()
cell.parameters.τᵢ
size(cell.Jac)
u(cell.Jac,ones(2032))
cell.Jac
(Matrix(cell.M))
(4)*cell.g.N+4+cell.g.Nₑ*2+cell.g.Nₕ

a= AnnABase.Parameters(N=400)
a=setproperties(a,N=300)
cell= AnnABase.Cell(a; op_flavor=:sparse,mode=:oc)
u= AnnABase.initial_conditions(cell)
Mdiag =[sum(cell.M[i,:]) for i in 1:size(cell.M)[1]]
using BenchmarkTools
@btime AnnABase.rhs!($du,$u,$cell,$0)
