using AnnA
using Unitful, UnitfulRecipes
parm = Parameters(light = t -> 1.0,
    vₙₕ= 10u"m/s" ,                 # electron surface recombination vel. at HTM
    vₚₕ= 0.01u"m/s" ,               # hole surface recombination vel. at HTM
    N=200,                          # grid size
    N₀=1e18u"cm^-3"                 # ionic concentration
)
prob_ocvd = OCVDProblem(
    parm,       # input parameter set
    50u"s",     # illumination time
    1e5u"s",    # time the decay will be simulated to
    alg_control = AnnA.AlgControl(
        dtmin = 1e-50,
        dt = 1e-4,
        ss_tol= 1e-5,
        maxiters = 5000,
        reltol = 1e-3,
        abstol = 1e-8,
    ),  
)
sol = @time solve(prob_ocvd)

using Plots
plot!(sol.t_decay,sol.V_decay,xscale=:log10,xlims=(1e-8,1e5))
sol.t_decay

td = AnnA.Tridiagonal(9,-1,2,1)
BlockBandedMatrix(td,[5,4],[5,4],(1,1))


l,u = 1,1          # block bandwidths
λ,μ = 1,1          # sub-block bandwidths: the bandwidths of each block
N = M = 8          # number of row/column blocks
rows = [5,5,5]  # block sizes
cols = [5,5,5]
BlockBandedMatrix(Zeros(sum(rows),sum(cols)), rows,cols, (l,u)) # creates a block-banded matrix of zeros
BlockBandedMatrix(Ones(sum(rows),sum(cols)), rows,cols, (l,u)) # creates a block-banded matrix with ones in the non-zero entries
BandedBlockBandedMatrix(I, rows,cols, (l,u), (λ,μ))                         # creates a block-banded  identity matrix

c=AnnA.Cell(parm,mode=:oc)
c.Jac

l,u = 3,3  
λ,μ = 1,1
N=c.g.N
Nₑ=c.g.Nₑ
Nₕ=c.g.Nₕ
rows = [N+1,N+1,N+1,N+1,Nₑ,Nₑ,Nₕ,Nₕ]  # block sizes
cols = [N+1,N+1,N+1,N+1,Nₑ,Nₑ,Nₕ,Nₕ]
bjac =BandedBlockBandedMatrix(c.Jac, rows,cols, (l,u), (λ,μ))                         # creates a block-banded  identity matrix

(bjac .!= c.Jac )|> any 

bicgstabl!(x, A, b, l; kwargs...) -> x, [history]

c.Jac\zeros(size(c.Jac,1))
c=AnnA.Cell(parm,mode=:oc)
A=c.Jac.+I(size(c.Jac,1))
b=ones(size(A,1))

xl=A\b

using AlgebraicMultigrid
ml = ruge_stuben(A)
p = aspreconditioner(ml)

@time _,itr =bicgstabl!(xi, A, b, 5,Pl=p,abstol=1e-6,reltol=1e-6,log=true)
@time _,itr =gmres!(xi, A, b,abstol=1e-6,reltol=1e-6,Pl=p,log=true,maxiter=10000)
@time _,itr =cg!(xi, A, b,abstol=1e-6,reltol=1e-6,Pl=p,log=true,maxiter=10000)
@time _,itr =chebyshev!(xi, A, b,-1000,1000,abstol=1e-6,reltol=1e-6,Pl=p,log=true,maxiter=10000)
@time _,itr =minres!(xi, A, b,abstol=1e-6,reltol=1e-6,log=true,maxiter=10000)
set_precision


xi=zero(b)

xi=Float128.(zero(b)) 
xi=Double64.(zero(b)) 
@time xi =jacobi!(xi, A, b; maxiter=10)
@time xi =gauss_seidel!(xi, A, b; maxiter=10) 
@time idrs!(xi, A, b,s=40,abstol=eps(Double64),reltol=eps(Double64),maxiter=10000)
@time  _,itr =lsmr!(xi, A, b,maxiter=10000,atol =eps(Double64), btol=eps(Double64),conlim=eps(Double64))
@time _,itr =lsqr!(xi, A, b,maxiter=10000,atol =eps(Double64), btol=eps(Double64),conlim=eps(Double64))



A*xi -b
A*xl -b
A*xl -b
using Plots
plot!(itr,yscale=:log10)
@time xl=A\b
b64=Double64.(b)
@time xl=A\b
xi

xi=zero(b)

#include("./src/linsolvers/idrs.jl")
itr = AnnA.idrs_iterable!(xi, A, b;s=40)
for i in itr
    @show itr.R[1]
end
norm(b)
eltype(b)(1)  isa typeof(norm(b))

for _ = itr

end
itr.X
A*itr.X -b
itr.R
