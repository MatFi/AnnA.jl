"""
Taken and modified from the IterativeSolvers.jl package
"""

export idrs, idrs!
import Base: iterate
using Random, UnPack,LinearAlgebra

"""
    idrs(A, b; s = 8, kwargs...) -> x, [history]
Same as [`idrs!`](@ref), but allocates a solution vector `x` initialized with zeros.
"""
idrs(A, b; kwargs...) = idrs!(zerox(A,b), A, b; kwargs...)

"""
    idrs!(x, A, b; s = 8, kwargs...) -> x, [history]
Solve the problem ``Ax = b`` approximately with IDR(s), where `s` is the dimension of the
shadow space.
# Arguments
- `x`: Initial guess, will be updated in-place;
- `A`: linear operator;
- `b`: right-hand side.
## Keywords
- `s::Integer = 8`: dimension of the shadow space;
- `abstol::Real = zero(real(eltype(b)))`,
  `reltol::Real = sqrt(eps(real(eltype(b))))`: absolute and relative
  tolerance for the stopping condition
  `|r_k| ≤ max(reltol * |r_0|, abstol)`, where `r_k = A * x_k - b`
  is the residual in the `k`th iteration;
- `maxiter::Int = size(A, 2)`: maximum number of iterations;
- `log::Bool`: keep track of the residual norm in each iteration;
- `verbose::Bool`: print convergence information during the iterations.
# Return values
**if `log` is `false`**
- `x`: approximate solution.
**if `log` is `true`**
- `x`: approximate solution;
- `history`: convergence history.
"""
function idrs!(x, A, b;
               s = 8,
               abstol::Real = zero(real(eltype(b))),
               reltol::Real = sqrt(eps(real(eltype(b)))),
               maxiter=size(A, 2),
               log::Bool=false,
               kwargs...)
    history = ConvergenceHistory(partial=!log)
    history[:abstol] = abstol
    history[:reltol] = reltol
    log && reserve!(history, :resnorm, maxiter)
    idrs_method!(history, x, A, b, s, abstol, reltol, maxiter; kwargs...)
    log && shrink!(history)
    log ? (x, history) : x
end


mutable struct IDRSIterable{UT,AT,US,MT,CS,T,EUT,S,O}
    X::UT
    A::AT
    C::UT
    R::UT
    Z::UT
    
    P::US
    U::US
    G::US
    Q::UT
    V::UT
    M::MT
    f::CS
    c::CS

    maxiter::T
    normR::EUT
    tol::EUT
    s::S
    om::O
end
function idrs_iterable!(x, A, b,args...;
                        s = 5,
                        abstol::Real =1e-6, #zero(real(eltype(b))),
                        reltol::Real =1e-6,# sqrt(eps(real(eltype(b)))),
                        maxiter=100,
                        log::Bool=false,
                        kwargs...)
    
    zb = zero(b)
    R= b-A*x
    T=typeof(b)
    x.=one.(x)
    IDRSIterable(
        x,A,b,
        R,
        zb,
        T[rand!(copy(b)) for k in 1:s],
        T[copy(zb) for k in 1:s],
        T[copy(zb) for k in 1:s],
        copy(zb),
        copy(zb),
        Matrix{eltype(b)}(I,s,s),
        zeros(eltype(b),s),
        zeros(eltype(b),s),  
        maxiter,
        norm(R),
        max(reltol * norm(R), abstol),
        s,
        eltype(b)(1)
    )
    #idrs_method!(history, x, A, b, s, abstol, reltol, maxiter; kwargs...)
    #log && shrink!(history)
    #log ? (x, history) : x
end

start(::IDRSIterable) = 0
converged(g::IDRSIterable) = g.normR < g.tol
done(g::IDRSIterable, iteration::Int) = iteration ≥ g.maxiter || converged(g)
#########################
# Method Implementation #
#########################

@inline function omega(t, s)
    angle = sqrt(2.)/2
    ns = norm(s)
    nt = norm(t)
    ts = dot(t,s)
    rho = abs(ts/(nt*ns))
    om = ts/(nt*nt)
    if rho < angle
        om = om*convert(typeof(om),angle)/rho
    end
    om
end

function iterate(g::IDRSIterable, iteration::Int=start(g))
    if done(g, iteration) return nothing end
    @unpack  X,A,C,R,Z,P,U, G,Q, V,M,f,c,normR,tol,s,om = g
    
   
    for i in 1:s
        f[i] = dot(P[i], R)
    end
   
    for k in 1:s
        # Solve small system and make v orthogonal to P
        
        c =  Matrix(LowerTriangular(M[k:s,k:s]))\f[k:s]
        V .= c[1] .* G[k]
        Q .= c[1] .* U[k]

        for i = k+1:s
            V .+= c[i-k+1] .* G[i]
            Q .+= c[i-k+1] .* U[i]
        end

        # Compute new U[:,k] and G[:,k], G[:,k] is in space G_j
        V .= R .- V

        U[k] .= Q .+ om .* V
        mul!(G[k], A, U[k])

        # Bi-orthogonalise the new basis vectors

        for i in 1:k-1
            alpha = dot(P[i],G[k])/M[i,i]
            G[k] .-= alpha .* G[i]
            U[k] .-= alpha .* U[i]
        end

        # New column of M = P'*G  (first k-1 entries are zero)

        for i in k:s
            M[i,k] = dot(P[i],G[k])
        end

        #  Make r orthogonal to q_i, i = 1..k

        beta = f[k]/M[k,k]
        R .-= beta .* G[k]
        X .+= beta .* U[k]

        g.normR = norm(R)
       
        
        #verbose && 
     #   @printf("%3d\t%1.2e\n",iter,normR)
      #  if normR < tol || iteration == maxiter
      #      setconv(log, 0<=normR<tol)
      #      return X
      #  end 
        if k < s
            f[k+1:s] .-=  beta*M[k+1:s,k]
        end
        
    end

    # Now we have sufficient vectors in G_j to compute residual in G_j+1
    # Note: r is already perpendicular to P so v = r
    copyto!(V, R)
    
    mul!(Q, A, V)
    g.om = omega(Q, R)
    R .-= om .* Q
    X .+= om .* V

    g.normR = norm(R)

   # nextiter!(log, mvps=1)
   # push!(log, :resnorm, normR)
    nothing, iteration+1
end

