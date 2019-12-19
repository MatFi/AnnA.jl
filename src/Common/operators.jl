
#This stuct can be rewritten in a nicer way
abstract type AbstractOperators end

struct Operators{A, L, D, AH, LH, DH,AE, LE, DE} <: AbstractOperators
    𝕴::A
    𝕴ₑ::AE
    𝕴ₕ::AH
    𝔏::L
    𝔏ₑ::LE
    𝔏ₕ::LH
    𝔇::D
    𝔇ₑ::DE
    𝔇ₕ::DH
end

"""
    Operators(g::Grid;dense::Bool=true)
Compute the operators based on a given `Grid`. Use the `dense` keyword to switch between sparse ord dens representation
"""
function Operators(g::Grid;flavor::Symbol=:dense)
    for n in fieldnames(Grid)
        @eval $n = $g.$(n)
    end

    #pack filednames to dict
    d=OrderedDict{Symbol,Any}()
    for key in fieldnames(Operators)
        d[key]=missing
    end

    d[:𝕴]  = (Tridiagonal(N+1,0,1,1)/2)[1:N,1:N+1]
    d[:𝕴ₑ] = (Tridiagonal(Nₑ+1,0,1,1)/2)[1:Nₑ,1:Nₑ+1]
    d[:𝕴ₕ] = (Tridiagonal(Nₕ+1,0,1,1)/2)[1:Nₕ,1:Nₕ+1]

    d[:𝔏]  = (Tridiagonal(vcat(dx[1:end-1]/6,0),vcat(0,(dx[1:end-1]+dx[2:end])/3,0),vcat(0,dx[2:end]/6)))
    d[:𝔏]  = d[:𝔏][2:N,1:N+1]
    d[:𝔏ₑ] = (Tridiagonal(vcat(dxₑ[1:end-1]/6,0),vcat(0,(dxₑ[1:end-1]+dxₑ[2:end])/3,0),vcat(0,dxₑ[2:end]/6)))
    d[:𝔏ₑ] =  d[:𝔏ₑ][2:Nₑ,1:Nₑ+1]
    d[:𝔏ₕ] = (Tridiagonal(vcat(dxₕ[1:end-1]/6,0),vcat(0,(dxₕ[1:end-1]+dxₕ[2:end])/3,0),vcat(0,dxₕ[2:end]/6)))
    d[:𝔏ₕ] =  d[:𝔏ₕ][2:Nₕ,1:Nₕ+1]
    d[:𝔇]  = Tridiagonal(0 ./dx,vcat(-1 ./dx,0), 1 ./dx)[1:N,1:N+1]
    d[:𝔇ₑ] = Tridiagonal(0 ./dxₑ,vcat(-1 ./dxₑ,0), 1 ./dxₑ)[1:Nₑ,1:Nₑ+1]
    d[:𝔇ₕ] = Tridiagonal(0 ./dxₕ,vcat(-1 ./dxₕ,0), 1 ./dxₕ)[1:Nₕ,1:Nₕ+1]

    if flavor == :matrix_free
        d[:𝕴]  = InterpolationOperator(size(d[:𝕴]))
        d[:𝕴ₑ] = InterpolationOperator(size(d[:𝕴ₑ]))
        d[:𝕴ₕ] = InterpolationOperator(size(d[:𝕴ₕ]))
        d[:𝔏]  = LinOperator(d[:𝔏], size(d[:𝔏]))
        d[:𝔏ₑ] = LinOperator(d[:𝔏ₑ], size(d[:𝔏ₑ]))
        d[:𝔏ₕ] = LinOperator(d[:𝔏ₕ], size(d[:𝔏ₕ]))
        d[:𝔇]  = DiffOperator(d[:𝔇], size(d[:𝔇]))
        d[:𝔇ₑ] = DiffOperator(d[:𝔇ₑ], size(d[:𝔇ₑ]))
        d[:𝔇ₕ] = DiffOperator(d[:𝔇ₕ], size(d[:𝔇ₕ]))
        Operators(collect(values(d))...)
    elseif flavor ==:dense
        Operators(Matrix.(collect(values(d)))...)
    elseif flavor ==:sparse
        Operators(collect(values(d))...)
    else
        error("unkown operator flavor $flavor")
    end
end



struct DiffOperator{T,S} <:AbstractOperators
    dx::T
    size::S
end
function LinearAlgebra.mul!(dfdx,dx::DiffOperator,f)
    if dx.size[2] != length(f)
        error("length of operator is $(dx.size[2]), while length(f) is $(length(f)) ")
    end

    for i in 1:dx.size[1]
        dfdx[i] = f[i+1]*dx.dx[i,i+1]+f[i]*dx.dx[i,i]
    end
    nothing
end


struct InterpolationOperator{S} <:AbstractOperators
    size::S
end
function LinearAlgebra.mul!(interp,o::InterpolationOperator,f)
    if o.size[2] != length(f)
        error("length of operator is $(o.size[2]), while length(f) is $(length(f)) ")
    end
    for i in 1:o.size[1]
        interp[i] = (f[i+1]+f[i])/2.0
    end
    nothing
end


"""
    struct LinOperator{T, S} <: AbstractOperators

"""
struct LinOperator{T,S} <: AbstractOperators
    dx::T
    size::S
end
function LinearAlgebra.mul!(lin,o::LinOperator,f)
    if o.size[2] != length(f)
        error("length of operator is $(o.size[2]), while length(f) is $(length(f)) ")
    end
    for i in 1:o.size[1]
        lin[i] = f[i]*o.dx[i,i]+f[i+1]*o.dx[i,i+1]+f[i+2]*o.dx[i,i+2]
    end
    nothing
end
