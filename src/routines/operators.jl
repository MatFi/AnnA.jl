
#This stuct can be rewritten in a nicer way
abstract type AbstractOperators end

struct Operators{A, L, D, AH, LH, DH,AE, LE, DE} <: AbstractOperators
    š“::A
    š“ā::AE
    š“ā::AH
    š::L
    šā::LE
    šā::LH
    š::D
    šā::DE
    šā::DH
end

"""
    Operators(g::Grid;dense::Bool=true)
Compute the operators based on a given `Grid`. Use the `dense` keyword to switch between sparse ord dens representation
"""
function Operators(g::Grid;flavor::Symbol=:dense,numtype=Float64)
    #for n in fieldnames(Grid)
    #    @eval $n = $g.$(n)
    #end
    p = unpac_struct(g)
    N=p[:N]
    Nā=p[:Nā]
    Nā=p[:Nā]
    dx=p[:dx]
    dxā=p[:dxā]
    dxā=p[:dxā]

    #pack filednames to dict
    d=OrderedDict{Symbol,Any}()
    for key in fieldnames(Operators)
        d[key]=missing
    end

    d[:š“]  = (Tridiagonal(N+1,0,1,1)/2)[1:N,1:N+1]
    d[:š“ā] = (Tridiagonal(Nā+1,0,1,1)/2)[1:Nā,1:Nā+1]
    d[:š“ā] = (Tridiagonal(Nā+1,0,1,1)/2)[1:Nā,1:Nā+1]

    d[:š]  = (Tridiagonal(vcat(dx[1:end-1]/6,0),vcat(0,(dx[1:end-1]+dx[2:end])/3,0),vcat(0,dx[2:end]/6)))
    d[:š]  = d[:š][2:N,1:N+1]
    d[:šā] = (Tridiagonal(vcat(dxā[1:end-1]/6,0),vcat(0,(dxā[1:end-1]+dxā[2:end])/3,0),vcat(0,dxā[2:end]/6)))
    d[:šā] =  d[:šā][2:Nā,1:Nā+1]
    d[:šā] = (Tridiagonal(vcat(dxā[1:end-1]/6,0),vcat(0,(dxā[1:end-1]+dxā[2:end])/3,0),vcat(0,dxā[2:end]/6)))
    d[:šā] =  d[:šā][2:Nā,1:Nā+1]
    d[:š]  = Tridiagonal(0 ./dx,vcat(-1 ./dx,0), 1 ./dx)[1:N,1:N+1]
    d[:šā] = Tridiagonal(0 ./dxā,vcat(-1 ./dxā,0), 1 ./dxā)[1:Nā,1:Nā+1]
    d[:šā] = Tridiagonal(0 ./dxā,vcat(-1 ./dxā,0), 1 ./dxā)[1:Nā,1:Nā+1]
    for k in keys(d)
        d[k] = dropzeros(sparse(d[k]))
    end
    if flavor == :matrix_free
        d[:š“]  = InterpolationOperator(size(d[:š“]))
        d[:š“ā] = InterpolationOperator(size(d[:š“ā]))
        d[:š“ā] = InterpolationOperator(size(d[:š“ā]))
        d[:š]  = LinOperator(d[:š])#, size(d[:š]))
        d[:šā] = LinOperator(d[:šā])#, size(d[:šā]))
        d[:šā] = LinOperator(d[:šā])#, size(d[:šā]))
        d[:š]  = DiffOperator(d[:š])#, size(d[:š]))
        d[:šā] = DiffOperator(d[:šā])#, size(d[:šā]))
        d[:šā] = DiffOperator(d[:šā])#, size(d[:šā]))
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
function DiffOperator(dx)
    for i in 1:size(dx,1)
        dx[i,1] = dx[i,i] 
        dx[i,2] = dx[i,i+1]
    end
    DiffOperator(dx,size(dx))

end

function LinearAlgebra.mul!(dfdx,dx::DiffOperator,f)
    @avx for i in 1:dx.size[1]
        dfdx[i] = f[i+1]*dx.dx[i,2]+f[i]*dx.dx[i,1]
    end
    nothing
end


struct InterpolationOperator{S} <:AbstractOperators
    size::S
end
function LinearAlgebra.mul!(interp,o::InterpolationOperator,f)
  #  if o.size[2] != length(f)
  #      error("length of operator is $(o.size[2]), while length(f) is $(length(f)) ")
  #  end
   @avx for i in 1:o.size[1]
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
function LinOperator(dx)
    for i in 1:size(dx,1)
        dx[i,1] = dx[i,i] 
        dx[i,2] = dx[i,i+1]
        dx[i,3] = dx[i,i+2]
    end
    LinOperator(dx,size(dx))

end
function LinearAlgebra.mul!(lin,o::LinOperator,f)
 #   if o.size[2] != length(f)
 #       error("length of operator is $(o.size[2]), while length(f) is $(length(f)) ")
 #   end
    @avx for i in 1:o.size[1]
        lin[i] = f[i]*o.dx[i,1]+f[i+1]*o.dx[i,2]+f[i+2]*o.dx[i,3]
    end
    nothing
end
