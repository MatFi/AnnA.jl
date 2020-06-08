data_ocvd = []
data_iv = []

polyfit(x::Vector, y::Vector, deg::Int) = collect(v ^ p for v in x, p in 0:deg) \ y

function load_raw(path::AbstractString; delim::Char='\t')
    return readdlm(path,delim,Float64, '\n'; comments = true, header = true)
end

function load(path::AbstractString; xrange::Tuple = (), sort::Bool=false,skipnan::Bool=false, delim::Char='\t')
    (data,header) = load_raw(path; delim=delim)
    if skipnan
        p=Array{Bool}(undef,size(data,1),1)
        any!(p, isnan.(data))
        data= data[.!vec(p), :]
    end
    if !isempty(xrange)
        data = data[ustrip(upreferred(xrange[1])) .<= data[:, 1] .<= ustrip(upreferred(xrange[2])), :]
    end

    sort && return data[sortperm(data[:,1]),:]
    return data
end


"""
    _resample(data::Array{T,2},N::Int,::Val{:xlog}) where T <: Number

resamples the data to a log scaled x-axis by appling peacwise quadratic fit arround the 'N' samle points ( this apporach is comparable to loesssmotthing with non-constant stepp size)
"""
function _resample(data::Array{T,2}, N::Int, ::Val{:xlog}) where {T<:Number}

    deg_base = 2  # mximum fit order

    ΔN = size(data)[1] / N
    xmin = log10(minimum(data))
    if isinf(xmin)
        @error "x-values all needs to be positive"
    end
    xmax = log10(maximum(data[end, 1]))
    Δx = (xmax - xmin) / N
    xi = i -> 10.0^(xmin + i * Δx)

    sampled = Array{T}(undef, N, 2)
    for i = 1:N
        x12 = (xi(i) + xi(i - 1)) / 2

        imin = argmin(abs.(data[:, 1] .- xi(i - 1)))        # find the index of min bound
        imin = data[imin, 1] < xi(i - 1) ? imin : imin - 1  # extend the min bound to the datapoint with lower value
        # treat imax vice versa
        imax = argmin(abs.(data[:, 1] .- xi(i)))
        imax = data[imax, 1] > xi(i - 1) ? imax : imax + 1

        # limit the index to stay in bounds
        imin = clamp(imin,1,size(data)[1]-1)

        imax = clamp(imax,2,size(data)[1])
        part = data[imin:imax, :]

        #adaptivly reduce polynomial order when sample points are low ( bonderies )
        deg = size(part)[1] > deg_base ? deg_base : size(part)[1]

        coefs = polyfit(part[:, 1], part[:, 2], deg)
        #sampled[i,:]=[x12,sum(coefs .*collect(x12 ^ p for p in 0:deg) )]
        sampled[i,:] = [x12 sum(coefs .* collect(x12^p for p = 0:deg))]

        # include the ends
        if i == N
            sampled[i,:] =  [xi(i) sum(coefs .* collect(xi(i)^p for p = 0:deg))]
        elseif i ==1
            sampled[i,:] =  [xi(i-1) sum(coefs .* collect(xi(i-1)^p for p = 0:deg))]
        end
    end
    return sampled
end

function _resample(data::Array{T,2}, N::Int, ::Val{:ylog}) where {T<:Number}
    data2=similar(data)
    data2[:,1]=data[:,2]
    data2[:,2]=data[:,1]
    ret = Array{Float64}(undef,N,2)
    data2 = resample(data2, N, :xlog) 
    ret[:,1]=data2[:,2]
    ret[:,2]=data2[:,1]
    return ret
end

function resample(data,N::Int,type::Symbol)
    _resample(data,N,Val(type))
end
