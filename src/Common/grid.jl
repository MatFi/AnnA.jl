struct Grid{A<:Integer,B<:AbstractArray} <: Function
    N::A
    x::B
    dx::B
    Nₑ::A
    xₑ::B
    dxₑ::B
    Nₕ::A
    xₕ::B
    dxₕ::B

    NN::B
    ddE::B
    ddH::B
end

function Grid(p::AbstractParameters)
    #pack filednames to dict
    d = OrderedDict{Symbol,Any}()
    for key in fieldnames(Grid)
        d[key] = missing
    end

    X = 0.2
    wₑ = convert(Float64, p.bₑ / p.b)
    wₕ = convert(Float64, p.bₕ / p.b)

    LD = sqrt(p.kB * p.T * p.εₚ / (p.q^2 * p.N₀))
    λ = convert(Float64, LD / p.b)
    #@show λ
    st = 1.0
    try
        tanfun = st -> λ - (tanh(st * (2 * X - 1)) / tanh(st) + 1) / 2
        st = find_zero(tanfun, 2, atol = 1e-3)
    catch e

    end
    st = p.freeze_ions ? 1 : st # are the ions loked (uniform grid)

    A = b -> (tanh(st * (1 - 2 / p.N)) - (1.0 - b) * tanh(st)) / b


    d[:N] = p.N
    d[:Nₑ] = trunc(Int, round(2 / (1 - atanh(A(wₑ)) / st)))
    d[:Nₕ] = trunc(Int, round(2 / (1 - atanh(A(wₕ)) / st)))


    #println(d[:Nₕ] ," ",d[:Nₑ] )
    #build the grids
    d[:x] = collect(range(0, 1; length = d[:N] + 1))
    d[:xₑ] = collect(range(-wₑ, 0; length = d[:Nₑ] + 1))
    d[:xₕ] = collect(range(1, 1 + wₕ; length = d[:Nₕ] + 1))

    d[:x] = @. (tanh(st * (2 * d[:x] - 1)) / tanh(st) + 1) / 2
    d[:xₑ] = @. wₑ * (tanh(st * (2 * d[:xₑ] / wₑ + 1)) / tanh(st) - 1) / 2
    d[:xₕ] = @. 1 + wₕ * (tanh(st * (2 * (d[:xₕ] - 1) / wₕ - 1)) / tanh(st) + 1) / 2

    d[:dx] = diff(d[:x])
    d[:dxₑ] = diff(d[:xₑ])
    d[:dxₕ] = diff(d[:xₕ])

    d[:xₑ] = d[:xₑ][1:end-1]
    d[:xₕ] = d[:xₕ][2:end]

    d[:NN] = (d[:dx][2:end] + d[:dx][1:end-1]) / 2
    d[:ddE] = (d[:dxₑ][2:end] + d[:dxₑ][1:end-1]) / 2
    d[:ddH] = (d[:dxₕ][2:end] + d[:dxₕ][1:end-1]) / 2
    Grid(collect(values(d))...)
end

"""
    length(g::Grid)
gives the length of the input-array which is ` 4*g.N+4+2*g.Nₑ+2*g.Nₕ`
"""
length(g::Grid) = 4 * g.N + 4 + 2 * g.Nₑ + 2 * g.Nₕ

#lagacy wrapper, better use `decompose(u::Array{Float64,1},g::Grid)`
function (decomp::Grid)(u::Array{Float64,1})
    decompose(u, decomp)
end

decompose(u::Array{Array{Float64,1},1}, g::Grid) = decompose.(u, g)

function decompose(u::Array{Float64,1}, g::Grid)
    N = g.N
    Nₑ = g.Nₑ
    Nₕ = g.Nₕ

    P = Array{Float64}(undef, N + 1)
    ϕ = Array{Float64}(undef, N + 1)
    n = Array{Float64}(undef, N + 1)
    p = Array{Float64}(undef, N + 1)
    @inbounds @simd for i = 1:N+1
        P[i] = u[i]
        ϕ[i] = u[N + 1 + i]
        n[i] = u[2 * N + 2 + i]
        p[i] = u[3 * N + 3 + i]
    end

    ϕₑ = Array{Float64}(undef, Nₑ)
    nₑ = Array{Float64}(undef, Nₑ)
    @inbounds @simd for i = 1:Nₑ
        ϕₑ[i] = u[4 * N + 4 + i]
        nₑ[i] = u[4 * N + Nₑ + 4 + i]
    end

    ϕₕ = Array{Float64}(undef, Nₕ)
    pₕ = Array{Float64}(undef, Nₕ)
    @inbounds @simd for i = 1:Nₕ
        ϕₕ[i] = u[4 * N + 2 * Nₑ + 4 + i]
        pₕ[i] = u[4 * N + 2 * Nₑ + Nₕ + 4 + i]
    end

    [[P], [ϕₑ, ϕ, ϕₕ], [nₑ, n], [p, pₕ]]
end
