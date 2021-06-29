struct Rhs{P,NP,G,OP,M,A,AT,AE,AH,AM,AEM,AHM,AP,AEP,AHP,GEN} <: Function
    parameters::P
    ndim::NP    # ndim parameters
    g::G        # the grid object
    o::OP       # the operators
    mode::M     # the mode

    mE::A   # neagtive electric field
    mEₑ::AE  # negative electric field in ETL
    mEₕ::AH  # negative electric field in HTL
    FP::A   # negative anion vavancy flux
    cd::AM   # charge density
    cdₑ::AEM  # charge density in ETL
    cdₕ::AHM  # charge density in HTL
    fn::A   # electron current
    fnₑ::AE  # electron current in ETL
    fp::A   # hole current
    fpₕ::AH  # hole current in HTL
    GRu::A   # generation recombination sum
    GRt::AT
    G::GEN    # generation rate of 1 sun 
  
    Buff_N::A
    Buff_Nₕ::AH
    Buff_Nₑ::AE
    Buff_N₋₁::AM
    Buff_Nₕ₋₁::AHM
    Buff_Nₑ₋₁::AEM

    P::AP
    ϕ::AP
    n::AP
    p::AP
    ϕₑ::AEP
    nₑ::AEP
    ϕₕ::AHP
    pₕ::AHP
end

function Rhs(parameters, g::Grid, ndim::NodimParameters, op::Operators, mode::Symbol, numtype=Float64)


    d = OrderedDict{Symbol,Any}()
    for key in fieldnames(Rhs)
        d[key] = missing
    end

    
    if numtype <: Union{Float64}
        cache = (u, n) -> DiffEqBase.dualcache(u, n)
    else
        cache = (u, n) -> lcache(length(u))
    end

    N = Val{ForwardDiff.pickchunksize(length(g))}
    d[:parameters] = parameters
    d[:ndim] = ndim
    d[:g] = g
    d[:o] = op
    d[:mode] = mode

    d[:mE]  = cache(zeros(g.N), N)
    d[:mEₑ] = cache(zeros(g.Nₑ), N)
    d[:mEₕ] = cache(zeros(g.Nₕ), N)

    d[:cd]  = cache(zeros(g.N - 1), N)
    d[:cdₑ] = cache(zeros(g.Nₑ - 1), N)
    d[:cdₕ] = cache(zeros(g.Nₕ - 1), N)
    d[:FP]  = cache(zeros(g.N), N)

    d[:fn]  = cache(zeros(g.N), N)
    d[:fp]  = cache(zeros(g.N), N)
    d[:fnₑ] = cache(zeros(g.Nₑ), N)
    d[:fpₕ] = cache(zeros(g.Nₕ), N)
    d[:GRu]  = cache(zeros(g.N), N)
    d[:GRt]  = cache(zeros(g.N), Val{1})
    G = zeros(g.N)
    mul!(G, op.𝕴, g.x)
    ndim.G(G, G, missing)
    d[:G] = G   
    d[:Buff_N]  = cache(zeros(g.N), N)
    d[:Buff_Nₑ]  = cache(zeros(g.Nₑ), N)
    d[:Buff_Nₕ]  = cache(zeros(g.Nₕ), N)
    d[:Buff_N₋₁]  = cache(zeros(g.N - 1), N)
    d[:Buff_Nₑ₋₁]  = cache(zeros(g.Nₑ - 1), N)
    d[:Buff_Nₕ₋₁]  = cache(zeros(g.Nₕ - 1), N)

    d[:ϕ]   = cache(zeros(g.N + 1), N)
    d[:P]   = cache(zeros(g.N + 1), N)
    d[:n]   = cache(zeros(g.N + 1), N)
    d[:p]   = cache(zeros(g.N + 1), N)
    d[:ϕₑ]  = cache(zeros(g.Nₑ + 1), N)
    d[:nₑ]  = cache(zeros(g.Nₑ + 1), N)
    d[:ϕₕ]  = cache(zeros(g.Nₕ + 1), N)
    d[:pₕ]  = cache(zeros(g.Nₕ + 1), N)

    Rhs(collect(values(d))...)
end

function (rhs!::Rhs)(du, u, pr, t)
    δ   = rhs!.ndim.δ
    χ   = rhs!.ndim.χ
    ϰ   = rhs!.ndim.ϰ
    λ   = rhs!.ndim.λ
    λ²  = rhs!.ndim.λ²
    λₑ² = rhs!.ndim.λₑ²
    λₕ² = rhs!.ndim.λₕ²
    κₙ  = rhs!.ndim.κₙ
    κₚ  = rhs!.ndim.κₚ
    κₑ  = rhs!.ndim.κₑ
    κₕ  = rhs!.ndim.κₕ
    rₑ  = rhs!.ndim.rₑ
    rₕ  = rhs!.ndim.rₕ

    N = rhs!.g.N
    Nₑ = rhs!.g.Nₑ
    Nₕ = rhs!.g.Nₕ

    P   = DiffEqBase.get_tmp(rhs!.P, u)
    ϕ   = DiffEqBase.get_tmp(rhs!.ϕ, u)
    n   = DiffEqBase.get_tmp(rhs!.n, u)
    p   =  DiffEqBase.get_tmp(rhs!.p, u)
    @avx for i in 1:N + 1
        P[i]   = u[i]
        ϕ[i]    = u[N + 1 + i]
        n[i]    = u[2 * N + 2 + i]
        p[i]    = u[3 * N + 3 + i]
    end

    ϕₑ  = DiffEqBase.get_tmp(rhs!.ϕₑ, u)
    nₑ  = DiffEqBase.get_tmp(rhs!.nₑ, u)

    @avx for i in 1:rhs!.g.Nₑ
        ϕₑ[i]   = u[4 * N + 4 + i]
        nₑ[i]   = u[4 * N + Nₑ + 4 + i]
    end
    ϕₑ[Nₑ + 1] = ϕ[1]
    nₑ[Nₑ + 1] = n[1] / rhs!.ndim.kₑ

    ϕₕ  = DiffEqBase.get_tmp(rhs!.ϕₕ, u)
    pₕ  = DiffEqBase.get_tmp(rhs!.pₕ, u)
    ϕₕ[1] = ϕ[N + 1]
    pₕ[1] = p[N + 1] / rhs!.ndim.kₕ
    # @inbounds @simd
    @avx for i in 1:rhs!.g.Nₕ
        ϕₕ[i + 1]  = u[4 * N + 2 * Nₑ + 4 + i]
        pₕ[i + 1]  = u[4 * N + 2 * Nₑ + Nₕ + 4 + i]
    end

    mE = DiffEqBase.get_tmp(rhs!.mE, u)
    mul!(mE, rhs!.o.𝔇, ϕ)
    mEₑ = DiffEqBase.get_tmp(rhs!.mEₑ, u)
    mul!(mEₑ, rhs!.o.𝔇ₑ, ϕₑ)
    mEₕ = DiffEqBase.get_tmp(rhs!.mEₕ, u)
    mul!(mEₕ, rhs!.o.𝔇ₕ, ϕₕ)

    Buff_N = DiffEqBase.get_tmp(rhs!.Buff_N, u)
    Buff_Nₑ = DiffEqBase.get_tmp(rhs!.Buff_Nₑ, u)
    Buff_Nₕ = DiffEqBase.get_tmp(rhs!.Buff_Nₕ, u)

    if rhs!.parameters.freeze_ions
        FP = DiffEqBase.get_tmp(rhs!.FP, u) .* 0
    else
        FP  = DiffEqBase.get_tmp(rhs!.FP, u)
        mul!(FP, rhs!.o.𝕴, P); FP .= FP .* mE;
        mul!(Buff_N, rhs!.o.𝔇, P)
        FP .= λ .* (Buff_N .+ FP)
    end

    fn  = DiffEqBase.get_tmp(rhs!.fn, u)
    mul!(fn, rhs!.o.𝕴, n); fn .= fn .* mE;
    mul!(Buff_N, rhs!.o.𝔇, n)
    fn .= κₙ .* (Buff_N .- fn)

    fp  = DiffEqBase.get_tmp(rhs!.fp, u)
    mul!(fp, rhs!.o.𝕴, p); fp .= fp .* mE;
    mul!(Buff_N, rhs!.o.𝔇, p)
    fp .= κₚ .* (Buff_N .+ fp)

    fnₑ  = DiffEqBase.get_tmp(rhs!.fnₑ, u)
    mul!(fnₑ, rhs!.o.𝕴ₑ, nₑ); fnₑ .= fnₑ .* mEₑ;
    mul!(Buff_Nₑ, rhs!.o.𝔇ₑ, nₑ)
    fnₑ .= κₑ .* (Buff_Nₑ .- fnₑ)

    fpₕ  = DiffEqBase.get_tmp(rhs!.fpₕ, pₕ)
    mul!(fpₕ, rhs!.o.𝕴ₕ, pₕ); fpₕ .= fpₕ .* mEₕ;
    mul!(Buff_Nₕ, rhs!.o.𝔇ₕ, pₕ)
    fpₕ .= κₕ .* (Buff_Nₕ .+ fpₕ)

    cd  = DiffEqBase.get_tmp(rhs!.cd, u)
    
    cd_buff = DiffEqBase.get_tmp(rhs!.Buff_N₋₁, u)
    mul!(cd_buff, rhs!.o.𝔏, P)
    cd .= rhs!.g.NN .- cd_buff .- ϰ
    mul!(cd_buff, rhs!.o.𝔏, n)
    cd .= cd .+ δ .* cd_buff
    mul!(cd_buff, rhs!.o.𝔏, p)
    cd .= cd .- δ .* χ .* cd_buff

    cdₑ  = DiffEqBase.get_tmp(rhs!.cdₑ, u)
    mul!(cdₑ, rhs!.o.𝔏ₑ, nₑ)
    cdₑ .= cdₑ .- rhs!.g.ddE

    cdₕ  = DiffEqBase.get_tmp(rhs!.cdₕ, u)
    mul!(cdₕ, rhs!.o.𝔏ₕ, pₕ)
    cdₕ .= rhs!.g.ddH .- cdₕ

    if t isa ForwardDiff.Dual
        GR = DiffEqBase.get_tmp(rhs!.GRt, t)
    else
        GR = DiffEqBase.get_tmp(rhs!.GRu, u)
    end
    mul!(Buff_N, rhs!.o.𝕴, n)   # Buff_N now contains the interpolatet n
    mul!(GR, rhs!.o.𝕴, p)       # GR contains interpolated p
    rhs!.ndim.R(GR, Buff_N, GR) # GR  contains now the recobination rate

   #     mul!(Buff_N, rhs!.o.𝕴, rhs!.g.x) # Buff_N ontains the interpolated x
    #    rhs!.ndim.G(Buff_N, Buff_N, t) # Buff_N contains generation rate

    # !!! her we have allocs !!!! but alternatives lead to Forward diff erro, or seem to be 25% slower s#####
    l =  rhs!.ndim.G.light(t)

    # This is a dirty hack to buffer the timegradient of genereration function 
    # temporarely in the return du vector
    GR   .= muladd(rhs!.G, l, -GR)
  
    @avx for i in 1:N - 1
        du[N + 2 + i] = mE[i + 1] - mE[i] - cd[i] / λ²;
    end
    @inbounds for i in 1:N - 1
        if !rhs!.parameters.freeze_ions
            du[i + 1] = FP[i + 1] - FP[i];
        else
            du[i + 1] = 0
        end
        genrec = (rhs!.g.dx[i + 1] * GR[i + 1] + rhs!.g.dx[i] * GR[i]) / 2
       
        du[2 * N + 3 + i] = fn[i + 1] - fn[i] + genrec;
        du[3 * N + 4 + i] = fp[i + 1] - fp[i] + genrec;
    end
    # P BC
    du[1] = FP[1];      # no Ion flux at bc
    du[N + 1] = -FP[N];   # no Ion flux at bc

    # ϕ BC
    du[N + 2] = mE[1] - rₑ * mEₑ[end] - rhs!.g.dx[1] * (1 / 2 - P[1] / 3 - P[2] / 6 + δ * (n[1] / 3 + n[2] / 6 - χ * (p[1] / 3 + p[2] / 6))) / λ² - rₑ * rhs!.g.dxₑ[end] * (nₑ[end - 1] / 6 + nₑ[end] / 3 - 1 / 2) / λₑ²; # continuity
    du[2 * N + 2] = rₕ * mEₕ[1] - mE[end] - rhs!.g.dx[end] * (1 / 2 - P[end - 1] / 6 - P[end] / 3 + δ * (n[end - 1] / 6 + n[end] / 3 - χ * (p[end - 1] / 6 + p[end] / 3))) / λ² - rₕ * rhs!.g.dxₕ[1] * (1 / 2 - pₕ[1] / 3 - pₕ[2] / 6) / λₕ²; # continuity

    # n BC
    du[2 * N + 3] = fn[1] - fnₑ[end] + (rhs!.g.dx[1] * GR[1]) / 2 - rhs!.ndim.Rₗ(nₑ[end], p[1]) # continuity
    du[3 * N + 3] = -fn[end] + rhs!.g.dx[N] * GR[N] / 2 - rhs!.ndim.Rᵣ(n[N + 1], pₕ[1])
#    println(du[3*N+3])

    # p BC
    du[3 * N + 4] = fp[1] + rhs!.g.dx[1] * GR[1] / 2 - rhs!.ndim.Rₗ(nₑ[end], p[1])
    du[4 * N + 4] = fpₕ[1] - fp[end] + (rhs!.g.dx[end] * GR[end]) / 2 - rhs!.ndim.Rᵣ(n[N + 1], pₕ[1]); # continuity
    
    
    ### ETM ###
    du[4 * N + 5] = ϕₑ[1]# -rhs!.ndim.V(t);
    du[4 * N + Nₑ + 5] = nₑ[1] - 1;
    @avx for i in 1:Nₑ - 1
        du[4 * N + 5 + i] = mEₑ[i + 1] - mEₑ[i] - cdₑ[i] / λₑ²;
        du[4 * N + Nₑ + 5 + i] = fnₑ[i + 1] - fnₑ[i]
    end
    ### HTM ###
   # du[4*N+2*Nₑ+Nₕ+4] = ϕₕ[end] + rhs!.ndim.Vbi - rhs!.ndim.V(t);  
    du[4 * N + 2 * Nₑ + 2 * Nₕ + 4] = pₕ[end] - 1;
    @avx for i in 1:Nₕ - 1
        du[4 * N + 2 * Nₑ + 4 + i] = mEₕ[i + 1] - mEₕ[i] - cdₕ[i] / λₕ²;
        du[4 * N + 2 * Nₑ + Nₕ + 4 + i] = fpₕ[i + 1] - fpₕ[i];
    end
    # Perform any additional step requested by the optional input argument rhs!.mode

    Vr = rhs!.ndim.V(t)
    σₛₕ = rhs!.ndim.σₛₕ
    σₛ = rhs!.ndim.σₛ(t)
    
    if rhs!.mode == :cc  # closed circuit is default
        ϕ_end = Vr / ( σₛₕ / σₛ + 1) + fnₑ[1] / ( σₛₕ + σₛ)
        # hackish fallback to fix initialisation
        if σₛₕ ≈ 0 && σₛ ≈ 0
            ϕ_end = ϕₕ[end] + rhs!.ndim.Vbi
        end
        du[4 * N + 2 * Nₑ + Nₕ + 4] = ϕₕ[end] + rhs!.ndim.Vbi - ϕ_end
    elseif rhs!.mode == :oc # open circuit
        ϕend = ϕₕ[end] + rhs!.ndim.Vbi
        du[4 * N + 5] = (fnₑ[1] - (ϕend) * σₛₕ - ( ϕend - Vr) * σₛ);  # shunt und series resist flux
        du[4 * N + 2 * Nₑ + Nₕ + 4] =  ϕₑ[1] # Potential reference at etl contact
        # du[N+1] = integrate(rhs!.g.x,P)-1;
        # du[N+1] = integrate(rhs!.g.x,P)-1)^4;
    else
        error("simulation mode $(rhs!.mode) is not recongnized")
    end
  
    return nothing
end


struct LCache{T,I}
    du::T
    length::I
end

lcache(n::Int) = LCache(Dict{DataType,AbstractArray}(), n)

function DiffEqBase.get_tmp(lc::LCache, u::T)::T where {T <: AbstractArray} 
    ltype = eltype(u)
    !(haskey(lc.du, ltype)) && (lc.du[ltype] = u[1:lc.length])
    return lc.du[ltype]
end

function DiffEqBase.get_tmp(lc::LCache, u::T)::Array{T,1} where {T <: Number} 
    
    !(haskey(lc.du, T)) && (lc.du[T] = Array{T,1}(undef, lc.length))
    return lc.du[T]
end
