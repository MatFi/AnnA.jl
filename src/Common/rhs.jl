struct Rhs{A,AE,AH,AM,AEM,AHM, AP, AEP, AHP} <: Function
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
    GR::A   # generation recombination sum

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


function Rhs(g::Grid)
    #pack filednames to dict
    d=OrderedDict{Symbol,Any}()
    for key in fieldnames(Rhs)
        d[key]=missing
    end

    N = Val{ForwardDiff.pickchunksize(length(g))}
    d[:mE]  = DiffEqBase.dualcache(zeros(g.N),N)
    d[:mEₑ] = DiffEqBase.dualcache(zeros(g.Nₑ),N)
    d[:mEₕ] = DiffEqBase.dualcache(zeros(g.Nₕ),N)

    d[:cd]  = DiffEqBase.dualcache(zeros(g.N-1),N)
    d[:cdₑ] = DiffEqBase.dualcache(zeros(g.Nₑ-1),N)
    d[:cdₕ] = DiffEqBase.dualcache(zeros(g.Nₕ-1),N)

    d[:FP]  = DiffEqBase.dualcache(zeros(g.N),N)
    d[:fn]  = DiffEqBase.dualcache(zeros(g.N),N)
    d[:fp]  = DiffEqBase.dualcache(zeros(g.N),N)
    d[:fnₑ] = DiffEqBase.dualcache(zeros(g.Nₑ),N)
    d[:fpₕ] = DiffEqBase.dualcache(zeros(g.Nₕ),N)
    d[:GR]  = DiffEqBase.dualcache(zeros(g.N),N)

    d[:Buff_N]  = DiffEqBase.dualcache(zeros(g.N),N)
    d[:Buff_Nₑ]  = DiffEqBase.dualcache(zeros(g.Nₑ),N)
    d[:Buff_Nₕ]  = DiffEqBase.dualcache(zeros(g.Nₕ),N)
    d[:Buff_N₋₁]  = DiffEqBase.dualcache(zeros(g.N-1),N)
    d[:Buff_Nₑ₋₁]  = DiffEqBase.dualcache(zeros(g.Nₑ-1),N)
    d[:Buff_Nₕ₋₁]  = DiffEqBase.dualcache(zeros(g.Nₕ-1),N)

    d[:ϕ]   = DiffEqBase.dualcache(zeros(g.N+1),N)
    d[:P]   = DiffEqBase.dualcache(zeros(g.N+1),N)
    d[:n]   = DiffEqBase.dualcache(zeros(g.N+1),N)
    d[:p]   = DiffEqBase.dualcache(zeros(g.N+1),N)
    d[:ϕₑ]  = DiffEqBase.dualcache(zeros(g.Nₑ+1),N)
    d[:nₑ]  = DiffEqBase.dualcache(zeros(g.Nₑ+1),N)
    d[:ϕₕ]  = DiffEqBase.dualcache(zeros(g.Nₕ+1),N)
    d[:pₕ]  = DiffEqBase.dualcache(zeros(g.Nₕ+1),N)

    Rhs(collect(values(d))...)
end

function (rhs!::Rhs)(du,u,c::Cell,t)

    δ   = c.ndim.δ
    χ   = c.ndim.χ
    λ   = c.ndim.λ
    λ²  = c.ndim.λ²
    λₑ² = c.ndim.λₑ²
    λₕ² = c.ndim.λₕ²
    κₙ  = c.ndim.κₙ
    κₚ  = c.ndim.κₚ
    κₑ  = c.ndim.κₑ
    κₕ  = c.ndim.κₕ
    rₑ  = c.ndim.rₑ
    rₕ  = c.ndim.rₕ
#=
#Array Partition
    P   = u.x[1]
    ϕ   = u.x[2]
    n   = u.x[3]
    p   = u.x[4]
    ϕₑ  = u.x[5]
    nₑ  = u.x[6]
    ϕₕ  = u.x[7]
    pₕ  = u.x[8]

VectorOfArray
    P   = u[:,1]
    ϕ   = u[:,2]
    n   = u[:,3]
    p   = u[:,4]
    ϕₑ  = u[:,5]
    nₑ  = u[:,6]
    ϕₕ  = u[:,7]
    pₕ  = u[:,8]
=#
    N= c.g.N
    Nₑ=c.g.Nₑ
    Nₕ=c.g.Nₕ

    P   = DiffEqBase.get_tmp(rhs!.P,u)
    ϕ   = DiffEqBase.get_tmp(rhs!.ϕ,u)
    n   = DiffEqBase.get_tmp(rhs!.n,u)
    p   =  DiffEqBase.get_tmp(rhs!.p,u)
    @inbounds @simd for i in 1:N+1
        P[i]   = u[i]
        ϕ[i]    = u[N+1+i]
        n[i]    = u[2*N+2+i]
        p[i]    = u[3*N+3+i]
    end

    ϕₑ  = DiffEqBase.get_tmp(rhs!.ϕₑ,u)
    nₑ  = DiffEqBase.get_tmp(rhs!.nₑ,u)
    ϕₑ[Nₑ+1] = ϕ[1]
    nₑ[Nₑ+1] = n[1]/c.ndim.kₑ
    @inbounds @simd for i in 1:c.g.Nₑ
        ϕₑ[i]   = u[4*N+4+i]
        nₑ[i]   = u[4*N+Nₑ+4+i]
    end
#    ϕₑ[Nₑ] = ϕ[1]  #nasty fix probably wrong Nₑ ind nodim
#    nₑ[Nₑ+1] = n[1]/c.ndim.kₑ

    ϕₕ  = DiffEqBase.get_tmp(rhs!.ϕₕ,u)
    pₕ  = DiffEqBase.get_tmp(rhs!.pₕ,u)
    ϕₕ[1]=ϕ[N+1]
    pₕ[1]=p[N+1]/c.ndim.kₕ
    #@inbounds @simd
    @inbounds @simd for i in 1:c.g.Nₕ
        ϕₕ[i+1]  = u[4*N+2*Nₑ+4+i]
        pₕ[i+1]  = u[4*N+2*Nₑ+Nₕ+4+i]
    end

    mE=DiffEqBase.get_tmp(rhs!.mE,u)
    mul!(mE,c.o.𝔇,ϕ)
    mEₑ=DiffEqBase.get_tmp(rhs!.mEₑ,u)
    mul!(mEₑ,c.o.𝔇ₑ,ϕₑ)
    mEₕ=DiffEqBase.get_tmp(rhs!.mEₕ,u)
    mul!(mEₕ,c.o.𝔇ₕ,ϕₕ)

    Buff_N=DiffEqBase.get_tmp(rhs!.Buff_N,u)
    Buff_Nₑ=DiffEqBase.get_tmp(rhs!.Buff_Nₑ,u)
    Buff_Nₕ=DiffEqBase.get_tmp(rhs!.Buff_Nₕ,u)

    FP  = DiffEqBase.get_tmp(rhs!.FP,u)
        mul!(FP,c.o.𝕴,P); FP .= FP.*mE;
        mul!(Buff_N ,c.o.𝔇,P)
        FP .= λ .* (Buff_N .+ FP)

    fn  = DiffEqBase.get_tmp(rhs!.fn,u)
        mul!(fn,c.o.𝕴,n); fn .= fn.*mE;
        mul!(Buff_N ,c.o.𝔇,n)
        fn .= κₙ.*(Buff_N .- fn)

    fp  = DiffEqBase.get_tmp(rhs!.fp,u)
        mul!(fp,c.o.𝕴,p); fp .= fp.*mE;
        mul!(Buff_N ,c.o.𝔇,p)
        fp .= κₚ.*(Buff_N .+ fp)

    fnₑ  = DiffEqBase.get_tmp(rhs!.fnₑ,u)
        mul!(fnₑ,c.o.𝕴ₑ,nₑ); fnₑ .= fnₑ.*mEₑ;
        mul!(Buff_Nₑ ,c.o.𝔇ₑ,nₑ)
        fnₑ .= κₑ.*(Buff_Nₑ .- fnₑ)

    fpₕ  = DiffEqBase.get_tmp(rhs!.fpₕ,pₕ)
        mul!(fpₕ,c.o.𝕴ₕ,pₕ); fpₕ .= fpₕ.*mEₕ;
        mul!(Buff_Nₕ ,c.o.𝔇ₕ,pₕ)
        fpₕ .= κₕ.*(Buff_Nₕ .+ fpₕ)

    cd  = DiffEqBase.get_tmp(rhs!.cd,u)
    cd_buff = DiffEqBase.get_tmp(rhs!.Buff_N₋₁,u)
        mul!(cd_buff,c.o.𝔏,P)
        cd .= c.g.NN .- cd_buff
        mul!(cd_buff,c.o.𝔏,n)
        cd .= cd .+ δ.*cd_buff
        mul!(cd_buff,c.o.𝔏,p)
        cd .= cd .- δ.*χ.*cd_buff
        #mlab cd = NN-Lo*P+delta*(Lo*n-chi*Lo*p); # charge density

    cdₑ  = DiffEqBase.get_tmp(rhs!.cdₑ,u)
        mul!(cdₑ,c.o.𝔏ₑ,nₑ)
        cdₑ .= cdₑ .- c.g.ddE
        #mlab   cdₑ = LoE*nE-ddE; # charge density in ETL
    cdₕ  = DiffEqBase.get_tmp(rhs!.cdₕ,u)
        mul!(cdₕ,c.o.𝔏ₕ,pₕ)
        cdₕ .= c.g.ddH .- cdₕ
        #mlab cdₕ = ddH-LoH*pH; # charge density in HTL

    GR  = DiffEqBase.get_tmp(rhs!.GR, u)
        mul!(Buff_N, c.o.𝕴, n)   # Buff_N now contains the interpolatet n
        mul!(GR, c.o.𝕴, p)       # GR contains interploated p
        c.ndim.R(GR ,Buff_N, GR) # GR  contains now the recobination rate

        mul!(Buff_N, c.o.𝕴, c.g.x) # Buff_N ontains the interpolated x
        c.ndim.G(Buff_N, Buff_N, t) # Buff_N contains generation rate

    #!!! her we have allocs !!!! but alternatives lead to Forward diff erro, or seem to be 25% slower s#####
    GR  = Buff_N * c.ndim.G.light(t) - GR
    #@inbounds @simd  for i in eachindex(GR)
    #    GR[i]  = Buff_N[i] * c.ndim.G.light(t) - GR[i]
    #end


    # P BC
    du[1] = FP[1];      # no Ion flux at bc
    du[N+1] = -FP[N];   # no Ion flux at bc

    # ϕ BC
    du[N+2] = mE[1]-rₑ*mEₑ[end] -c.g.dx[1]*(1/2-P[1]/3-P[2]/6+δ*(n[1]/3+n[2]/6-χ*(p[1]/3+p[2]/6)))/λ² -rₑ*c.g.dxₑ[end]*(nₑ[end-1]/6+nₑ[end]/3-1/2)/λₑ²; # continuity
    du[2*N+2] = rₕ*mEₕ[1]-mE[end] -c.g.dx[end]*(1/2-P[end-1]/6-P[end]/3+δ*(n[end-1]/6+n[end]/3-χ*(p[end-1]/6+p[end]/3)))/λ² -rₕ*c.g.dxₕ[1]*(1/2-pₕ[1]/3-pₕ[2]/6)/λₕ²; #continuity

    # n BC
    du[2*N+3] = fn[1]-fnₑ[end]+(c.g.dx[1]*GR[1])/2#-c.ndim.Rₗ(nₑ[end],p[1]); # continuity
    du[3*N+3] = -fn[end]+c.g.dx[N]*GR[N]/2#-c.ndim.Rᵣ(n[N+1],pₕ[1]);
#    println(du[3*N+3])

    # p BC
    du[3*N+4] = fp[1]+c.g.dx[1]*GR[1]/2#-c.ndim.Rₗ(nₑ[end],p[1]);
    du[4*N+4] = fpₕ[1]-fp[end]+(c.g.dx[end]*GR[end])/2#-c.ndim.Rᵣ(n[N+1],pₕ[1]); # continuity
    @inbounds @simd for i in 1:N-1
        du[i+1] = FP[i+1] - FP[i];
        du[N+2+i] = mE[i+1]-mE[i]-cd[i]/λ²;
        du[2*N+3+i] = fn[i+1]-fn[i]+(c.g.dx[i+1]*GR[i+1]+c.g.dx[i]*GR[i])/2;
        du[3*N+4+i] = fp[i+1]-fp[i]+(c.g.dx[i+1]*GR[i+1]+c.g.dx[i]*GR[i])/2;
    end
    ### ETM ###
    du[4*N+5] = ϕₑ[1]#-c.ndim.V(t);
    du[4*N+Nₑ+5] = nₑ[1]-1;
    @inbounds @simd for i in 1:Nₑ-1
        du[4*N+5+i] = mEₑ[i+1]-mEₑ[i]-cdₑ[i]/λₑ²;
        du[4*N+Nₑ+5+i] = fnₑ[i+1]-fnₑ[i]
    end
    ### HTM ###
    du[4*N+2*Nₑ+Nₕ+4] = ϕₕ[end] + c.ndim.Vbi - c.ndim.V(t);
    du[4*N+2*Nₑ+2*Nₕ+4] = pₕ[end]-1;
    @inbounds @simd for i in 1:Nₕ-1
        du[4*N+2*Nₑ+4+i] = mEₕ[i+1]-mEₕ[i]-cdₕ[i]/λₕ²;
        du[4*N+2*Nₑ+Nₕ+4+i] = fpₕ[i+1]-fpₕ[i];
    end
    # Perform any additional step requested by the optional input argument c.mode
    if c.mode == :cc  # cosed circuit is default
    elseif c.mode == :oc #oopen circuit
        du[4*N+5] = (fnₑ[1]-( ϕₕ[end] +c.ndim.Vbi) *c.ndim.σₛₕ) ;  # no flux and shunt
        du[4*N+2*Nₑ+Nₕ+4] =  ϕₑ[1] # Potential reference at etl contact
        #du[N+1] = integrate(c.g.x,P)-1)^4;
    elseif c.mode == :precondition
        # Overwrite right-hand BC to ensure conservation of ion vacancies
        du[N+1] = integrate(c.g.x,P)-1;
        #du[4*N+2*Nₑ+Nₕ+4] = ϕₕ[end] +c.ndim.Vbi - c.ndim.V(t)
    else
        error("simulation mode $(c.mode) is not recongnized")
    end
    return nothing
end
