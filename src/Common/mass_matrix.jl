function mass_matrix(g::Grid, p::NodimParameters;mode=:oc)

    for n in fieldnames(Grid)
        @eval $n = $g.$(n)
    end
    for n in fieldnames(NodimParameters)
        @eval $n = $p.$(n)
    end

    M11 = Tridiagonal(vcat(dx[1:end-1]/6,dx[end]/6), vcat(dx[1]/3,(dx[1:end-1]+dx[2:end])/3, dx[end]/3), vcat(dx[1]/6,dx[2:end]/6))
    M12 = spzeros(N+1,N+1)
    M15 = spzeros(N+1,Nₑ)
    M17 = spzeros(N+1,Nₕ)
    M33 = σ*copy(M11); M33[1,1] = σ*(dxₑ[end]/kₑ+dx[1])/3
    M36 = copy(M15); M36[1,end] = σ*dxₑ[end]/6
    M44 = σ*χ*copy(M11); M44[end,end]= σ*χ*(dx[end]+dxₕ[1]/kₕ)/3
    M48 = copy(M17); M48[end,1]= σ*χ*dxₕ[1]/6

    M51 = spzeros(Nₑ,N+1)
    M63 = copy(M51); M63[end,1] = σ*dxₑ[end]/kₑ/6
    M55 = spzeros(Nₑ,Nₑ)
    M57 = spzeros(Nₑ,Nₕ)
    M66 = σ*Tridiagonal(dxₑ[1:end-1]/6, vcat(0,(dxₑ[1:end-1]+dxₑ[2:end])/3),vcat(0,dxₑ[2:end-1])/6)

    M71 = spzeros(Nₕ,N+1)
    M75 = spzeros(Nₕ,Nₑ)
    M77 = spzeros(Nₕ,Nₕ)
    M84 = copy(M71); M84[1,end] = σ*χ*dxₕ[1]/kₕ/6
    M88 = σ*χ*Tridiagonal(vcat(dxₕ[2:end-1]/6,0), vcat((dxₕ[1:end-1]+dxₕ[2:end])/3, 0),dxₕ[2:end]/6)

#    M88[end,end-2:end] = zeros(3)
#    M66[1,1:3] =zeros(3)
    M= [M11 M12 M12 M12 M15 M15 M17 M17     # P equation
        M12 M12 M12 M12 M15 M15 M17 M17     # ϕ equation
        M12 M12 M33 M12 M15 M36 M17 M17     # n equation
        M12 M12 M12 M44 M15 M15 M17 M48     # p equation
        M51 M51 M51 M51 M55 M55 M57 M57     # ϕₑ equation
        M51 M51 M63 M51 M55 M66 M57 M57     # nₑ equation
        M71 M71 M71 M71 M75 M75 M77 M77     # ϕₕ equation
        M71 M71 M71 M84 M75 M75 M77 M88]    # pₕ equation


    if mode==:cc
        M
    elseif mode==:oc
        M[4*N+5, 4*N+Nₑ+5 : 4*N+Nₑ+6] = σ*dxₑ[1]*[1/3, 1/6]
        M
    elseif mode==:precondition
        M[1:N+1,:] = σ*M[1:N+1, :]
        M
    else
        error("Simulation mode $mode is not known, use :oc, :cc or :precondition")
    end


end
