function mass_matrix(g::Grid, p::NodimParameters;mode=:oc)

    #=for n in fieldnames(Grid)
        @eval $n = $g.$(n)
    end

    for n in fieldnames(NodimParameters)
        @eval $n = $p.$(n)
    end
    =#
    p = unpac_struct(g,p)

    M11 = Tridiagonal(vcat(p[:dx][1:end-1]/6,p[:dx][end]/6), vcat(p[:dx][1]/3,(p[:dx][1:end-1]+p[:dx][2:end])/3, p[:dx][end]/3), vcat(p[:dx][1]/6,p[:dx][2:end]/6))
    M12 = spzeros(p[:N]+1,p[:N]+1)
    M15 = spzeros(p[:N]+1,p[:Nₑ])
    M17 = spzeros(p[:N]+1,p[:Nₕ])
    M33 = p[:σ]*copy(M11); M33[1,1] = p[:σ]*(p[:dxₑ][end]/p[:kₑ]+p[:dx][1])/3
    M36 = copy(M15); M36[1,end] = p[:σ]*p[:dxₑ][end]/6
    M44 = p[:σ]*p[:χ]*copy(M11); M44[end,end]= p[:σ]*p[:χ]*(p[:dx][end]+p[:dxₕ][1]/p[:kₕ])/3
    M48 = copy(M17); M48[end,1]= p[:σ]*p[:χ]*p[:dxₕ][1]/6

    M51 = spzeros(p[:Nₑ],p[:N]+1)
    M63 = copy(M51); M63[end,1] = p[:σ]*p[:dxₑ][end]/p[:kₑ]/6
    M55 = spzeros(p[:Nₑ],p[:Nₑ])
    M57 = spzeros(p[:Nₑ],p[:Nₕ])
    M66 = p[:σ]*Tridiagonal(p[:dxₑ][1:end-1]/6, vcat(0,(p[:dxₑ][1:end-1]+p[:dxₑ][2:end])/3),vcat(0,p[:dxₑ][2:end-1])/6)

    M71 = spzeros(p[:Nₕ],p[:N]+1)
    M75 = spzeros(p[:Nₕ],p[:Nₑ])
    M77 = spzeros(p[:Nₕ],p[:Nₕ])
    M84 = copy(M71); M84[1,end] = p[:σ]*p[:χ]*p[:dxₕ][1]/p[:kₕ]/6
    M88 = p[:σ]*p[:χ]*Tridiagonal(vcat(p[:dxₕ][2:end-1]/6,0), vcat((p[:dxₕ][1:end-1]+p[:dxₕ][2:end])/3, 0),p[:dxₕ][2:end]/6)

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
        
    elseif mode==:oc
        M[4*p[:N]+5, 4*p[:N]+p[:Nₑ]+5 : 4*p[:N]+p[:Nₑ]+6] = p[:σ]*p[:dxₑ][1]*[1/3, 1/6]
    #    M[1:p[:N]+1,:] = p[:σ]*M[1:p[:N]+1, :]
        
    elseif mode==:precondition
        M[1:p[:N]+1,:] = p[:σ]*M[1:p[:N]+1, :]
        
    else
        error("Simulation mode $mode is not known, use :oc, :cc or :precondition")
    end
    dropzeros(sparse(M))

end
