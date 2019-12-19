function get_jac_sparse_pattern(g::Grid;mode::Symbol=:oc)
    N=g.N
    Nₑ=g.Nₑ
    Nₕ=g.Nₕ
    J11 = Tridiagonal(N+1,1,1,1)
    J13 = spzeros(N+1,N+1)
    J15 = spzeros(N+1,Nₑ)
    J17 = spzeros(N+1,Nₕ)
    J25 = copy(J15); J25[1,Nₑ] = 1
    J27 = copy(J17); J27[N+1,1] = 1
    J51 = spzeros(Nₑ,N+1)
    J52 = copy(J51); J52[Nₑ,1] = 1
    J55 = Tridiagonal(Nₑ,1,1,1); J55[1,2] = 0
    J56 = copy(J55); J56[1,1] = 0
    J57 = spzeros(Nₑ,Nₕ)
    J71 = spzeros(Nₕ,N+1)
    J72 = copy(J71); J72[1,N+1]=1
    J75 = spzeros(Nₕ,Nₑ)
    J77 = Tridiagonal(Nₕ,1,1,1); J77[Nₕ,Nₕ-1] = 0
    J78 = copy(J77); J78[Nₕ,Nₕ] = 0

    J   = [ J11 J11 J13 J13 J15 J15 J17 J17
            J11 J11 J11 J11 J25 J25 J27 J27
            J13 J11 J11 J11 J25 J25 J17 J17
            J13 J11 J11 J11 J15 J15 J27 J27
            J51 J52 J52 J51 J55 J56 J57 J57
            J51 J52 J52 J51 J56 J55 J57 J57
            J71 J72 J71 J72 J75 J75 J77 J78
            J71 J72 J71 J72 J75 J75 J78 J77 ]

    if mode == :oc
        J[4*N+5,:] = [zeros(1,4*N+4) [1] [1] zeros(1,Nₑ-2) [1] [1] zeros(1,Nₑ-2) zeros(1,2*Nₕ)]
        #Shunt resitance dependant current
        J[4*N+5,4*N+4+2*Nₑ+Nₕ]=1
        J[4*N+2*Nₑ+Nₕ+4,:] =  [zeros(1,4*N+4) [1] zeros(1,Nₑ-1) zeros(1,Nₑ) zeros(1,Nₕ-1) [0] zeros(1,Nₕ)]
    elseif mode == :precondition
        # enforces the conservation on ions
        J[N+1,:] = [ones(1,N+1) zeros(1,3*N+2*Nₑ+2*Nₕ+3)]
    elseif mode == :cc
    else error("simulation mode $mode is not recognised")
    end

    return sparse(J)
end
