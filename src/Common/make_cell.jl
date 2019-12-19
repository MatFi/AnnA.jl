mutable struct Cell{V,M,J,G,O,P,NP,C,S,A<:AlgControl}
    M::M        # Mass Matrix
    Jac::J      # Jacobian pattern
    g::G        # Grid
    o::O        # Operators
    parameters::P   # User choosen input parameters
    ndim::NP    # Non dimensionalized parameters
    rhs::C      # Rhs! function
    mode::S     # :cc  closes circuit, :oc open circuit, :precondition
    alg_ctl::A  # Alg control type
    u0::V       # Initial_conditions
end


function Cell(
    p::AbstractParameters;
    op_flavor::Symbol=:sparse,
    mode::Symbol=:cc,
    alg_ctl::AlgControl=AlgControl()
)
    if mode == :oc && p.V(0)>0
        @warn "nonzero potential at t=0 in :oc mode
            will be ignored during the initialisation"
    end
    grid   = Grid(p)
    ndim= NonDimensionalize(p)
    operators  = Operators(grid;flavor=op_flavor)
    rhs = Rhs(grid)
    massMatrix   = mass_matrix(grid, ndim;mode=mode)
    Jac = get_jac_sparse_pattern(grid;mode=mode)
    u0 = init_guess(grid,ndim)
    Cell(massMatrix,Jac,grid,operators,p,ndim,rhs,mode,alg_ctl,u0)
end
