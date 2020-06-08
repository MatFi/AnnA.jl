using AnnABase,  Unitful

parm = AnnABase.Parameters(
    ε = 50,
    T = 300 * u"K",
    b = 500 * u"nm",
    bₑ = 50 * u"nm",
    bₕ = 50 * u"nm",
    dₑ = 2.4e18 * u"cm^-3",
    dₕ = 1.9e18 * u"cm^-3",
    N₀ = 1.4e22 * u"m^-3",
    Dᵢ₀ = 8e-3 * u"m^2/s",
    τₚ = 4e-6 * u"s", # AnnAOptim.OptimVariable(3e-5u"nm",3e-7u"nm",:τₙ;logspaced=true)
    τₙ = 0.35e-6 * u"s", # for good iv do 0.1e-6
    Ecₑ = -4.06 * u"eV", #[-4.2,-3.7]u"m^-3"
    Evₕ = (-5.00 + 0.005) * u"eV", #[-5.3,-5]u"m^-3"
    light = (t) -> 0,
    vₚₑ = 6e-1 * u"m/s",
    vₙₑ = 1.0e-4 * u"m/s",
    k₂ₑ = 0 * u"m^4/s",
    vₚₕ = 0.5e-1 * u"m/s",
    vₙₕ = 1e-0 * u"m/s",
    Rₛₕ = Inf * u"V/A*mm^2",
    gc = 6.98e18 * u"cm^-3", #10.1021/acs.jpcc.6b10914
    gv = 2.49e18 * u"cm^-3", #10.1021/acs.jpcc.6b10914
    gcₑ = 2.0e25 * u"m^-3",
    gvₕ = 2.0e25 * u"m^-3",
    Dₑ = 1e-8 * u"m^2/s",
    Dₕ = 1e-8 * u"m^2/s",
   # Dₙ = 9e-6 * u"m^2/s",
   # Dₚ = 1.0e-6 * u"m^2/s",
    Dₙ = 5.568e-6 * u"m^2/s", #10.1021/acsenergylett.9b02310
    Dₚ = 6.2641e-5* u"m^2/s", #10.1021/acsenergylett.9b02310
    k₂ = 4.78e-11 * u"cm^3/s", # 10.1103/PhysRevApplied.6.044017 w photonrecyceling
    Ec = -3.7 * u"eV",       # Perovskite Conduction band energy
    Ev = -5.3 * u"eV",        # Perovskite Valence band energy
    freeze_ions = false,
    N = 500,
);

prob_ocvd = AnnABase.OCVDProblem.([parm],
    50 * u"s",
    1e6 * u"s";
    alg_control = AnnABase.AlgControl(ss_tol = 1e-3, reltol = 1e-5, abstol = 1e-5),
);
sol_ocvd = AnnABase.solve.(prob_ocvd);

prob_iv =
    AnnABase.IVProblem(
        parm,
        (-0.2 * u"V", 1.2 * u"V"),
        0.01 * u"mV/s";
        double_sweep = true,
        alg_control = AnnABase.AlgControl(
            ss_tol = 1e-4,
            reltol = 1e-8,
            abstol = 1e-8,
        ),
    )

;
sol_iv = AnnABase.solve(prob_iv);


using LinearAlgebra, AnnABase.DiffEqBase,OrdinaryDiffEq, DiffEqBase
c=AnnABase.Cell(parm)
M=c.M
u0=c.u0
fu=c.rhs
tspan = (0.0,1.0)


N = 5
M = collect(I((5)))
M[1:2, 1:2] .= 0
M[end, 1] = 1
u0 = ones(N)
M = [0 0 0 0 0
     0 0 0 0 0
     0 0 1 0 0
     0 0 0 1 0
     1 0 0 0 1]

function f(du, u, p, t)
    du .= 0 * u
    du[1] = 1- u[1]
    du[2] = 1 - u[2]

end
odefun = ODEFunction(f;mass_matrix = M,)
prob = ODEProblem(odefun, u0, tspan)
sol = solve(prob, Rodas5())

OrdinaryDiffEq. ArrayInterface.issingular(M)
M != I

#=
ii
=
