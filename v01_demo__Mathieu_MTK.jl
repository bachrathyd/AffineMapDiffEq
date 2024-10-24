
#TODO: Working simulation, 
#TODO: But non-working Affine mapping, in DDE_mapping, because the problems with the types!!!

using ModelingToolkit, DifferentialEquations
#using SymbolicIndexingInterface: is_markovian
using ModelingToolkit: t_nounits as t, D_nounits as D


using StaticArrays
using Plots
@parameters ζ = 0.02  δ = 1.5 [tunable=true] ϵ = 0.15 τ = 2pi b = 0.5 T = 2pi
@variables x(..)
tau = 2pi
F = 0.1 * (cos(2pi * t / T) .^ 10)
eqs = [D(D(x(t))) + 2 * ζ *D(x(t))+(δ + ϵ * cos(2pi * t / T)) * x(t) ~ b* x(t - tau)+F]
      
@mtkbuild sys = System(eqs, t)
ModelingToolkit.is_dde(sys)
#@test !is_markovian(sys)

tspan = (0.0,2pi*100.0)
probM = DDEProblem(sys,
    SVector(x(t)=> 1.0, D(x(t)) => 0.0),
    tspan,
    constant_lags = [tau], jac = true; u0_constructor = x -> SVector(x...))


Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-4)#
@time sol_mtk = solve(probM; Solver_args...);
#@benchmark sol_mtk = solve(prob, alg, reltol = 1e-4)

plot(sol_mtk,xlim=tspan,linestyle=:dash, linewidth=3)





# ---------------- Affine mapping ---------------------------
using DDE_mapping
using KrylovKit
Neig=8#number of required eigen values
Krylov_arg=(Neig,:LM, KrylovKit.Arnoldi(tol=1e-12,krylovdim=8+5,verbosity=0));

τmax=tau #maximal timedelay in the mapping
Nstep = 100 # discretization number of the mapping
Timeperiod=2pi # timeperiod of the mapping

#Creating the problem
dpMathieu = dynamic_problemSampled(probM, Solver_args, τmax,
Timeperiod; Historyresolution=Nstep,
    zerofixpont=false,    affineinteration=1,
    Krylov_arg=Krylov_arg);

#solvig the problem (mu: Floquet multiplier, saff: discretized foxed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
@time mu, saff, sol0 = affine(dpMathieu;p=(0.02, 1.5, 0.15, 0.5, 6.283185307179586, 6.283185307179586))
