
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D

using StaticArrays
using Plots
@parameters ζ = 0.02  δ = 1.5 [tunable=true] ϵ = 0.15 τ = 2pi b = 0.5 T = 2pi
@variables  x(..)=1.0 y(t)=0.0 
tau = 2pi
T_per=2pi
F = 0.1 * (cos(2pi * t / T) .^ 10)
eqs = [D(x(t)) ~ y,
    D(y) + 2 * ζ *D(x(t))+(δ + ϵ * cos(2pi * t / T)) * x(t) ~ b* x(t - tau)+F]
      
@mtkbuild sys = System(eqs, t)
ModelingToolkit.is_dde(sys)

tspan = (0.0,T_per*1000.0)
prob_MTK = DDEProblem{false}(sys,
    [],
    tspan,
    constant_lags = [tau],jac=true,u0_constructor = SVector{2})


Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-5)#
@time sol_mtk = solve(prob_MTK; Solver_args...);

plot(sol_mtk,xlim=tspan,linestyle=:dash, linewidth=3)


# ---------------- Affine mapping ---------------------------

Base.:+(a::SVector, b::Bool) = a .+ b
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?

using DDE_mapping
using KrylovKit
Neig=8#number of required eigen values
Krylov_arg=(Neig,:LM, KrylovKit.Arnoldi(tol=1e-12,krylovdim=8+5,verbosity=0));

τmax=tau #maximal timedelay in the mapping
Nstep = 1000 # discretization number of the mapping
Timeperiod=T_per # timeperiod of the mapping

#Creating the problem
dpMathieu_MTK = dynamic_problemSampled(prob_MTK, Solver_args, τmax,
Timeperiod; Historyresolution=Nstep,
    zerofixpont=false,    affineinteration=1,
    Krylov_arg=Krylov_arg);

#solvig the problem (mu: Floquet multiplier, saff: discretized foxed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
@time mu, saff, sol0 = affine(dpMathieu;p=(0.02, 1.5, 0.15, 0.5, 6.283185307179586, 6.283185307179586));


# Comparing the solutions:

t_select_period=(-T_per:0.01:0) .+ sol_mtk.t[end]

plot(getindex.(sol(t_select_period).u,1), getindex.(sol(t_select_period).u,2), marker=:circle,markersize=6,lab="")
plot!(getindex.(saff,1), getindex.(saff,2), marker=:circle,markersize=4,lab="")#
plot!(sol0[1, :], sol0[2, :], marker=:circle,lw = 0,lab="")#marker=:cross,markersize=2)#
##
#in log scale
plot(log.(abs.(mu[1])))
#in complex plane
scatter(mu[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))