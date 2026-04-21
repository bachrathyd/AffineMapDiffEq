using Pkg
Pkg.activate(".")
using AffineMapDiffEq
using StaticArrays
using DifferentialEquations
using LinearAlgebra
using KrylovKit
using InteractiveUtils

# 1. Setup a simple DDE problem for testing
function DelayMathieu!(du, u, h, p, t)
    ζ, δ, ϵ, b, τ, T = p
    du[1] = u[2]
    du[2] = -(δ + ϵ * cos(2π * t / T)) * u[1] - 2ζ * u[2] + b * h(p, t - τ)[1]
end

ζ = 0.02
δ_init = 1.5
ϵ_init = 0.15
τ = 2π
b = 0.5
T = 2π
p_init = (ζ, δ_init, ϵ_init, b, τ, T)

u0 = @SVector [1.0, 0.0]
h(p, t) = @SVector [0.0, 0.0]
probMapping = DDEProblem{true}(DelayMathieu!, u0, h, (0.0, T), p_init; constant_lags=[τ])
Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-6)
Krylov_arg = (2, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=10, verbosity=0))

dpMathieu = dynamic_problemSampled(probMapping, Solver_args, τ;
    Historyresolution=10, zerofixpont=false, affineinteration=1, Krylov_arg=Krylov_arg)

println("--- @code_warntype affine(dpMathieu; p=p_init) ---")
@code_warntype affine(dpMathieu; p=p_init)

println("\n--- @code_warntype LinMap(dpMathieu, [u0 for _ in 1:10]; p=p_init) ---")
s_test = [u0 for _ in 1:10]
@code_warntype LinMap(dpMathieu, s_test; p=p_init)
