# Example 03: Benchmarking Mapping Performance with Different State Types
# Investigating StaticArrays vs. Mutating Vectors vs. Dual Numbers

using Pkg
Pkg.activate(".")
using AffineMapDiffEq
using StaticArrays
using DifferentialEquations
using LinearAlgebra
using KrylovKit
using BenchmarkTools

# 1. Standard Vector (Out-of-place)
function Mathieu_OOP_Vector(u, h, p, t)
    ζ, δ, ϵ, b, τ, T = p
    F = 0.1 * (cos(2π * t / T)^10)
    dx = u[2]
    ddx = -(δ + ϵ * cos(2π * t / T)) * u[1] - 2ζ * u[2] + b * h(p, t - τ)[1] + F
    [dx, ddx]
end

# 2. Standard Vector (In-place)
function Mathieu_IIP_Vector(du, u, h, p, t)
    ζ, δ, ϵ, b, τ, T = p
    F = 0.1 * (cos(2π * t / T)^10)
    du[1] = u[2]
    du[2] = -(δ + ϵ * cos(2π * t / T)) * u[1] - 2ζ * u[2] + b * h(p, t - τ)[1] + F
end

# 3. Mutable MVector (Out-of-place) - MVector is better for KrylovKit internal mutation
function Mathieu_OOP_MVector(u, h, p, t)
    ζ, δ, ϵ, b, τ, T = p
    F = 0.1 * (cos(2π * t / T)^10)
    dx = u[2]
    ddx = -(δ + ϵ * cos(2π * t / T)) * u[1] - 2ζ * u[2] + b * h(p, t - τ)[1] + F
    @MArray [dx, ddx]
end

# Parameters
p = (0.02, 1.5, 0.15, 0.5, 2π, 2π)
τ = 2π
Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-8)
Krylov_arg = (5, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=15, verbosity=0))

function run_benchmark(name, func, u0, h, iip)
    println("\n=== Testing: $name ===")
    prob = DDEProblem{iip}(func, u0, h, (0.0, 2π), p; constant_lags=[τ])
    
    dp = dynamic_problemSampled(prob, Solver_args, 2π; 
                                Historyresolution=100, zerofixpont=false, 
                                affineinteration=1, Krylov_arg=Krylov_arg)
    
    # Warmup and verification
    try
        @time mu, saff, sol0 = affine(dp; p=p)
        println("Success! Max multiplier: ", maximum(abs.(mu[1])))
        
        # Benchmark
        println("Benchmarking $name...")
        b = @benchmark affine($dp; p=$p) samples=3 seconds=10
        display(b)
    catch e
        @error "Failed $name" exception=(e, catch_backtrace())
    end
end

println("Starting benchmarks...")

run_benchmark("Vector (Out-of-place)", Mathieu_OOP_Vector, [1.0, 0.0], (p, t) -> [0.0, 0.0], false)
run_benchmark("Vector (In-place)", Mathieu_IIP_Vector, [1.0, 0.0], (p, t) -> [0.0, 0.0], true)
run_benchmark("MVector (Out-of-place)", Mathieu_OOP_MVector, @MVector([1.0, 0.0]), (p, t) -> @MVector([0.0, 0.0]), false)
