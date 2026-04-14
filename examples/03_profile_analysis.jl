using Pkg
Pkg.activate(".")
using AffineMapDiffEq
using StaticArrays
using DifferentialEquations
using KrylovKit
using Profile
# using ProfileView

# Use the In-place Vector version which was fastest
function Mathieu_IIP_Vector(du, u, h, p, t)
    ζ, δ, ϵ, b, τ, T = p
    F = 0.1 * (cos(2π * t / T)^10)
    du[1] = u[2]
    du[2] = -(δ + ϵ * cos(2π * t / T)) * u[1] - 2ζ * u[2] + b * h(p, t - τ)[1] + F
end

p = (0.02, 1.5, 0.15, 0.5, 2π, 2π)
τ = 2π
u0 = [1.0, 0.0]
h(p, t) = [0.0, 0.0]
prob = DDEProblem{true}(Mathieu_IIP_Vector, u0, h, (0.0, 2π), p; constant_lags=[τ])

Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-8)
Krylov_arg = (5, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=15, verbosity=0))

dp = dynamic_problemSampled(prob, Solver_args, 2π; 
                            Historyresolution=100, zerofixpont=false, 
                            affineinteration=1, Krylov_arg=Krylov_arg)

println("Warming up...")
affine(dp; p=p)

println("Profiling...")
Profile.clear()
@profile for i in 1:10
    affine(dp; p=p)
end

# Save profile data for analysis (ProfileView is interactive, so we'll use a text-based summary here)
Profile.print(format=:flat, sortedby=:count, mincount=100)
