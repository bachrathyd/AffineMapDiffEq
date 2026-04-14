# AffineMapDiffEq.jl

`AffineMapDiffEq.jl` is a Julia package for the stability analysis of Differential Equations (ODEs/DDEs) and non-smooth/nonlinear systems using **affine mapping** and **spectral analysis**.

## Key Features
- **Affine Mapping:** Efficiently find fixed points and Floquet multipliers of periodic solutions.
- **Spectral Analysis:** Stability analysis using Krylov-based methods (`KrylovKit.jl`) or Subspace Iteration (ISSI).
- **AD Compatible:** Seamlessly works with `ForwardDiff.jl` for sensitivity analysis and exact Jacobians.
- **Performance:** Optimized for type-stability and low allocations, supporting `StaticArrays`.

## Installation
Currently, this package is in development. You can use it by cloning the repository and activating the project.

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Quick Example: Mathieu Equation
The following example demonstrates the stability analysis of a delayed Mathieu equation with external forcing.

```julia
using AffineMapDiffEq
using StaticArrays
using DifferentialEquations
using KrylovKit

# 1. Define the Governing Equation
function DelayMathieu!(du, u, h, p, t)
    ζ, δ, ϵ, b, τ, T = p
    F = 0.1 * (cos(2π * t / T)^10)
    du[1] = u[2]
    du[2] = -(δ + ϵ * cos(2π * t / T)) * u[1] - 2ζ * u[2] + b * h(p, t - τ)[1] + F
end

# 2. Setup Parameters and Problem
p = (0.02, 1.5, 0.15, 0.5, 2π, 2π)
u0 = [1.0, 0.0]
h(p, t) = [0.0, 0.0]
prob = DDEProblem{true}(DelayMathieu!, u0, h, (0.0, 2π), p; constant_lags=[2π])

# 3. Solver and Mapping Configuration
Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :reltol => 1e-3)
Krylov_arg = (8, :LM, KrylovKit.Arnoldi(tol=1e-12))

dp = dynamic_problemSampled(prob, Solver_args, 2π; 
    Historyresolution=100,
    zerofixpont=false, 
    affineinteration=2,
    Krylov_arg=Krylov_arg)

# 4. Calculate Fixed Point and Floquet Multipliers
mus, saff, sol0 = affine(dp; p=p)
mu = mus[1] # Eigenvalues
println("Max Floquet multiplier: ", maximum(abs.(mu)))
```

## Folder Structure
- `src/`: Core implementation (types, mapping, spectrum, interpolation).
- `examples/`: Usage examples and demos.
- `docs/`: Documentation sources.
- `paper/`: Research publication sources.

## License
MIT
