# AffineMapDiffEq.jl

`AffineMapDiffEq.jl` is a Julia package for the stability analysis of Differential Equations (ODEs/DDEs) and non-smooth/nonlinear systems using **affine mapping** and **spectral analysis**.

## Key Features
- **Affine Mapping:** Efficiently find fixed points and Floquet multipliers of periodic solutions.
- **Spectral Analysis:** Stability analysis using Krylov-based methods (`KrylovKit.jl`).
- **AD Compatible:** Seamlessly works with `ForwardDiff.jl` for sensitivity analysis.
- **DDE Support:** Built on top of `DifferentialEquations.jl` for robust delay equation solving.

## Installation
Currently, this package is in development. You can use it by cloning the repository and activating the project.

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Quick Example: Mathieu Equation
```julia
using AffineMapDiffEq
using StaticArrays
using DifferentialEquations
using KrylovKit

# Define the DDE
function DelayMathieu(u, h, p, t)
    ζ, δ, ϵ, b, τ, T = p
    dx = u[2]
    ddx = -(δ + ϵ * cos(2π * t / T)) * u[1] - 2ζ * u[2] + b * h(p, t - τ)[1] + 0.1*(cos(2π*t/T)^10)
    @MArray [dx, ddx]
end

# Parameters
p = (0.02, 1.5, 0.15, 0.5, 2π, 2π)
u0 = @MArray [1.0, 0.0]
h(p, t) = @MArray [0.0, 0.0]
prob = DDEProblem(DelayMathieu, u0, h, (0.0, 2π), p; constant_lags=[2π])

# Setup Sampling Problem
dp = dynamic_problemSampled(prob, Dict(:alg => MethodOfSteps(BS3())), 2π, 2π)

# Calculate stability
mu, saff, sol0 = affine(dp; p=p)
println("Max Floquet multiplier: ", maximum(abs.(mu[1])))
```

## Folder Structure
- `src/`: Core implementation.
- `examples/`: Usage examples and demos.
- `docs/`: Documentation.
- `paper/`: Research publication sources.

## License
MIT
