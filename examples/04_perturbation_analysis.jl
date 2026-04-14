using AffineMapDiffEq
using DifferentialEquations
using StaticArrays
using LinearAlgebra
using KrylovKit
using BenchmarkTools
using GLMakie

# --- Problem Setup (Mathieu Equation) ---
# System parameters
δ = 1.0
ε = 0.2
a1 = 0.1
ω = 1.0
τ = 2π

# DDE definition: x'' + a1*x' + (δ + ε*cos(ω*t))*x = 0
function mathieu_dde(du, u, h, p, t)
    x = u[1]
    dx = u[2]
    # In this simple ODE-like DDE, we don't even use h, 
    # but we define it as a DDE for testing the mapping.
    du[1] = dx
    du[2] = -(δ + ε * cos(ω * t)) * x - a1 * dx
end

u0 = [1.0, 0.0]
tspan = (0.0, τ)
h(p, t) = [1.0, 0.0] # Constant history
probMapping = DDEProblem(mathieu_dde, u0, h, tspan)

# Solver and Mapping Arguments
Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-8)
τmax = τ
Nstep = 20
Neig = 2
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=Neig + 4, verbosity=0));

# --- Fast Test ---
println("=== Fast Test ===")
dp_exact = dynamic_problemSampled(probMapping, Solver_args, τmax; 
    Historyresolution=Nstep, perturbation_size=0.0, Krylov_arg=Krylov_arg)

dp_fd = dynamic_problemSampled(probMapping, Solver_args, τmax; 
    Historyresolution=Nstep, perturbation_size=1e-6, Krylov_arg=Krylov_arg)

println("Running affine (exact)...")
mu_exact, _, _ = affine(dp_exact)
println("Running affine (FD 1e-6)...")
mu_fd, _, _ = affine(dp_fd)

println("Exact multipliers (first 2): ", mu_exact[1][1:2])
println("FD 1e-6 multipliers (first 2): ", mu_fd[1][1:2])
println("Diff: ", norm(mu_exact[1][1:2] - mu_fd[1][1:2]))

# --- Benchmarking ---
println("\n=== Benchmarking Mapping Turn (TheMapping) ===")
# We benchmark the mapping call inside schursolve implicitly by timing the whole affine call 
# but for a more precise look we could benchmark LinMap.
println("Timing affine (exact):")
@btime affine($dp_exact) samples=5

println("Timing affine (FD):")
@btime affine($dp_fd) samples=5

# --- Error Analysis ---
println("\n=== Error Analysis vs Perturbation Size ===")
epsilons = 10.0 .^ (-1:-1:-16)
errors = Float64[]

# Get the "true" multipliers (exact AD)
true_mu = mu_exact[1]

for eps in epsilons
    dp_loc = dynamic_problemSampled(probMapping, Solver_args, τmax; 
        Historyresolution=Nstep, perturbation_size=eps, Krylov_arg=Krylov_arg)
    
    mu_loc, _, _ = affine(dp_loc)
    
    # Handle potentially different number of eigenvalues
    n_comp = min(length(true_mu), length(mu_loc[1]))
    err = norm(true_mu[1:n_comp] - mu_loc[1][1:n_comp])
    push!(errors, err)
    println("eps = 1e-$(Int(abs(log10(eps)))), error = $err, n_mu = $(length(mu_loc[1]))")
end

# --- Plotting ---
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1], 
    xlabel = "Perturbation Size (epsilon)", 
    ylabel = "Error (norm)",
    xscale = log10, 
    yscale = log10,
    title = "Finite Difference Error vs Perturbation Size")

scatterlines!(ax, epsilons, errors, color = :blue, markersize = 10)
# Add a line for machine epsilon reference
vlines!(ax, [sqrt(eps(Float64))], color = :red, linestyle = :dash, label = "sqrt(eps)")

save("examples/perturbation_error_analysis.png", fig)
println("\nPlot saved to examples/perturbation_error_analysis.png")

display(fig)
