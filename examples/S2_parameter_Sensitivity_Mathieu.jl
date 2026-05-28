# Example 02: Delayed Mathieu Equation - Stability and Amplitude Analysis
# Structured Visualization: Time Series, Periodic Orbit, and Spectrum

using Pkg;
Pkg.activate(".");
using AffineMapDiffEq
using StaticArrays
using DifferentialEquations
using LinearAlgebra
using KrylovKit
using GLMakie
using MDBM

# 1. Define the Governing Equation (In-place)
function DelayMathieu!(du, u, h, p, t)
    ζ, δ, ϵ, b, τ, T = p
    F = 0.1 * (cos(2π * t / T)^10)
    du[1] = u[2]
    du[2] = -(δ + ϵ * cos(2π * t / T)) * u[1] - 2ζ * u[2] + b * h(p, t - τ)[1] + F
end

# 2. Setup Parameters
const ζ =0.01
const δ_init = 1.5
const ϵ_init = 0.15
const τ = 2π
const b = 0.05
const T = 2π
const p_init = (ζ, δ_init, ϵ_init, b, τ, T)

# 3. Create Figure Layout
fig = Figure(size=(1200, 900))
ax_time = GLMakie.Axis(fig[1, 1:2], title="Long Simulation Time Series", xlabel="Time (t)", ylabel="x(t)")
ax_orbit = GLMakie.Axis(fig[2, 1], title="Periodic Orbit Comparison", xlabel="x", ylabel="dx/dt")
ax_complex = GLMakie.Axis(fig[2, 2], title="Floquet Multipliers (Complex Plane)", xlabel="Re(μ)", ylabel="Im(μ)", aspect=GLMakie.DataAspect())

# 4. Long Simulation
u0 = @MArray [1.0, 0.0]
h(p, t) = @MArray [0.0, 0.0]
prob_long = DDEProblem{true}(DelayMathieu!, u0, h, (0.0, T * 200.0), p_init; constant_lags=[τ])
Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-6)

println("Running long simulation...")
@time sol_long = solve(prob_long; Solver_args...)

# Plot Time Series (Top Row)
lines!(ax_time, sol_long.t, [u[1] for u in sol_long.u], color=:blue, linewidth=1)

# Extract steady-state periodic orbit (last period)
t_end = sol_long.t[end]
t_orbit = range(t_end - T, t_end, length=200)
sol_steady = [sol_long(t) for t in t_orbit]
lines!(ax_orbit, getindex.(sol_steady, 1), getindex.(sol_steady, 2), color=:blue, linewidth=3, label="Long Sim (Steady State)")

display(fig) # Show intermediate progress

# 5. Affine Mapping Analysis
println("Calculating Affine results...")
Neig = 6#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-15, krylovdim=Neig * 2 + 2, verbosity=0, maxiter=10))

probMapping = DDEProblem{true}(DelayMathieu!, u0, h, (0.0, T), p_init; constant_lags=[τ])

dpMathieu = dynamic_problemSampled(probMapping, Solver_args, τ;
    Historyresolution=100, zerofixpont=false, affineinteration=2, Krylov_arg=Krylov_arg)

#Warm up"
@time mus, saff, sol0 = affine(dpMathieu; p=p_init);
@time mus, saff, sol0 = affine(dpMathieu; p=p_init);
mu_vals = mus[1]

#@code_warntype  affine(dpMathieu; p=p_init);

# Plot Affine Fixed Point (Bottom Left)
lines!(ax_orbit, getindex.(saff, 1), getindex.(saff, 2), color=:red, linestyle=:dash, linewidth=2, label="Affine Fixed Point")
scatter!(ax_orbit, getindex.(saff, 1), getindex.(saff, 2), color=:red, markersize=6)

lines!(ax_orbit, [u[1] for u in sol0.u], [u[2] for u in sol0.u], color=:magenta, linewidth=1, label="Affine Periodic Orbit")
axislegend(ax_orbit)

# Plot Complex Plane (Bottom Right)
# Unit Circle
θ = range(0, 2π, length=100)
lines!(ax_complex, cos.(θ), sin.(θ), color=:black, linestyle=:dash)
# Multipliers
scatter!(ax_complex, real.(mu_vals), imag.(mu_vals), color=:red, markersize=10, label="μ")
vlines!(ax_complex, [0], color=:gray, linewidth=0.5)
hlines!(ax_complex, [0], color=:gray, linewidth=0.5)

display(fig)
save("examples/02_mathieu_analysis.png", fig)



## ─── Robustness Margin C ─────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "src", "robust_control.jl"))

# Predefined ranges for parameter variations (as a retio of the nominal values)
const ζ_var = 0.02
const δ__var = 0.02
const ϵ_var = 0.02
const b_var = 0.4
const τ_var = 0.1 
const T_var = 0.0
# Uncertainty tuple aligned with p_init = (ζ, δ, ϵ, b, τ, T)
const delta_p = (ζ_var, δ__var, ϵ_var, b_var, τ_var, T_var)

println("\nComputing robustness margin C at nominal parameters...")
@time C, C_vec, mu_C, sens = robustness_margin(dpMathieu, p_init, delta_p)

println("Floquet multipliers |μ|     = ", round.(abs.(mu_C); digits=6))
println("Per-root margin     C_j     = ", round.(C_vec; digits=4))
println("Worst-case margin   C       = ", round(C; digits=4))
println("\nSensitivity matrix  ∂|μ_j|/∂p_i:")
println("  (rows = multiplier index, cols = ζ  δ  ϵ  b  τ  T)")
display(round.(sens; digits=4))


## ─── Robustness Map over (δ, ϵ) ──────────────────────────────────────────────
# Spectral radius (nominal stability) + worst-case robustness margin C
# computed on a grid in (δ, ϵ), using a Threaded for-loop like in
# 02_delayed_mathieu_amplitude.jl. The uncertainty ratios above are interpreted
# as |Δp_i| = ratio · |p_i_nominal|, so the uncertainty box scales with the
# local nominal δ, ϵ values at each grid point.

println("\nStarting Robustness Map (threaded over (δ, ϵ))...")

const δv = LinRange(-1.0, 15.0, 80)
const ϵv = LinRange(-0.001, 5.0, 41)
const b_chart = 0.05

# Fast dynamic problem for the sweep. perturbation_size defaults to 0.0, so the
# robustness margin uses the exact ForwardDiff (Dual) sensitivity branch. Neig is
# both the number of tracked multipliers and the Schur-subspace size used to
# compress the operator for the Dual sensitivity; a larger Neig captures the
# dominant mode's left eigenvector more fully → more accurate C. ~10 balances
# accuracy and speed.
Neig_fast = 10
Krylov_arg_fast = (Neig_fast, :LM,
    KrylovKit.Arnoldi(tol=1e-12, krylovdim=Neig_fast + 20, verbosity=0, maxiter=10))
dp_fast = dynamic_problemSampled(probMapping, Solver_args, τ;
    Historyresolution=25, zerofixpont=false, affineinteration=1, Krylov_arg=Krylov_arg_fast)

Spec_map = zeros(length(ϵv), length(δv))
C_map    = fill(NaN, length(ϵv), length(δv))

@time Threads.@threads for j in 1:length(δv)
    println("  Processing δ = ", round(δv[j]; digits=4), "  (", j, "/", length(δv), ")")
    δ_v = δv[j]
    for i in 1:length(ϵv)
        ϵ_v = ϵv[i]
        p_loc = (ζ, δ_v, ϵ_v, b_chart, τ, T)

        # Uncertainty box scaled by the local nominal parameter values
        delta_p_loc = (ζ_var  * abs(ζ),
                       δ__var * abs(δ_v),
                       ϵ_var  * abs(ϵ_v),
                       b_var  * abs(b_chart),
                       τ_var  * abs(τ),
                       T_var  * abs(T))

        # robustness_margin early-exits with C=-Inf if |μ| ≥ 1, so it does
        # not pay the FD sensitivity cost in the unstable region.
        C_loc, _, mu_loc, _ = redirect_stdout(devnull) do
            robustness_margin(dp_fast, p_loc, delta_p_loc)
        end

        ρ = maximum(abs.(mu_loc))
        Spec_map[i, j] = ρ
        if isfinite(C_loc) && ρ < 1.0
            C_map[i, j] = C_loc
        end
    end
end

println("  Robustness map: ",
        count(isfinite, C_map), " stable / ",
        length(C_map), " total grid points")

## ── Plot: spectral radius everywhere + robustness margin C on stable region ──
fig_C = Figure(size=(1400, 650))
ax_rho    = GLMakie.Axis(fig_C[1, 1],
    title="Spectral radius ρ = max|μ|  (stability boundary ρ=1)", xlabel="δ", ylabel="ϵ")
ax_stable = GLMakie.Axis(fig_C[1, 2],
    title="Robustness margin C (stable region only)", xlabel="δ", ylabel="ϵ")

C_plot = copy(C_map)
C_plot[.!isfinite.(C_plot)] .= NaN
finite_vals = filter(isfinite, C_plot)
crange = isempty(finite_vals) ? (0.0, 1.0) : (minimum(finite_vals), maximum(finite_vals))

# Left: spectral radius on the stable region only (unstable → NaN)
Spec_stable = copy(Spec_map)
Spec_stable[Spec_stable .>= 1.0] .= NaN
hm_rho = heatmap!(ax_rho, δv, ϵv, log.(Spec_stable'), colormap=:inferno)
Colorbar(fig_C[2, 1], hm_rho, vertical=false, label="log(ρ)  (unstable → NaN)")

# Right: stable region only (unstable → NaN)
hm_stable = heatmap!(ax_stable, δv, ϵv, C_plot',
    colormap=:viridis, colorrange=crange)
Colorbar(fig_C[2, 2], hm_stable, vertical=false, label="C")

save("examples/parameter_sensitivity_robustness_map.png", fig_C)
display(fig_C)

println("Robustness map complete. Figure saved.")

