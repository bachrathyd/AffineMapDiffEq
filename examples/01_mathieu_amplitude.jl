# Example 01: Delayed Mathieu Equation - Stability and Amplitude Analysis
# Structured Visualization: Time Series, Periodic Orbit, and Spectrum

#using Pkg;
#Pkg.activate(".");
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
ζ = 0.02
δ_init = 1.5
ϵ_init = 0.15
τ = 2π
b = 0.5
T = 2π
p_init = (ζ, δ_init, ϵ_init, b, τ, T)

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

@time mus, saff, sol0 = affine(dpMathieu; p=p_init);
mu_vals = mus[1]

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
#save("examples/01_mathieu_analysis.png", fig)

## 6. Stability Chart (Brute Force + MDBM Overlay)
println("Starting Stability Chart analysis...")
δv = LinRange(-1.0, 5.0, 60)
ϵv = LinRange(-0.001, 10.0, 61)
b_chart = 0.05
Spec_chart = zeros(length(ϵv), length(δv))
Amp_chart = zeros(length(ϵv), length(δv))


Neig = 2#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=Neig + 10, verbosity=0, maxiter=10));
dp_fast = dynamic_problemSampled(probMapping, Solver_args, τ;
    Historyresolution=25, zerofixpont=false, affineinteration=1, Krylov_arg=Krylov_arg)

@time Threads.@threads for j in 1:length(δv)
    println("Calculating for δ = $(δv[j])...")# percentage: $(round(j/length(δv)*100, digits=2))%")
    δ_v = δv[j]
    for i in 1:length(ϵv)
        ϵ_v = ϵv[i]
        mu_loc, s0_loc, _ = affine(dp_fast; p=(ζ, δ_v, ϵ_v, b_chart, τ, T))
        Spec_chart[i, j] = maximum(abs.(mu_loc[1]))
        Amp_chart[i, j] = norm(getindex.(s0_loc, 1))
    end

end

println("MDBM calculation...")

δv = LinRange(-1.0, 5.0, 10)
ϵv = LinRange(-0.001, 10.0, 11)
function foo_stab(δ_loc, ϵ_loc)
    return log(spectralradius(dp_fast; p=(ζ, δ_loc, ϵ_loc, b_chart, τ, T)))
end
mymdbm = MDBM_Problem(foo_stab, [MDBM.Axis(δv, "δ"), MDBM.Axis(ϵv, "ϵ")])
MDBM.solve!(mymdbm, 4, interpolationorder=0)
MDBM.solve!(mymdbm, 0, interpolationorder=1)
x_mdbm, y_mdbm = getinterpolatedsolution(mymdbm)

# Final Stability Chart Figure
fig_chart = Figure(size=(1200, 600))
ax_amp = GLMakie.Axis(fig_chart[1, 1], title="Amplitude + MDBM Boundary", xlabel="δ", ylabel="ϵ")
ax_rho = GLMakie.Axis(fig_chart[1, 2], title="Spectral Radius + MDBM Boundary", xlabel="δ", ylabel="ϵ")

# Mask unstable regions for amplitude
Amp_masked = copy(Amp_chart)
Amp_masked[Spec_chart.>1.0] .= NaN

hm1 = heatmap!(ax_amp, δv, ϵv, log10.(Amp_masked' .+ 1e-6), colormap=:viridis)
scatter!(ax_amp, x_mdbm, y_mdbm, color=:red, markersize=4)
Colorbar(fig_chart[2, 1], hm1, vertical=false, label="log10(Amplitude)")

Spec_sat = copy(Spec_chart)
Spec_sat[Spec_sat.>1.0] .= 1.0
hm2 = heatmap!(ax_rho, δv, ϵv, log.(Spec_sat'), colormap=:inferno)
scatter!(ax_rho, x_mdbm, y_mdbm, color=:black, markersize=4)
Colorbar(fig_chart[2, 2], hm2, vertical=false, label="log(ρ)")

#save("examples/01_mathieu_stability_chart.png", fig_chart)
display(fig_chart)

println("Analysis complete.")
