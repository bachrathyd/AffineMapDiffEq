#]add https://github.com/fodorgera/SddeToOde.jl.git
#include("02_delayed_mathieu_amplitude.jl")
5 + 5
using DifferentialEquations
using StaticArrays


using GLMakie
# MATHIEU TEST

# using Pkg
# Pkg.develop("SddeToOde")

# develop SddeToOde
# activate SddeToOde
# add StaticArrays
using SddeToOde
begin

    const a1 = 0.02# 2ζ
    const δ = 1.5
    const ε = 0.15
    const τ = 2π
    const b0 = 0.05#-0.2
    const c0 = 0.15
    const σ0 = 0.1
    const ω = 1.0

    n = 2
    AM(t, yy, xx) = [0 1.0; -(δ + ε * cos(ω * t)) -a1]
    BM(t, yy, xx) = [0 0; b0 0]
    cM(t, yy, xx) = [0.0, c0]

    αM(t, yy, xx) = [0 0.0; -σ0*(δ+ε*cos(ω * t)) -σ0*a1]
    βM(t, yy, xx) = σ0 * BM(t, yy, xx)
    γM(t, yy, xx) = [0.0, σ0]

end;

# 3. Create Figure Layout
fig = Figure(size=(1200, 900))
ax_time = GLMakie.Axis(fig[1, 1:2], title="Long Simulation Time Series", xlabel="Time (t)", ylabel="x(t)")
ax_orbit = GLMakie.Axis(fig[2, 1], title="Periodic Orbit Comparison", xlabel="x", ylabel="dx/dt")
ax_complex = GLMakie.Axis(fig[2, 2], title="Floquet Multipliers (Complex Plane)", xlabel="Re(μ)", ylabel="Im(μ)", aspect=GLMakie.DataAspect())

# deterministic history φ(t) for t≤0
φ(t) = [1.0, 1.0];
T_sim = τ * 100;
mval = 15

println("Running SddeToOde with m = $mval...")
mathieu_stoch, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=T_sim, φ=φ, m=mval, dde=false)
Solver_args = Dict(:alg => Tsit5(), :verbose => false, :reltol => 1e-4)#, :adaptive => false)#
#Solver_args = Dict(:alg => Tsit5(), :verbose => false, :adaptive => false, :dt => τ / 50)#, )#

@time mathieu_sol = solve(mathieu_stoch; Solver_args...);

mathieu_res = [SddeToOde.get_x_moments(mathieu_sol, mathieu_dde_meta, t) for t in mathieu_sol.t]
mathieu_avgs = [r[1][1] for r in mathieu_res]
mathieu_vars = [r[2][1, 1] for r in mathieu_res]

lines!(ax_time, mathieu_sol.t, mathieu_avgs, linewidth=1)#, color=:blue
band!(ax_time, mathieu_sol.t, mathieu_avgs .- sqrt.(mathieu_vars), mathieu_avgs .+ sqrt.(mathieu_vars), color=(:blue, 0.2), label="±1 std")
lines!(ax_time, mathieu_sol.t, mathieu_avgs .+ sqrt.(mathieu_vars), linewidth=1)#, color=:blue
lines!(ax_time, mathieu_sol.t, mathieu_avgs .- sqrt.(mathieu_vars), linewidth=1)#, color=:blue

# Extract steady-state periodic orbit (last period)
t_end = mathieu_sol.t[end]
t_orbit = range(t_end - τ, t_end, length=2000)
sol_steady = [mathieu_sol(t) for t in t_orbit]
lines!(ax_orbit, getindex.(sol_steady, 1), getindex.(sol_steady, 2), linewidth=3, label="Long Sim (Steady State)")#
lines!(ax_orbit, getindex.(sol_steady, 1), getindex.(sol_steady, 2), linewidth=3, label="Long Sim (Steady State)")#

display(fig) # Show intermediate progress

##

using BenchmarkTools
using LinearAlgebra

mvals_collected = Observable(Int[])
solvers_to_test = [
    ("Tsit5", Dict(:alg => Tsit5(), :verbose => false, :reltol => 1e-4)),
    ("BS3", Dict(:alg => BS3(), :verbose => false, :reltol => 1e-4)),
    ("RK4", Dict(:alg => RK4(), :verbose => false, :reltol => 1e-4))
]#,
#   ("Rosenbrock23", Dict(:alg => Rosenbrock23(autodiff=false), :verbose => false, :reltol => 1e-4))
#]

times_collected = [Observable(Float64[]) for _ in solvers_to_test]
trend_collected = [Observable(Float64[]) for _ in solvers_to_test]
labels_measured = [Observable("$name") for (name, _) in solvers_to_test]
labels_trend = [Observable("$name O(?)") for (name, _) in solvers_to_test]

fig_perf = Figure(size=(800, 600))
ax_perf = GLMakie.Axis(fig_perf[1, 1], title="Solve Time vs m", xlabel="m", ylabel="Time (s)", yscale=log10, xscale=log10)

colors = [:blue, :green, :red, :orange, :purple]
for i in 1:length(solvers_to_test)
    scatterlines!(ax_perf, mvals_collected, times_collected[i], color=colors[i], markersize=8, label=labels_measured[i])
    lines!(ax_perf, mvals_collected, trend_collected[i], color=colors[i], linestyle=:dash, linewidth=2, label=labels_trend[i])
end
axislegend(ax_perf, position=:lt)
display(fig_perf)

for mval in 10:1:30
    println("Running SddeToOde with m = $mval...")
    mathieu_stoch, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=T_sim, φ=φ, m=mval, dde=false)

    push!(mvals_collected[], mval)

    for (i, (sname, Solver_args)) in enumerate(solvers_to_test)
        # Perform time measurements with max 1 seconds
        #t_elapsed = @belapsed solve($mathieu_stoch; $(Solver_args)...) seconds=0.25
        t_elapsed = @elapsed solve(mathieu_stoch; (Solver_args)...)
        println("  Elapsed time for $sname, m=$mval: $t_elapsed s")

        # Collect and plot the data
        push!(times_collected[i][], max(1e-12, t_elapsed))

        if length(mvals_collected[]) >= 2
            log_m = log.(mvals_collected[])
            log_t = log.(times_collected[i][])
            A = hcat(ones(length(log_m)), log_m)
            b, a = A \ log_t

            trend_collected[i][] = exp.(b .+ a .* log_m)
            labels_trend[i][] = "$sname O(m^$(round(a, digits=2)))"
            println("  Fitted trend for $sname: O(m^$(round(a, digits=2)))")
        else
            trend_collected[i][] = times_collected[i][]
        end

        notify(times_collected[i])
        notify(trend_collected[i])
        notify(labels_trend[i])
    end

    notify(mvals_collected)
    autolimits!(ax_perf)
    display(fig_perf)
    yield()
end

## ---------------- Affine mapping ---------------------------
using AffineMapDiffEq
using KrylovKit

Neig = 10#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-3, krylovdim=(Neig + 25), verbosity=6, maxiter=10));
Solver_args = Dict(:alg => Tsit5(), :verbose => false, :reltol => 1e-4)#,:abstol => 1e-4
Solver_args = Dict(:alg => BS3(), :verbose => false, :reltol => 1e-3)#,:abstol => 1e-4
Solver_args = Dict(:alg => RK4(), :verbose => false, :reltol => 1e-3)#,:abstol => 1e-4
@show mval
#Solver_args = Dict(:alg => Tsit5(), :verbose => false, :adaptive => false, :dt => τ / 50)

τmax = NaN#τ/20000 #maximal timedelay in the mapping
Nstep = 1 # discretization number of the mapping
Timeperiod = τ *3# timeperiod of the mapping
mval=30
mathieu_stoch, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=Timeperiod, φ=φ, m=mval, dde=false);

#Creating the problem
dpMathieu = dynamic_problemSampled(mathieu_stoch, Solver_args, τmax,
    Historyresolution=Nstep,
    zerofixpont=false, affineinteration=1,
    Krylov_arg=Krylov_arg, perturbation_size=1e-4)

p = (0.0,)
#solvig the problem (mu: Floquet multiplier, saff: discretized fixed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
println("Calculating Affine results...")
@time mus, saff, sol0 = affine(dpMathieu; p=p);

mu_vals = mus[1]
@show abs.(mu_vals)
#@profview  affine(dpMathieu; p=p)



#@code_warntype  affine(dpMathieu; p=p)
#Solver_args = Dict(:alg => Tsit5(), :verbose => false, :adaptive => false, :dt => τ / 50)
sol0 = solve(remake(mathieu_stoch, u0=saff[1], tspan=(0.0, Timeperiod * 2)); Solver_args...)

# Plot Affine Fixed Point (Bottom Left)
lines!(ax_orbit, getindex.(saff, 1), getindex.(saff, 2), color=:red, linestyle=:dash, linewidth=2, label="Affine Fixed Point")
GLMakie.scatter!(ax_orbit, getindex.(saff, 1), getindex.(saff, 2), markersize=6)#, color=:red

lines!(ax_orbit, [u[1] for u in sol0.u], [u[2] for u in sol0.u], linewidth=1, label="Affine Periodic Orbit")#, color=:magenta
axislegend(ax_orbit)

# Plot Complex Plane (Bottom Right)
# Unit Circle
θ = range(0, 2π, length=100)
lines!(ax_complex, cos.(θ), sin.(θ), color=:black, linestyle=:dash)
# Multipliers
GLMakie.scatter!(ax_complex, real.(mu_vals), imag.(mu_vals), markersize=10, label="μ_max = $(maximum(abs.(mu_vals)))")#, color=:red
vlines!(ax_complex, [0], color=:gray, linewidth=0.5)
hlines!(ax_complex, [0], color=:gray, linewidth=0.5)
axislegend(ax_complex)
display(fig)




## ------------------ Stability boundary - MDBM -----------------


using AffineMapDiffEq
using KrylovKit
Neig = 1#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-2, krylovdim=Neig + 11, verbosity=0, maxiter=10));
reltol = 1e-4
Solver_args = Dict(:alg => Tsit5(), :verbose => false, :reltol => reltol)#,:abstol => 1e-4
Solver_args = Dict(:alg => BS3(), :verbose => false, :reltol => reltol)#,:abstol => 1e-4
Solver_args = Dict(:alg => RK4(), :verbose => false, :reltol => reltol)#,:abstol => 1e-4
#Solver_args = Dict(:alg => Tsit5(), :verbose => false, :adaptive => false, :dt => τ / 50)#, )#
τmax = NaN#τ/20000 #maximal timedelay in the mapping
Nstep = 1 # discretization number of the mapping
Timeperiod = τ*3 # timeperiod of the mapping

mathieu_stoch, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=Timeperiod, φ=φ, m=mval, dde=false);

#Creating the problem
dpMathieu = dynamic_problemSampled(mathieu_stoch, Solver_args, τmax,
    Historyresolution=Nstep,
    zerofixpont=false, affineinteration=1,
    Krylov_arg=Krylov_arg, perturbation_size=1e-4)


function fooDelay(δ_loc, ϵ_loc)
    begin

        n = 2
        AM(t, yy, xx) = [0 1.0; -(δ_loc + ϵ_loc * cos(ω * t)) -a1]
        BM(t, yy, xx) = [0 0; b0 0]
        cM(t, yy, xx) = [0.0, c0]

        αM(t, yy, xx) = [0 0.0; -σ0*(δ_loc+ϵ_loc*cos(ω * t)) -σ0*a1]
        βM(t, yy, xx) = σ0 * BM(t, yy, xx)
        γM(t, yy, xx) = [0.0, σ0]

    end


    # deterministic history φ(t) for t≤0
    φ(t) = [1.0, 1.0]
    T_sim = τ * 1
    #mval = 10

    τmax = τ
    mathieu_stoch, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=T_sim, φ=φ, m=mval, dde=false)



    #Creating the problem
    dpMathieu = dynamic_problemSampled(mathieu_stoch, Solver_args, τmax,
        Historyresolution=Nstep,
        zerofixpont=false, affineinteration=1,
        Krylov_arg=Krylov_arg, perturbation_size=1e-4)

    p = (0.0,)
    #solvig the problem (mu: Floquet multiplier, saff: discretized fixed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
    mus, saff, sol0 = affine(dpMathieu; p=p)
    #@time 
    #mu_vals = mus[1]
    #@show abs.(mu_vals)
    return mus, saff, sol0
end

using LinearAlgebra

fooDelay(1.0, 2.0)
#Brute-force grid evaluation
const δv = LinRange(0.0, 5.0, 67*2)
const ϵv = LinRange(-0.001, 5.0, 66*2)
Spec_chart_stoc = zeros(length(ϵv), length(δv))
Amp_chart_stoc = zeros(length(ϵv), length(δv))
#Threads.@threads 
@time Threads.@threads for j in 1:length(δv)
    δ_v = δv[j]
    # println("Evaluating percentage: $(round(j/length(δv)*100, digits=2))% , at δ = $δ_v")
    for i in 1:length(ϵv)
        ϵ_v = ϵv[i]

        # println("Evaluating percentage: $(round(((j-1)*length(ϵv)+i)/(length(δv)*length(ϵv))*100, digits=2))% , at δ = $δ_v, ϵ = $ϵ_v ")
        println("Evaluating at δ = $δ_v, ϵ = $ϵ_v ")

        mu_loc, s0_loc, _ = fooDelay(δ_v, ϵ_v)
        Spec_chart_stoc[i, j] = maximum(abs.(mu_loc[1]))
        Amp_chart_stoc[i, j] = norm(s0_loc[1]) # norm of the fixed point
    end
end
#
#for _ in 1:100
# while true
println("Plot the actual state")
# Final Stability Chart Figure
# fig_chart = Figure(size=(1200, 600))
fig_chart = Figure(size=(2200, 1200))
ax_amp = GLMakie.Axis(fig_chart[1, 1], title="Amplitude + MDBM Boundary (ODE)", xlabel="δ", ylabel="ϵ")
ax_rho = GLMakie.Axis(fig_chart[1, 2], title="Spectral Radius + MDBM Boundary (ODE)", xlabel="δ", ylabel="ϵ")

Amp_masked_stoc = copy(Amp_chart_stoc)
Amp_masked_stoc[Spec_chart_stoc.>1.0] .= NaN
Amp_masked_stoc[Amp_masked_stoc.>100.0] .= 100.0

Spec_sat_stoc = copy(Spec_chart_stoc)
Spec_sat_stoc[Spec_sat_stoc.>1.0] .= 1.0


hm1 = heatmap!(ax_amp, δv, ϵv, log10.(Amp_masked_stoc' .+ 1e-6), colormap=:viridis)
Colorbar(fig_chart[2, 1], hm1, vertical=false, label="log10(Amplitude)")

contour!(ax_amp, δv, ϵv, log10.(Amp_masked_stoc' .+ 1e-6), levels=-1:0.25:2, labels=true, color=:black)

hm2 = heatmap!(ax_rho, δv, ϵv, log.(Spec_sat_stoc'), colormap=:inferno)
Colorbar(fig_chart[2, 2], hm2, vertical=false, label="log(ρ)")



#try
#    scatter!(ax_amp, x_mdbm, y_mdbm, color=:red, markersize=4)
#    scatter!(ax_rho, x_mdbm, y_mdbm, color=:blue, markersize=4)
#catch e
#end
display(fig_chart)

#     sleep(2.51)
# end

#save("stoch_mathieu_analysis_sig_$(σ0)_$(mval)_additive only.png", fig_chart)
save("stoch_mathieu_analysis_sig_$(σ0)_$(mval)_multip.png", fig_chart)
##


using MDBM

ax1 = MDBM.Axis(LinRange(0.0, 5.0, 10), "δ") # initial grid in x direction
ax2 = MDBM.Axis(LinRange(-0.001, 5.0, 11), "ϵ") # initial grid in y direction

function fooDelay_mdbm_wrapper(δ_l::Float64, ϵ_l::Float64)::Float64
    mu_loc, _, _ = fooDelay(δ_l, ϵ_l)
    mu_vals = mu_loc[1]
    return log(maximum(abs.(mu_vals)))
end
fooDelay_mdbm_wrapper(4.0, 1.0)
mymdbm_stoc = MDBM_Problem(fooDelay_mdbm_wrapper, [ax1, ax2])
#for _ in 1:4
@time MDBM.solve!(mymdbm_stoc, 3)
#points where the function foo was evaluated
x_mdbm_stoc, y_mdbm_stoc = getinterpolatedsolution(mymdbm_stoc)

scatter!(ax_amp, x_mdbm_stoc, y_mdbm_stoc, color=:blue, markersize=6)
scatter!(ax_rho, x_mdbm_stoc, y_mdbm_stoc, color=:blue, markersize=6)

display(fig_chart)
#end


display(fig_chart)


save("stoch_mathieu_analysis_sig_$(σ0)_$(mval)_multip.png", fig_chart)