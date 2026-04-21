#]add https://github.com/fodorgera/SddeToOde.jl.git
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
    a1 = 0.2#0.17
    δ = 1.0
    b0 = 0.1#-0.2
    ε = 1.0
    c0 = 0.15
    σ0 = 0.05
    ω = 1.0

    n = 2
    AM(t) = [0 1.0; -(δ + ε * cos(ω * t)) -a1]
    BM(t) = [0 0; b0 0]
    cM(t) = [0.0, c0]

    αM(t) = [0 0.0; -σ0*(δ+ε*cos(ω * t)) -σ0*a1]
    βM(t) = σ0 * BM(t)
    γM(t) = [0.0, σ0]

    τ = 2 * π
end;


# 3. Create Figure Layout
fig = Figure(size=(1200, 900))
ax_time = GLMakie.Axis(fig[1, 1:2], title="Long Simulation Time Series", xlabel="Time (t)", ylabel="x(t)")
ax_orbit = GLMakie.Axis(fig[2, 1], title="Periodic Orbit Comparison", xlabel="x", ylabel="dx/dt")
ax_complex = GLMakie.Axis(fig[2, 2], title="Floquet Multipliers (Complex Plane)", xlabel="Re(μ)", ylabel="Im(μ)", aspect=GLMakie.DataAspect())


# deterministic history φ(t) for t≤0
φ(t) = [1.0, 1.0];
T_sim = τ * 1300;
T_sim = τ * 40;
mval = 10
for mval in 10#6:2:10
    println

    mathieu_stoch, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=T_sim, φ=φ, m=mval, dde=false)
    Solver_args = Dict(:alg => Tsit5(), :verbose => false, :reltol => 1e-2, :abstol => 1e-2)#, :adaptive => false)#
    #Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-3)#
    @time mathieu_sol = solve(mathieu_stoch; Solver_args...)
    mathieu_res = [SddeToOde.get_x_moments(mathieu_sol, mathieu_dde_meta, t) for t in mathieu_sol.t]


    mathieu_avgs = [r[1][1] for r in mathieu_res]
    mathieu_vars = [r[2][1, 1] for r in mathieu_res]

    lines!(ax_time, mathieu_sol.t, [u[1] for u in mathieu_sol.u], linewidth=1)#, color=:blue

    # Extract steady-state periodic orbit (last period)
    t_end = mathieu_sol.t[end]
    t_orbit = range(t_end - τ, t_end, length=2000)
    sol_steady = [mathieu_sol(t) for t in t_orbit]
    lines!(ax_orbit, getindex.(sol_steady, 1), getindex.(sol_steady, 2), linewidth=3, label="Long Sim (Steady State)")#

    display(fig) # Show intermediate progress
end
## ---------------- Affine mapping ---------------------------

using AffineMapDiffEq
using KrylovKit
Neig = 12#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-2, krylovdim=Neig + 20, verbosity=2, maxiter=10));
Solver_args = Dict(:alg => Tsit5(), :verbose => false, :reltol => 1e-3)#,:abstol => 1e-4
τmax = NaN#τ/20000 #maximal timedelay in the mapping
Nstep = 1 # discretization number of the mapping
Timeperiod = τ # timeperiod of the mapping

mathieu_stoch, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=Timeperiod, φ=φ, m=mval, dde=false);

#Creating the problem
dpMathieu = dynamic_problemSampled(mathieu_stoch, Solver_args, τmax,
    Historyresolution=Nstep,
    zerofixpont=false, affineinteration=1,
    Krylov_arg=Krylov_arg, perturbation_size=1e-4)

p = (0.0,)
#solvig the problem (mu: Floquet multiplier, saff: discretized fixed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
@time mus, saff, sol0 = affine(dpMathieu; p=p);





#@code_warntype  affine(dpMathieu; p=p)
Solver_args = Dict(:alg => Tsit5(), :verbose => false, :reltol => 1e-4, :abstol => 1e-4)
sol0 = solve(remake(mathieu_stoch, u0=saff[1], tspan=(0.0, Timeperiod * 2)); Solver_args...)

mu_vals = mus[1]
@show abs.(mu_vals)
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




# ##
# @warn "eddig működik csak"
# 
# 
# sol0_res = [SddeToOde.get_x_moments(sol0, mathieu_dde_meta, t) for t in sol0.t]
# 
# sol0_avgs = [r[1][1] for r in sol0_res];
# sol0_vars = [r[2][1, 1] for r in sol0_res];
# 
# line!(sol0.t .+ T_sim, [sol0_avgs, sol0_avgs .+ sqrt.(sol0_vars)], lw=5, label=["avg" "avg + std"])
# #plot!(sol0.t .+ T_sim, sol0_avgs, lw = 3,label=["avg"])
# #plot!(sol0.t .+ T_sim, [sol0_avgs, sol0_avgs .+(sol0_vars)], lw = 5,label=["avg" "avg + std"])
# 
# ##--------------------------------------------------








## ------------------ Stability boundary - MDBM -----------------


using AffineMapDiffEq
using KrylovKit
Neig = 2#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-2, krylovdim=Neig + 15, verbosity=0, maxiter=10));
Solver_args = Dict(:alg => Tsit5(), :verbose => false, :reltol => 1e-3)#,:abstol => 1e-4
τmax = NaN#τ/20000 #maximal timedelay in the mapping
Nstep = 1 # discretization number of the mapping
Timeperiod = τ # timeperiod of the mapping

mathieu_stoch, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=Timeperiod, φ=φ, m=mval, dde=false);

#Creating the problem
dpMathieu = dynamic_problemSampled(mathieu_stoch, Solver_args, τmax,
    Historyresolution=Nstep,
    zerofixpont=false, affineinteration=1,
    Krylov_arg=Krylov_arg, perturbation_size=1e-4)


function fooDelay(δ_loc, ϵ_loc)
    begin
        a1 = 0.2#0.17

        #δ = 1.0
        b0 = 0.1#-0.2
        #ε = 1.0
        c0 = 0.15
        σ0 = 0.05
        ω = 1.0

        n = 2
        AM(t) = [0 1.0; -(δ_loc + ϵ_loc * cos(ω * t)) -a1]
        BM(t) = [0 0; b0 0]
        cM(t) = [0.0, c0]

        αM(t) = [0 0.0; -σ0*(δ_loc+ϵ_loc*cos(ω * t)) -σ0*a1]
        βM(t) = σ0 * BM(t)
        γM(t) = [0.0, σ0]

        τ = 2 * π
    end


    # deterministic history φ(t) for t≤0
    φ(t) = [1.0, 1.0]
    T_sim = τ * 1
    mval = 10

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
δ = 1.0
ε = 2.0
fooDelay(δ, ε)
#Brute-force grid evaluation
δv = LinRange(0.0, 5.0, 37)
ϵv = LinRange(-0.001, 5.0, 36)
Spec_chart = zeros(length(ϵv), length(δv))
Amp_chart = zeros(length(ϵv), length(δv))
#Threads.@threads 
@time Threads.@threads for j in 1:length(δv)
    δ_v = δv[j]
    println("Evaluating percentage: $(round(j/length(δv)*100, digits=2))% , at δ = $δ_v")
    for i in 1:length(ϵv)
        ϵ_v = ϵv[i]

        # println("Evaluating percentage: $(round(((j-1)*length(ϵv)+i)/(length(δv)*length(ϵv))*100, digits=2))% , at δ = $δ_v, ϵ = $ϵ_v ")
        println("Evaluating at δ = $δ_v, ϵ = $ϵ_v ")

        mu_loc, s0_loc, _ = fooDelay(δ_v, ϵ_v)
        Spec_chart[i, j] = maximum(abs.(mu_loc[1]))
        Amp_chart[i, j] = norm(s0_loc[1]) # norm of the fixed point
    end
end


# Final Stability Chart Figure
fig_chart = Figure(size=(1200, 600))
ax_amp = GLMakie.Axis(fig_chart[1, 1], title="Amplitude + MDBM Boundary (ODE)", xlabel="δ", ylabel="ϵ")
ax_rho = GLMakie.Axis(fig_chart[1, 2], title="Spectral Radius + MDBM Boundary (ODE)", xlabel="δ", ylabel="ϵ")

Amp_masked = copy(Amp_chart)
Amp_masked[Spec_chart.>1.0] .= NaN

Spec_sat = copy(Spec_chart)
Spec_sat[Spec_sat.>1.0] .= 1.0


hm1 = heatmap!(ax_amp, δv, ϵv, log10.(Amp_masked' .+ 1e-6), colormap=:viridis)
Colorbar(fig_chart[2, 1], hm1, vertical=false, label="log10(Amplitude)")


hm2 = heatmap!(ax_rho, δv, ϵv, log.(Spec_sat'), colormap=:inferno)
Colorbar(fig_chart[2, 2], hm2, vertical=false, label="log(ρ)")
display(fig_chart)
##


using MDBM
ax1 = MDBM.Axis(0.0:1.0:5.0, "δ") # initial grid in x direction
ax2 = MDBM.Axis(-0.001:2.0:5.0, "ϵ") # initial grid in y direction

function fooDelay_mdbm_wrapper(δ_l::Float64, ϵ_l::Float64)::Float64
    mu_loc, _, _ = fooDelay(δ_l, ϵ_l)
    mu_vals = mu_loc[1]
    return log(maximum(abs.(mu_vals)))
end
fooDelay_mdbm_wrapper(4.0, 1.0)
mymdbm = MDBM_Problem(fooDelay_mdbm_wrapper, [ax1, ax2])
for _ in 1:2
    @time MDBM.solve!(mymdbm, 1)
    #points where the function foo was evaluated
    x_mdbm, y_mdbm = getinterpolatedsolution(mymdbm)

    scatter!(ax_amp, x_mdbm, y_mdbm, color=:red, markersize=4)
    scatter!(ax_rho, x_mdbm, y_mdbm, color=:black, markersize=4)

    display(fig_chart)
end


display(fig_chart)