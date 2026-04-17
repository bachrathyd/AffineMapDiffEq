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
    σ0 = 0.2
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
T_sim = τ * 20;
mval = 20

# mathieu_ode, mathieu_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=T_sim, φ=φ, m=mval);
# @time mathieu_sol = solve(mathieu_ode, Tsit5(), dt=τ/1000, adaptive=false);
# mathieu_res = [SddeToOde.get_x_moments(mathieu_sol, mathieu_meta, t) for t in mathieu_sol.t];
# 
# mathieu_avgs = [r[1][1] for r in mathieu_res];
# mathieu_vars = [r[2][1,1] for r in mathieu_res];
# 
# plot(mathieu_sol.t, [mathieu_avgs, mathieu_avgs .+ sqrt.(mathieu_vars)], label=["avg" "avg + std"])


# if you want DDEProblem outpu, use dde=true
mathieu_dde_stoch, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=T_sim, φ=φ, m=mval, dde=true);

Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-5)#, :adaptive => false)#
#Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-3)#
@time mathieu_sol = solve(mathieu_dde_stoch; Solver_args...);
mathieu_res = [SddeToOde.get_x_moments(mathieu_sol, mathieu_dde_meta, t) for t in mathieu_sol.t];


mathieu_avgs = [r[1][1] for r in mathieu_res];
mathieu_vars = [r[2][1, 1] for r in mathieu_res];

lines!(ax_time, mathieu_sol.t, [u[1] for u in mathieu_sol.u], color=:blue, linewidth=1)

# Extract steady-state periodic orbit (last period)
t_end = mathieu_sol.t[end]
t_orbit = range(t_end - τ, t_end, length=2000)
sol_steady = [mathieu_sol(t) for t in t_orbit]
lines!(ax_orbit, getindex.(sol_steady, 1), getindex.(sol_steady, 2), color=:blue, linewidth=3, label="Long Sim (Steady State)")

display(fig) # Show intermediate progress

## ---------------- Affine mapping ---------------------------

using AffineMapDiffEq
using KrylovKit
Neig = 4#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-3, krylovdim=Neig + 10, verbosity=4, maxiter=10));

τmax = τ #maximal timedelay in the mapping
Nstep = 100 # discretization number of the mapping
Timeperiod = τ # timeperiod of the mapping
mathieu_dde_stoch, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=Timeperiod, φ=φ, m=mval, dde=true);

#Creating the problem
dpMathieu = dynamic_problemSampled(mathieu_dde_stoch, Solver_args, τmax,
    Historyresolution=Nstep,
    zerofixpont=false, affineinteration=2,
    Krylov_arg=Krylov_arg, perturbation_size=1e-1)

p = 0.5
#solvig the problem (mu: Floquet multiplier, saff: discretized fixed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
@time mus, saff, sol0 = affine(dpMathieu; p=p);

mu_vals = mus[1]

# Plot Affine Fixed Point (Bottom Left)
lines!(ax_orbit, getindex.(saff, 1), getindex.(saff, 2), color=:red, linestyle=:dash, linewidth=2, label="Affine Fixed Point")
GLMakie.scatter!(ax_orbit, getindex.(saff, 1), getindex.(saff, 2), color=:red, markersize=6)

lines!(ax_orbit, [u[1] for u in sol0.u], [u[2] for u in sol0.u], color=:magenta, linewidth=1, label="Affine Periodic Orbit")
axislegend(ax_orbit)

# Plot Complex Plane (Bottom Right)
# Unit Circle
θ = range(0, 2π, length=100)
lines!(ax_complex, cos.(θ), sin.(θ), color=:black, linestyle=:dash)
# Multipliers
GLMakie.scatter!(ax_complex, real.(mu_vals), imag.(mu_vals), color=:red, markersize=10, label="μ")
vlines!(ax_complex, [0], color=:gray, linewidth=0.5)
hlines!(ax_complex, [0], color=:gray, linewidth=0.5)

display(fig)




##
@warn "eddig működik csak"


sol0_res = [SddeToOde.get_x_moments(sol0, mathieu_dde_meta, t) for t in sol0.t]

sol0_avgs = [r[1][1] for r in sol0_res];
sol0_vars = [r[2][1, 1] for r in sol0_res];

line!(sol0.t .+ T_sim, [sol0_avgs, sol0_avgs .+ sqrt.(sol0_vars)], lw=5, label=["avg" "avg + std"])
#plot!(sol0.t .+ T_sim, sol0_avgs, lw = 3,label=["avg"])
#plot!(sol0.t .+ T_sim, [sol0_avgs, sol0_avgs .+(sol0_vars)], lw = 5,label=["avg" "avg + std"])

##--------------------------------------------------








## ------------------ Stability boundary - MDBM -----------------



using MDBM
ax1 = MDBM.Axis(0.0:1.0:5.0, "δ") # initial grid in x direction
ax2 = MDBM.Axis(-0.001:2.0:10.0, "ϵ") # initial grid in y direction
function fooDelay(δ, ϵ)
    begin
        a1 = 0.2#0.17
        #δ = 1.0
        b0 = 0.1#-0.2
        #ε = 1.0
        c0 = 0.15
        σ0 = 0.2
        ω = 1.0

        n = 2
        AM(t) = [0 1.0; -(δ + ε * cos(ω * t)) -a1]
        BM(t) = [0 0; b0 0]
        cM(t) = [0.0, c0]

        αM(t) = [0 0.0; -σ0*(δ+ε*cos(ω * t)) -σ0*a1]
        βM(t) = σ0 * BM(t)
        γM(t) = [0.0, σ0]

        τ = 2 * π
    end


    # deterministic history φ(t) for t≤0
    φ(t) = [1.0, 1.0]
    T_sim = τ * 1
    mval = 20

    τmax = τ
    mathieu_dde_stoch, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=T_sim, φ=φ, m=mval, dde=true)



    #Creating the problem
    dpMathieu = dynamic_problemSampled(mathieu_dde_stoch, Solver_args, τmax,
        Timeperiod; Historyresolution=Nstep,
        zerofixpont=false, affineinteration=1,
        Krylov_arg=Krylov_arg)

    p = 0.5
    #solvig the problem (mu: Floquet multiplier, saff: discretized fixed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
    @time mu, saff, sol0 = affine(dpMathieu; p=p)
end

δ = 1.0
ε = 1.0
fooDelay(δ, ε)

mymdbm = MDBM_Problem(fooDelay, [ax1, ax2])
@time MDBM.solve!(mymdbm, 1)
#points where the function foo was evaluated
x_eval, y_eval = getevaluatedpoints(mymdbm)
x_sol, y_sol = getinterpolatedsolution(mymdbm)
#scatter(x_eval,y_eval,markersize=1)
#scatter!(x_sol,y_sol,markersize=2)
scatter!(x_sol, y_sol, markersize=1)
