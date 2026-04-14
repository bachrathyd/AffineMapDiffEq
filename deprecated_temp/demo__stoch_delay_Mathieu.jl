#]add https://github.com/fodorgera/SddeToOde.jl.git
5 + 5
using DifferentialEquations
using Plots

plotly()
using StaticArrays
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

# deterministic history φ(t) for t≤0
φ(t) = [1.0, 1.0];
T_sim = τ * 1300;
T_sim = τ * 20;
mval = 10

# mathieu_ode, mathieu_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=T_sim, φ=φ, m=mval);
# @time mathieu_sol = solve(mathieu_ode, Tsit5(), dt=τ/1000, adaptive=false);
# mathieu_res = [SddeToOde.get_x_moments(mathieu_sol, mathieu_meta, t) for t in mathieu_sol.t];
# 
# mathieu_avgs = [r[1][1] for r in mathieu_res];
# mathieu_vars = [r[2][1,1] for r in mathieu_res];
# 
# plot(mathieu_sol.t, [mathieu_avgs, mathieu_avgs .+ sqrt.(mathieu_vars)], label=["avg" "avg + std"])


# if you want DDEProblem outpu, use dde=true
mathieu_dde, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=T_sim, φ=φ, m=mval, dde=true);

Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-6, :dt => τ / 1000, :adaptive => false)#
Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-2)#
@time mathieu_sol = solve(mathieu_dde; Solver_args...);
mathieu_res = [SddeToOde.get_x_moments(mathieu_sol, mathieu_dde_meta, t) for t in mathieu_sol.t];


mathieu_avgs = [r[1][1] for r in mathieu_res];
mathieu_vars = [r[2][1, 1] for r in mathieu_res];

plot(mathieu_sol.t, [mathieu_avgs, mathieu_avgs .+ sqrt.(mathieu_vars)], label=["avg" "avg + std"])

## ---------------- Affine mapping ---------------------------

using AffineMapDiffEq
using KrylovKit
Neig = 10#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-42, krylovdim=Neig + 4, verbosity=0, maxiter=10));

τmax = τ #maximal timedelay in the mapping
Nstep = 3 # discretization number of the mapping
Timeperiod = τ # timeperiod of the mapping
mathieu_dde, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=Timeperiod, φ=φ, m=mval, dde=true);

#Creating the problem
dpMathieu = dynamic_problemSampled(mathieu_dde, Solver_args, τmax,
    Historyresolution=Nstep,
    zerofixpont=false, affineinteration=1,
    Krylov_arg=Krylov_arg)

p = 0.5
#solvig the problem (mu: Floquet multiplier, saff: discretized fixed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
@time mu, saff, sol0 = affine(dpMathieu; p=p);
@time mu, saff, sol0 = affine(dpMathieu; p=p);




sol0_res = [SddeToOde.get_x_moments(sol0, mathieu_dde_meta, t) for t in sol0.t]

sol0_avgs = [r[1][1] for r in sol0_res];
sol0_vars = [r[2][1, 1] for r in sol0_res];

plot!(sol0.t .+ T_sim, [sol0_avgs, sol0_avgs .+ sqrt.(sol0_vars)], lw=5, label=["avg" "avg + std"])
#plot!(sol0.t .+ T_sim, sol0_avgs, lw = 3,label=["avg"])
#plot!(sol0.t .+ T_sim, [sol0_avgs, sol0_avgs .+(sol0_vars)], lw = 5,label=["avg" "avg + std"])

##--------------------------------------------------




# Plotting the Floquet multipliers
#in log scale
plot(log.(abs.(mu[1])))
#in complex plane
scatter(mu[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))






## ------------------ Stability boundary - MDBM -----------------



using MDBM
ax1 = MDBM.Axis(0.0:1.0:5.0, "δ") # initial grid in x direction
ax2 = MDBM.Axis(-0.001:2.0:10.0, "ϵ") # initial grid in y direction
function fooDelay(δ, ϵ)
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
    end


    # deterministic history φ(t) for t≤0
    φ(t) = [1.0, 1.0]
    T_sim = τ * 1
    mval = 20

    τmax = τ
    mathieu_dde, mathieu_dde_meta = SddeToOde.get_ode_from_sdde(AM, BM, cM, αM, βM, γM; τ=τ, T=T_sim, φ=φ, m=mval, dde=true)



    #Creating the problem
    dpMathieu = dynamic_problemSampled(mathieu_dde, Solver_args, τmax,
        Timeperiod; Historyresolution=Nstep,
        zerofixpont=false, affineinteration=1,
        Krylov_arg=Krylov_arg)

    p = 0.5
    #solvig the problem (mu: Floquet multiplier, saff: discretized fixed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
    @time mu, saff, sol0 = affine(dpMathieu; p=p)
end

    δ = 1.0
    ε = 1.0
fooDelay(δ,ε )

mymdbm = MDBM_Problem(fooDelay, [ax1, ax2])
@time MDBM.solve!(mymdbm, 1)
#points where the function foo was evaluated
x_eval, y_eval = getevaluatedpoints(mymdbm)
x_sol, y_sol = getinterpolatedsolution(mymdbm)
#scatter(x_eval,y_eval,markersize=1)
#scatter!(x_sol,y_sol,markersize=2)
scatter!(x_sol, y_sol, markersize=1)
