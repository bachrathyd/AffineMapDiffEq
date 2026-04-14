# Example 01: Delayed Mathieu Equation with Damping
# Stability analysis using Affine Mapping (Optimized Version)

using Pkg; Pkg.activate(".")
using AffineMapDiffEq
using StaticArrays
using DifferentialEquations
using LinearAlgebra
using KrylovKit
using Plots

# 1. Define the Governing Equation (In-place for maximum performance)
function DelayMathieu!(du, u, h, p, t)
    # Parameters
    ζ, δ, ϵ, b, τ, T = p
    # External forcing
    F = 0.1 * (cos(2π * t / T)^10)
    # Components of the delayed differential equation
    du[1] = u[2]
    du[2] = -(δ + ϵ * cos(2π * t / T)) * u[1] - 2ζ * u[2] + b * h(p, t - τ)[1] + F
end

# 2. Setup Parameters
ζ = 0.02          # damping coefficient
δ = 1.5           # natural frequency
ϵ = 0.15          # excitation coefficient
τ = 2π            # time delay
b = 0.5           # delay gain
T = 2π            # time period
p = (ζ, δ, ϵ, b, τ, T)

# 3. Initial Simulation
u0 = [1.0, 0.0]
h(p, t) = [0.0, 0.0]
probMathieu = DDEProblem{true}(DelayMathieu!, u0, h, (0.0, T * 500.0), p; constant_lags=[τ])

## Solver configuration
Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-3)
println("Running long simulation...")
@time sol = solve(probMathieu; Solver_args...)
plot(sol, title="Mathieu Equation Time Series")

# 4. Affine Mapping for Fixed Point and Stability
println("Calculating Affine Mapping...")
Neig = 8
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=18+5, verbosity=0))

τmax = τ
Nstep = 100
probMapping = DDEProblem{true}(DelayMathieu!, u0, h, (0.0, T), p; constant_lags=[τ])

dpMathieu = dynamic_problemSampled(probMapping, Solver_args, τmax; 
    Historyresolution=Nstep,
    zerofixpont=false, affineinteration=2,
    Krylov_arg=Krylov_arg)

# Solve for Floquet multipliers (mu) and the discretized fixed point (saff)
# The first call includes compilation time.
@time mu, saff, sol0 = affine(dpMathieu; p= (ζ, δ, ϵ, b, τ, T));
#@profview affine(dpMathieu; p= (ζ, δ, ϵ, b+0.2, τ, T));
println("Max Floquet multiplier: ", maximum(abs.(mu[1])))

# 5. Plot Results
# Plot Floquet multipliers in the complex plane
p_mu = scatter(mu[1], label="Floquet Multipliers", aspect_ratio=:equal)
plot!(p_mu, sin.(0:0.01:2π), cos.(0:0.01:2π), label="Unit Circle", color=:black, linestyle=:dash)
title!(p_mu, "Floquet Multipliers")
display(p_mu)


#last period of the long simulation:
t_select_period=0.0:0.01:T
t_select_delay=eriod=0.0:0.01:τ
sol_period=sol(sol.t[end] .- t_select_period)
sol_delay=sol(sol.t[end] .- t_select_delay)
# Comparing the solutions:
plot(sol_period[1, :], sol_period[2, :], marker=:circle,markersize=6,lab="")
plot!(getindex.(saff,1), getindex.(saff,2), marker=:circle,markersize=4,lab="")#
plot!(sol0[1, :], sol0[2, :], marker=:circle,lw = 0,lab="")#marker=:cross,markersize=2)#





## ----------------------- creating stability chart -------------------------

Krylov_arg=(1,:LM, KrylovKit.Arnoldi(tol=1e-12,krylovdim=8+5,verbosity=0));
Nstep=25
dpMathieu = dynamic_problemSampled(probMapping, Solver_args, τmax; Historyresolution=Nstep,
    zerofixpont=false,    affineinteration=1,
    Krylov_arg=Krylov_arg)

## ------------------ Peak-2-Peak ampNorm color, and spectral radius - BRUTE FORCE in the plane of δ-ϵ  ----------------

b=0.05 # delay parameter, Note: b=0 provide the traditional stability chart of the Mathieu equations
δv =LinRange( -1.0,5.0,40) # initial grid in x direction
ϵv = LinRange(-0.001,10.0,41) # initial grid in y direction
Aaff = zeros(size(ϵv, 1), size(δv, 1))
Spek_aff = zeros(size(ϵv, 1), size(δv, 1))
@time Threads.@threads  for j in 1:size(δv, 1)
    @show j
    @inbounds δ = δv[j]
    Threads.@threads         for i in 1:size(ϵv, 1)
        @inbounds ϵ = ϵv[i]
       muaff, s0aff, sol0 = affine(dpMathieu; p=(ζ, δ, ϵ, b, τ, T))
            Aaff[i, j] = norm(getindex.(s0aff, 1)) # norm of the motion
        # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1))) #maximum amplitude of the orbit
        Spek_aff[i, j] = maximum(abs.(muaff[1]))
    end
end

#Plotting the maximal amplitud on the stable domain only
Aaffsat = deepcopy(Aaff);
Aaffsat[Spek_aff.>1.0] .= 1.0;#eliminate the positions of instable case
heatmap(δv, ϵv, log.(Aaffsat))



#Plotting the maximal amplitud on the stable domain only
Spek_affsat = deepcopy(Spek_aff);
Spek_affsat[Spek_affsat.>1.0] .= 1.0;
heatmap(δv, ϵv, log.(Spek_affsat))




## ------------------ Stability boundary - MDBM -----------------
using MDBM
b=0.05;
ax1 = MDBM.Axis(0.0:1.0:5.0, "δ") # initial grid in x direction
ax2 = MDBM.Axis(-0.001:2.0:10.0, "ϵ") # initial grid in y direction
function fooDelay(δ, ϵ)
    ABSmuMax = spectralradius(dpMathieu; p=(ζ, δ, ϵ, b, τ, T))
    #return ABSmuMax - 1.0
    return log(ABSmuMax)
end
# @benchmark fooDelay(2.0, 2.0)
mymdbm = MDBM_Problem(fooDelay, [ax1, ax2])
@time MDBM.solve!(mymdbm, 5)
#points where the function foo was evaluated
x_eval, y_eval = getevaluatedpoints(mymdbm)
x_sol, y_sol = getinterpolatedsolution(mymdbm)
#scatter(x_eval,y_eval,markersize=1)
#scatter!(x_sol,y_sol,markersize=2)
scatter!(x_sol, y_sol, markersize=1)
