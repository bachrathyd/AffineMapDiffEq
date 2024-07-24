#Demo: Delayed Mathieu equation with damping

using Revise
using DDE_mapping

using BenchmarkTools
using Plots
plotly()
using Profile
using StaticArrays
using DifferentialEquations

using LinearAlgebra

using MDBM


# Governing equation

function DelayMathieu(u, h, p, t)
    # Parameters
    ζ, δ, ϵ, b, τ, T = p
    #External forcing
    F = 0.1 * (cos(2pi * t / T) .^ 10)
    # Components of the delayed differential equation
    dx = u[2]
    ddx = -(δ + ϵ * cos(2pi * t / T)) * u[1] - 2 * ζ * u[2] + b * h(p, t - τ)[1] + F
    # Update the derivative vector
    SA[dx, ddx]
end
Base.:+(a::SVector, b::Bool) = a .+ b
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?


## parameters 
ζ = 0.02          # damping coefficient
δ = 1.5#0.2          # nat. freq
ϵ = 0.15#4#5#8;#5         # cut.coeff
τ = 2pi          # Time delay
b = 0.5
T = 2pi
p = ζ, δ, ϵ, b, τ, T
#p = (ζ, ωn, k, τ,10.0)

# test simulation ---------------
#initial condition
u0 = SA[1.0, 0.0]
#history function
h(p, t) = SA[0.0; 0.0]
probMathieu = DDEProblem(DelayMathieu, u0, h, (0.0, T * 1000.0), p; constant_lags=[τ])

#Parameters for the solver as a Dict (it is necessary to collect it for later use)
#Solver_args = Dict(:alg => MethodOfSteps(BS3()), :callback => cb, :adaptive => true, :dt => 0.01, :verbose => false, :reltol => 1e-9)#, save_everystep=false)#abstol,reltol)
Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-4)#
Solver_args = Dict()#
sol = solve(probMathieu; Solver_args...)#abstol,reltol

plot(sol)

#last period of the long simulation:
t_select_period=0.0:0.01:T
t_select_delay=eriod=0.0:0.01:τ
sol_period=sol(sol.t[end] .- t_select_period)
sol_delay=sol(sol.t[end] .- t_select_delay)

#plot the state phase of the last segment
plot(sol_period.t,getindex.(sol_period.u,1))
plot!(sol_period.t,getindex.(sol_period.u,2))
plot!(sol_delay.t,getindex.(sol_delay.u,1))
plot!(sol_delay.t,getindex.(sol_delay.u,2))
#plot the phase space (u - du)
plot(getindex.(sol_delay.u,1),getindex.(sol_delay.u,2))
plot!(getindex.(sol_period.u,1),getindex.(sol_period.u,2))


# ---------------- Affine mapping ---------------------------
using KrylovKit
Neig=8#number of required eigen values
Krylov_arg=(Neig,:LM, KrylovKit.Arnoldi(tol=1e-12,krylovdim=8+5,verbosity=0));

τmax=τ #maximal timedelay in the mapping
Nstep = 100 # discretization number of the mapping
Timeperiod=T # timeperiod of the mapping

#Creating the problem
dpMathieu = dynamic_problemSampled(probMathieu, Solver_args, τmax,
Timeperiod; Historyresolution=Nstep,
    zerofixpont=false,    affineinteration=1,
    Krylov_arg=Krylov_arg)


#solvig the problem (mu: Floquet multiplier, saff: discretized foxed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
@time mu, saff, sol0 = affine(dpMathieu; p=p);


# Comparing the solutions:
plot(sol_period[1, :], sol_period[2, :], marker=:circle,markersize=6,lab="")
plot!(getindex.(saff,1), getindex.(saff,2), marker=:circle,markersize=4,lab="")#
plot!(sol0[1, :], sol0[2, :], marker=:circle,lw = 0,lab="")#marker=:cross,markersize=2)#
##


# Plotting the Floquet multipliers
#in log scale
plot(log.(abs.(mu[1])))
#in complex plane
scatter(mu[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))








# ----------------------- creating stability chart -------------------------

Krylov_arg=(1,:LM, KrylovKit.Arnoldi(tol=1e-12,krylovdim=8+5,verbosity=0));

dpMathieu = dynamic_problemSampled(probMathieu, Solver_args, τmax,
Timeperiod; Historyresolution=Nstep,
    zerofixpont=false,    affineinteration=1,
    Krylov_arg=Krylov_arg)

## ------------------ Peak-2-Peak ampNorm color, and spectral radius - BRUTE FORCE in the plane of δ-ϵ  ----------------

b=0.05 # delay parameter, Note: b=0 provide the traditional stability chart of the Mathieu equations
δv = -1.0:0.051:5.0 # initial grid in x direction
ϵv = -0.001:0.05:10.0 # initial grid in y direction
Aaff = zeros(size(ϵv, 1), size(δv, 1))
Spek_aff = zeros(size(ϵv, 1), size(δv, 1))

@time Threads.@threads  for j in 1:size(δv, 1)
    @inbounds δ = δv[j]
    Threads.@threads     for i in 1:size(ϵv, 1)
        @inbounds ϵ = ϵv[i]
        muaff, s0aff, sol0 = affine(dpMathieu; p=(ζ, δ, ϵ, b, τ, T))
        Aaff[i, j] = norm(getindex.(s0aff, 1)) # norm of the motion
        # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1))) #maximum amplitude of the orbit
        Spek_aff[i, j] = maximum(abs.(muaff[1]))
    end
end

#Plotting the maximal amplitud on the stable domain only
Aaffsat = deepcopy(Aaff);
Aaffsat[Spek_aff.>1.0] .= 0.0;#eliminate the positions of instable case
heatmap(δv, ϵv, log.(Aaffsat))



#Plotting the maximal amplitud on the stable domain only
Spek_affsat = deepcopy(Spek_aff);
Spek_affsat[Spek_affsat.>1.0] .= 0.0;
heatmap(δv, ϵv, log.(Spek_affsat))


#------------------ Stability boundary - MDBM -----------------

b=0.05;
ax1 = Axis(0:1:5, "δ") # initial grid in x direction
ax2 = Axis(-0.001:2.0:10.0, "ϵ") # initial grid in y direction
function fooDelay(δ, ϵ)
    ABSmuMax = spectralradius(dpMathieu; p=(ζ, δ, ϵ, b, τ, T))
    #return ABSmuMax - 1.0
    return log(ABSmuMax)
end
mymdbm = MDBM_Problem(fooDelay, [ax1, ax2])
@time MDBM.solve!(mymdbm, 5)
#points where the function foo was evaluated
x_eval, y_eval = getevaluatedpoints(mymdbm)
x_sol, y_sol = getinterpolatedsolution(mymdbm)
#scatter(x_eval,y_eval,markersize=1)
#scatter!(x_sol,y_sol,markersize=2)
scatter!(x_sol, y_sol, markersize=1)



#--------------------------



## ------------------ Peak-2-Peak ampNorm color, and spectral radius - BRUTE FORCE  in the plane of δ-b-----------------

Krylov_arg=(1,:LM, KrylovKit.Arnoldi(tol=1e-12,krylovdim=8+5,verbosity=0));

dpMathieu = dynamic_problemSampled(probMathieu, Solver_args, τmax,
Timeperiod; Historyresolution=Nstep,
    zerofixpont=false,    affineinteration=1,
    Krylov_arg=Krylov_arg)


δv = 0:0.051:10 # initial grid in x direction
bv = -1.501:0.05:1.5 # initial grid in y direction
Aaff = zeros(size(bv, 1), size(δv, 1))
Spek_aff = zeros(size(bv, 1), size(δv, 1))

@time Threads.@threads  for j in 1:size(δv, 1)
    @inbounds δ = δv[j]
    Threads.@threads     for i in 1:size(bv, 1)
        @inbounds b = bv[i]
        muaff, s0aff, sol0 = affine(dpMathieu; p=(ζ, δ, ϵ, b, τ, T))
        Aaff[i, j] = norm(getindex.(s0aff, 1)) # norm of the motion
        # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1))) #maximum amplitude of the orbit
        Spek_aff[i, j] = maximum(abs.(muaff[1]))
    end
end

#Plotting the maximal amplitud on the stable domain only
Aaffsat = deepcopy(Aaff);
Aaffsat[Spek_aff.>1.0] .= 0.0;#eliminate the positions of instable case
heatmap(δv, bv, log.(Aaffsat))



#Plotting the maximal amplitud on the stable domain only
Spek_affsat = deepcopy(Spek_aff);
Spek_affsat[Spek_affsat.>1.0] .= 0.0;
heatmap(δv, bv, log.(Spek_affsat))


#------------------ Stability boundary - MDBM -----------------


ax1 = Axis(0:2:10, "δ") # initial grid in x direction
ax2 = Axis(-1.5:1.4:1.5, "b") # initial grid in y direction
function fooDelay(δ, b)
    ABSmuMax = spectralradius(dpMathieu; p=(ζ, δ, ϵ, b, τ, T))
    return ABSmuMax - 1.0
end
mymdbm = MDBM_Problem(fooDelay, [ax1, ax2])
@time MDBM.solve!(mymdbm, 4)
#points where the function foo was evaluated
x_eval, y_eval = getevaluatedpoints(mymdbm)
x_sol, y_sol = getinterpolatedsolution(mymdbm)
#scatter(x_eval,y_eval,markersize=1)
#scatter!(x_sol,y_sol,markersize=2)
scatter!(x_sol, y_sol, markersize=1)


