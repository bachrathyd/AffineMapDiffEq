#Demo: neutral delay systems, tau1-tau2 space
5+5
using DDE_mapping

using Plots
plotly()
using StaticArrays
using DifferentialEquations

using LinearAlgebra

using MDBM


# Governing equation

function mixedDelay_oscill(u, h, p, t)
    # Parameters
    a, b, c,zeta, τ1, τ2 = p
    # Components of the delayed differential equation
    dx = u[2]    
    ddx = -a * u[1] -zeta * u[2] + b * h(p, t - τ1)[2] +
       c * h(p, t - τ2,Val{1})[2] 
    # Update the derivative vector
    @MArray[dx, ddx]
    #SA[dx, ddx]
end
Base.:+(a::SVector, b::Bool) = a .+ b
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?


## parameters 
a = 1.0#
zeta=0.2
b = 0.4#0.7
c = 0.3#0.5
τ1=4.0
τ2=5.0
τmax=maximum([τ1,τ2])
T=τmax

p= (a, b, c, zeta, τ1, τ2) 
# test simulation ---------------
#initial condition
u0 = @MArray [1.5, 0.5]
#u0 = SA[1.0, 0.0]
#history function
h(p, t) = @MArray [0.8; -1.0]
h(p, t, ::Type{Val{1}}) = @MArray [0.5; 1.5]
#h(p, t, ::Type{Val{2}}) = @MArray [0.0; 1.0]
    
#h(p, t) = SA[0.0; 0.0]
probMixtau = DDEProblem(mixedDelay_oscill, u0, h, (0.0, T * 55), p) #, neutral=true)#; constant_lags=[τ]

#Parameters for the solver as a Dict (it is necessary to collect it for later use)
#Solver_args = Dict(:alg => MethodOfSteps(BS3()), :callback => cb, :adaptive => true, :dt => 0.01, :verbose => false, :reltol => 1e-9)#, sτ1ve_everystep=false)#abstol,reltol)
#Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-12)#
Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-15, :dtmax => 1e-1)#; constrained = true
sol = solve(probMixtau; Solver_args...)#abstol,reltol

plot(sol)
plot!(sol(sol.t, Val{1}, idxs=2),lw=2)

#last period of the long simulation:
t_select_period = 0.0:0.01:T
t_select_delay = eriod = 0.0:0.001:τmax
sol_period = sol(sol.t[end] .- t_select_period)
sol_delay = sol(sol.t[end] .- t_select_delay)

#plot the state phase of the last segment
plot(sol_period.t, getindex.(sol_period.u, 1))
plot!(sol_period.t, getindex.(sol_period.u, 2))
plot!(sol_delay.t, getindex.(sol_delay.u, 1))
plot!(sol_delay.t, getindex.(sol_delay.u, 2))
#plot the phase space (u - du)
plot(getindex.(sol_delay.u, 1), getindex.(sol_delay.u, 2))
plot!(getindex.(sol_period.u, 1), getindex.(sol_period.u, 2))


# ---------------- Affine mapping ---------------------------
using KrylovKit
Neig = 8#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-32, krylovdim=18 + 5, verbosity=0));

Nstep = 125 # discretization number of the mapping
Timeperiod =5*T# dh*2#copy(dh) # timeperiod of the mapping

#Creating the problem
probMixtau = DDEProblem(mixedDelay_oscill, u0, h, (0.0, Timeperiod), p)#; constant_lags=[τ]

dp_mix_delay = dynamic_problemSampled(probMixtau, Solver_args, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=true, affineinteration=0,
    Krylov_arg=Krylov_arg)


#solvig the problem (mu: Floquet multiplier, saff: discretized foxed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
@time mu, saff, sol0 = affine(dp_mix_delay; p=(a, b, c,zeta, τ1, τ2) );


# Comparing the solutions:
plot(sol_period[1, :], sol_period[2, :], marker=:circle, markersize=6, lab="")
plot!(getindex.(saff, 1), getindex.(saff, 2), marker=:circle, markersize=4, lab="")#
plot!(sol0[1, :], sol0[2, :], marker=:circle, lw=0, lab="")#marker=:cross,markersize=2)#
#


# Plotting the Floquet multipliers
#in log scale
plot(log.(abs.(mu[1])))
#in complex plane
scatter(mu[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))
## ----------------------- creating stability chart -------------------------

Krylov_arg = (1, :LM, KrylovKit.Arnoldi(tol=1e-8, krylovdim=7, verbosity=0));

τmax=32
Nstep = 100
Timeperiod=5*τmax
dp_mix_delay = dynamic_problemSampled(probMixtau, Solver_args, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=true, affineinteration=0,
    Krylov_arg=Krylov_arg)

# ------------------ Peak-2-Peak ampNorm color, and spectral radius - BRUTE FORCE in the plane of δ-ϵ  ----------------
τ1v = 0.05:0.2:15.0 # initial grid in x direction
τ2v = 0.05:0.2:15.0 # initial grid in y direction
Aaff = zeros(size(τ1v, 1), size(τ2v, 1))
Spek_aff = zeros(size(τ1v, 1), size(τ2v, 1))

#Threads.@threads
@time Threads.@threads for j in 1:size(τ2v, 1)
    @inbounds τ2 = τ2v[j]
    @show j
    Threads.@threads for i in 1:size(τ1v, 1)
        @inbounds τ1 = τ1v[i]
        muaff, s0aff, sol0 = affine(dp_mix_delay; p=(a, b, c,zeta, τ1, τ2))
        Aaff[i, j] = norm(getindex.(s0aff, 1)) # norm of the motion
        # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1))) #maximum amplitude of the orbit
        Spek_aff[i, j] = maximum(abs.(muaff[1]))
    end
end


#Plotting the maximal amplitud on the stable domain only
Spek_affsat = deepcopy(Spek_aff);
Spek_affsat[Spek_affsat.>1.0] .= 1.0;
#heatmap(τ2v, τ1v, log.(Spek_affsat))
heatmap(τ1v, τ2v, log.(Spek_affsat'))


#------------------ Stability boundary - MDBM -----------------

ax1 = Axis(0.0:3:15.0, "τ1") # initial grid in x direction
ax2 = Axis(0.0:3:15.0, "τ2") # initial grid in y direction
function fooDelay(τ1, τ2)
    ABSmuMax = spectralradius(dp_mix_delay; p=(a, b, c,zeta, τ1, τ2))
    #return ABSmuMax - 1.0
    return log(ABSmuMax)
end
mymdbm = MDBM_Problem(fooDelay, [ax1, ax2])
@time MDBM.solve!(mymdbm, 4)
x_sol, y_sol = getinterpolatedsolution(mymdbm)
scatter!(x_sol, y_sol, markersize=1)

##--------------------------