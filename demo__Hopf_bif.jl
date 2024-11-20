#Demo: Delayed Nonline Oscill with nonlinearity
5 + 5

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

function DelayedNonlineOscill(u, h, p, t)
    # Parameters
    ζ, δ, b, τ, μ = p
    # Components of the delayed differential equation
    dx = u[2]
    ddx = -δ * u[1] - 2 * ζ * u[2] + b * h(p, t - τ)[1] - μ * u[2]^3
    # Update the derivative vector
    SA[dx, ddx]
end

Base.:+(a::SVector, b::Bool) = a .+ b
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?
ζ = 0.01         # damping coefficient
δ = 1.#0.2          # nat. freq
b = -0.06#0.3
τ = 0.5#2pi          # Time delay
μ = 5
p = ζ, δ, b, τ, μ
#p = (ζ, ωn, k, τ,10.0)

# test simulation ---------------
#initial condition
u0 = SA[0.001, 0.0]
#history function
h(p, t) = SA[0.0; 0.0]

Tlongsim = 2000.2
Tend = 7.0
prob_long = DDEProblem(DelayedNonlineOscill, u0, h, (0.0, Tlongsim), p; constant_lags=[τ])


#Parameters for the solver as a Dict (it is necessary to collect it for later use)
Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-7)#

sol = solve(remake(prob_long, p=(ζ, δ, b, τ, μ)); Solver_args...)#abstol,reltol
plot(sol)
##
#last period of the long simulation:
t_select_period = 0.0:0.01:Tend
t_select_delay = LinRange(0, τ, 200)
sol_period = sol(sol.t[end] .- t_select_period)
sol_delay = sol(sol.t[end] .- t_select_delay)

#plot the state phase of the last segment
plot(sol_period.t, sol_period[1, :])
plot!(sol_period.t, sol_period[2, :])
plot!(sol_delay.t, sol_delay[1, :])
plot!(sol_delay.t, sol_delay[2, :])
#plot the phase space (u - du)
plot(sol_delay[1, :], sol_delay[2, :])
plot!(sol_period[1, :], sol_period[2, :])


## ---------------- simulation max amplitude ----------------------
## parameters 
norm_solperiod = zeros(Float64, 0)
plot()
bv = LinRange(-0.1, 0.05, 50)
@time for bloc in bv
    println(bloc)
    #prob_long = DDEProblem(DelayedNonlineOscill, u0, h, (0.0, Tlongsim), (ζ, δ, bloc, τ, μ); constant_lags=[τ])
    #sol = solve(prob_long; Solver_args...)#abstol,reltol

    sol = solve(remake(prob_long, p=(ζ, δ, bloc, τ, μ)); Solver_args...)#abstol,reltol

    plot(sol)
    #last period of the long simulation:
    t_select_period = 0.0:0.01:Tend
    t_select_delay = eriod = 0.0:0.01:τ
    sol_period = sol(sol.t[end] .- t_select_period)
    sol_delay = sol(sol.t[end] .- t_select_delay)
    #push!(norm_solperiod, norm(sol_period.u))
    push!(norm_solperiod, maximum(getindex.(sol_period.u,1)))
end

plot(bv, norm_solperiod)
##





## ---------------- Affine mapping ---------------------------
using KrylovKit
Neig = 10#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-9, krylovdim=3 + 50, verbosity=0));

τmax = τ #maximal timedelay in the mapping
Nstep = 100 # discretization number of the mapping
Timeperiod = Tend # timeperiod of the mapping

#Creating the problem
dp_0_Tfix = dynamic_problemSampled(prob_long, Solver_args, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=true, affineinteration=0,
    Krylov_arg=Krylov_arg)

μ₀ = []
bv_affine = LinRange(-0.1, 0.05, 250)
@time for bloc in bv_affine
    mu, saff, sol0 = affine(dp_0_Tfix; p=(ζ, δ, bloc, τ, μ))
    push!(μ₀, mu[1])
end

plot(bv, norm_solperiod)
#scatter()
for k in 1:Neig
    plot!(bv_affine, (abs.(getindex.(μ₀, k)) .- 1) .* 10)
end
plot!()


## ---------------------- Hopf solution -------------------
b = -0.06




using ForwardDiff
u01_Dual = ForwardDiff.Dual{Float64}(0.05, 2.0) 
u01_Dual.partials
u01_Dual.value

#initial condition
u0 = SA[0.05, 0.0]
u0 = SA[ ForwardDiff.Dual{Float64}(0.05, 2.0) ,  ForwardDiff.Dual{Float64}(0.05, 2.0) ]
#history function
h(p, t) = SA[0.5; 0.0]
h(p, t) = SA[ ForwardDiff.Dual{Float64}(0.05, 2.0) ,  ForwardDiff.Dual{Float64}(0.05, 2.0) ]

TlongsimDUA=ForwardDiff.Dual{Float64}(Tlongsim, 0.0) 
prob_DUAL = DDEProblem(DelayedNonlineOscill, u0, h, (ForwardDiff.Dual{Float64}(0.0, 0.0) , TlongsimDUA), p; constant_lags=[τ])


function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
    u[2]
end
function affect_long!(integrator)
    #println(typeof(integrator))
    if integrator.t > 200.0#1500.0
        terminate!(integrator)
    end
end
function affect_short!(integrator)
    if integrator.t > 5.0
        terminate!(integrator)
    end
end
using DifferentialEquations
cb_long = ContinuousCallback(condition, affect_long!)
cb_short = ContinuousCallback(condition, affect_short!)

Solver_args_T_long = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-6,
    :callback => cb_long)#
Solver_args_T_short = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-6,
    :callback => cb_short)#

Nstep = 500
sol = solve(remake(prob_DUAL, p=(ζ, δ, b, τ, μ)); Solver_args_T_long...)#abstol,reltol
sol = solve(remake(prob_long, p=(ζ, δ, b, τ, μ)); Solver_args_T_long...)#abstol,reltol
plot(sol)



t_select_period = LinRange(-Tend, 0, 500)
sol_period = sol(sol.t[end] .+ t_select_period)

t_select_delay = LinRange(-τ, 0, Nstep)
sol_delay = sol(sol.t[end] .+ t_select_delay)

plot(t_select_period, sol_period[1, :])
plot!(t_select_period, sol_period[2, :])
plot!(t_select_delay, sol_delay[1, :],lw=3)
plot!(t_select_delay, sol_delay[2, :],lw=3)


plot(sol[1, :],sol[2, :])
plot!(sol_period[1, :],sol_period[2, :],lw=4)
plot!(sol_delay[1, :],sol_delay[2, :],lw=4)


#Creating the problem
Timeperiod = 20.0;#Tend
dp_Hopf_callback = dynamic_problemSampled(prob_long, Solver_args_T_short, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=false, affineinteration=2,
    Krylov_arg=Krylov_arg)
#@time mu, saff, sol0 = affine(dp_Hopf_callback, sol_delay.u; p=(ζ, δ, b, τ, μ))

@time mu, saff, sol0 = affine(dp_Hopf_callback, sol_delay.u; p=(ζ, δ, b, τ, μ))

plot(sol0)
plot()
# Comparing the solutions:
plot(sol_period[1, :], sol_period[2, :], marker=:circle, markersize=6, lab="")
plot!(getindex.(saff, 1), getindex.(saff, 2), marker=:circle, markersize=4, lab="")#
plot!(sol0[1, :], sol0[2, :], marker=:circle, lw=0, lab="")#marker=:cross,markersize=2)#
##


# Plotting the Floquet multipliers
#in log scale
plot(log.(abs.(mu[1])))
#in complex plane
scatter(mu[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))


using ForwardDiff
a = ForwardDiff.Dual{Float64}(00.0, 2.0) 
a.partials
a.value

a==-0.0


#TODO: ez nem hiszem, hogy jó megoldás
Base.:convert(::Type{Float64}, x::ForwardDiff.Dual{Float64, Float64, 1}) = x.value










dp=dp_Hopf_callback
s0 = sol_delay.u
s0=a0
plot()
plot!(dp.StateSmaplingTime, getindex.(s0, 1),lw=1)
plot!(dp.StateSmaplingTime, getindex.(s0, 2),lw=1)

v0, sol_loc = DDE_mapping.LinMap(dp, s0; p=(ζ, δ, b, τ, μ));
plot!(sol_loc.t,sol_loc[1,:])
plot!(sol_loc.t,sol_loc[2,:])
Ti = sol_loc.t[end]
#plot!(xlim=(-τ, Ti))

plot!(dp.StateSmaplingTime.+Ti, getindex.(v0, 1),lw=2)
plot!(dp.StateSmaplingTime.+Ti, getindex.(v0, 2),lw=2)
plot!(dp.StateSmaplingTime, getindex.(v0, 1),lw=2)
plot!(dp.StateSmaplingTime, getindex.(v0, 2),lw=2)
@show norm(s0-v0)
s0=v0
plot!(legend=false)



#println(norm(s0-v0))
Nstep = size(dp.StateSmaplingTime, 1)
s_start = rand(typeof(dp.Problem.u0), Nstep)



dp=dp_Hopf_callback
s0 = sol_delay.u
s0=a0
v0 = LinMap(dp, s0; p=p)[1]

#println("Float perturbation")
#EPSI_TODO_REMOVE = 1e-6;
#s_start .*= EPSI_TODO_REMOVE
#TheMapping(s) = (LinMap(dp, s + s0; p=p)[1] - v0)

println("Dual perturbation - it seems to be faster! ;-)")
one_espilon_Dual = ForwardDiff.Dual{Float64}(0.0, 1.0) 
TheMapping(s) = DDE_mapping.partialpart.(LinMap(dp, s * one_espilon_Dual + s0; p=p)[1] - v0)

vstart=TheMapping(s_start)


#dp=dp_Hopf_callback
#s0 = s_start
plot()
plot!(dp.StateSmaplingTime, getindex.(s0, 1),lw=1)
plot!(dp.StateSmaplingTime, getindex.(s0, 2),lw=1)

v0 = TheMapping(s0)

plot!(dp.StateSmaplingTime.+Ti, getindex.(v0, 1),lw=4)
plot!(dp.StateSmaplingTime.+Ti, getindex.(v0, 2),lw=4)
plot!(dp.StateSmaplingTime, getindex.(v0, 1),lw=5)
plot!(dp.StateSmaplingTime, getindex.(v0, 2),lw=5)
@show norm(s0-v0)
s0=v0
plot!(legend=false)

##

mus = getindex(schursolve(TheMapping, s_start, dp.Krylov_arg...), [3, 2, 1])



plot(log.(abs.(mus[1])))
#in complex plane
scatter(mus[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))

a0 = real.(DDE_mapping.find_fix_pont(s0, v0, mus[1], mus[2]))

plot(sol_delay.t .- sol.t[end], getindex.(a0, 1))
plot!(sol_delay.t .- sol.t[end], getindex.(a0, 2))

s0=a0




































































# ----------------------- creating stability chart -------------------------

Krylov_arg = (1, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=8 + 5, verbosity=0));

dpMathieu = dynamic_problemSampled(probMathieu, Solver_args, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=false, affineinteration=1,
    Krylov_arg=Krylov_arg)

## ------------------ Peak-2-Peak ampNorm color, and spectral radius - BRUTE FORCE in the plane of δ-ϵ  ----------------

b = 0.05 # delay parameter, Note: b=0 provide the traditional stability chart of the Mathieu equations
δv = -1.0:0.051:5.0 # initial grid in x direction
ϵv = -0.001:0.05:10.0 # initial grid in y direction
Aaff = zeros(size(ϵv, 1), size(δv, 1))
Spek_aff = zeros(size(ϵv, 1), size(δv, 1))

@time Threads.@threads for j in 1:size(δv, 1)
    @inbounds δ = δv[j]
    Threads.@threads for i in 1:size(ϵv, 1)
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

b = 0.05;
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

Krylov_arg = (1, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=8 + 5, verbosity=0));

dpMathieu = dynamic_problemSampled(probMathieu, Solver_args, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=false, affineinteration=1,
    Krylov_arg=Krylov_arg)


δv = 0:0.051:10 # initial grid in x direction
bv = -1.501:0.05:1.5 # initial grid in y direction
Aaff = zeros(size(bv, 1), size(δv, 1))
Spek_aff = zeros(size(bv, 1), size(δv, 1))

@time Threads.@threads for j in 1:size(δv, 1)
    @inbounds δ = δv[j]
    Threads.@threads for i in 1:size(bv, 1)
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


