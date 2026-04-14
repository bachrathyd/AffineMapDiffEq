#Demo: Duffing Chain
# from: Fabia Bayer
# DUFFING_CHAIN Dynamics of a forced Duffing chain 
# given by a series of n = length(x)/2 mass-spring-damper elements 
# with localized nonlinear stiffness and periodic forcing of the last element.
# 
# Fabia Bayer, INM, Uni Stuttgart, July 2024
#
# PARAMETERS
# ----------
# t : time (appears in the forcing F(t) = gamma*cos(omega*t) ). 
#     May be a row vector, then dfdx is a 3D array.
# x : state x = [q; qdot] where q(k) denotes the absolute coordinate of the k-th mass.
#     Must have as many columns as t.
#
# k1: linear spring constant
# k3: cubic spring constant
# l0: initial elongation of springs
# gamma, omega: Amplitude and frequency of the forcing
#
# idx_nonlin: Vector of indices indicating which connections have cubic
# spring stiffness (e.g. [1, 2, 4, 6])

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

using MAT

# Governing equation

function duffing_chain(x::T, h, p, t)::T where T
    k1, k3, d, l0, idx_nonlin, gamma, omega = p
    n = size(x, 1) ÷ 2
    q = x[1:n]
    qdot = x[n+1:end]

    # Relative coordinates: distance between masses
    q_rel = q .- vcat(zeros(1, length(t)), q[1:end-1, :]) .- l0
    qdot_rel = qdot .- vcat(zeros(1, length(t)), qdot[1:end-1, :])

    #dqrel_dq = toeplitz([1, -1, zeros(1, n-2)], [1, zeros(1, n-1)])

    # Interaction forces. F(k) is between masses k-1 and k. Pulling force is positive.
    F = zeros(typeof(x[1]), n + 1, length(t))
    F[end, :] .= gamma .* cos.(omega .* t)
    F[1:n, :] .= d .* qdot_rel .+ k1 .* q_rel
    F[idx_nonlin, :] .= F[idx_nonlin, :] .+ k3 .* q_rel[idx_nonlin, :] .^ 3

    #dF_dqrel = [k1 * I(n) for _ in 1:length(t)]
    #index_diags = [i:(n+1):n^2 for i in 1:(length(idx_nonlin))]
    #for (i, idx) in enumerate(idx_nonlin)
    #    for j in 1:length(t)
    #        dF_dqrel[j][idx, idx] += 3 * k3 * q_rel[idx, j]^2
    #    end
    #end
    #
    #dFdq = [dF_dqrel[j] * dqrel_dq for j in 1:length(t)]
    #dFdq = vcat(dFdq, [zeros(1, n) for _ in 1:length(t)]) # F_{n+1} = gamma*cos(omega*t) does not depend on q
    #
    #dFdqdot = [d * dqrel_dq for _ in 1:length(t)]
    #dFdqdot = vcat(dFdqdot, [zeros(1, n) for _ in 1:length(t)])

    # Assemble xdot
    xdot = SA[qdot..., (F[2:n+1] .- F[1:n])...]

    # Assemble derivative
    #dfdx = hcat([zeros(n, n) repmat(I(n), 1, length(t)); dFdq[2:n+1, :] .- dFdq[1:n, :], dFdqdot[2:n+1, :] .- dFdqdot[1:n, :]])

    return xdot #, dfdx
end


Base.:+(a::SVector, b::Bool) = a .+ b
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?


## parameters 
m = 1;
k1 = 0.5;
k3 = 1;
L0 = 1;
d = 0.2;


gamma = 0.3;

#n = 10;
#omega=0.364068217896367

n = 50;
omega = 0.115953733056267

idx_nonlin = 1:n; # all springs are nonlinear
T = 2 * pi / omega
τ = T / 10


p = k1, k3, d, L0, idx_nonlin, gamma, omega

u0 = SA[1.3 * L0 * (1:n)...; 0.0 * (1:n)...]

using ForwardDiff
one_espilon_Dual = ForwardDiff.Dual{Float64}(0.0, 1.0)
xx = u0 * one_espilon_Dual

h(p, t) = u0

Tsolve = T * 400.0
probDC = DDEProblem(duffing_chain, u0, h, (0.0, Tsolve), p; constant_lags=[τ])
t_select_period = 0.0:0.01:T
#Parameters for the solver as a Dict (it is necessary to collect it for later use)
#Solver_args = Dict(:alg => MethodOfSteps(BS3()),  :adaptive => false, :dt => 0.2, :verbose => false, :reltol => 1e-5)#, save_everystep=false)#abstol,reltol)
Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-5, :abstol => 1e-5, :saveat => Tsolve .- t_select_period)#
#Solver_args = Dict()#
duffing_chain(xx, h, p, 0.0)
duffing_chain(u0, h, p, 0.0)

@time sol = solve(probDC; Solver_args...);#abstol,reltol
#plot(sol)
@profview solve(probDC; Solver_args...);#abstol,reltol

#last period of the long simulation:
t_select_period = 0.0:0.01:T
t_select_delay = eriod = 0.0:0.01:τ
sol_period = sol(Tsolve .- t_select_period)
sol_delay = sol(Tsolve .- t_select_delay)

#plot the state phase of the last segment
plot(sol_period.t, getindex.(sol_period.u, 1))
plot!(sol_period.t, getindex.(sol_period.u, 1 + n))
plot!(sol_delay.t, getindex.(sol_delay.u, 1))
plot!(sol_delay.t, getindex.(sol_delay.u, 1 + n))
#plot the phase space (u - du)
plot(getindex.(sol_delay.u, 1), getindex.(sol_delay.u, 1 + n))
plot!(getindex.(sol_period.u, 1), getindex.(sol_period.u, 1 + n))


# ---------------- Affine mapping ---------------------------

Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-5, :abstol => 1e-5)#
using KrylovKit
Neig = 4#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-5, krylovdim=Neig + 10, verbosity=0));

τmax = τ #maximal timedelay in the mapping
Nstep = 10 # discretization number of the mapping
Timeperiod = T # timeperiod of the mapping

#Creating the problem
dpDC = dynamic_problemSampled(probDC, Solver_args, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=false, affineinteration=0,
    Krylov_arg=Krylov_arg)


s0 = [u0 for _ in 1:Nstep]
#solvig the problem (mu: Floquet multiplier, saff: discretized foxed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
@time mu, saff, sol0 = affine(dpDC, s0; p=p);


#using Profile
#@profview affine(dpDC, s0; p=p);


plot(sol0.t, getindex.(sol0.u, 1))
plot!(dpDC.StateSmaplingTime .+ sol0.t[end], getindex.(saff, 1))
# Comparing the solutions:
plot(sol_period[1, :], sol_period[1+n, :], lw=3,lab="long sim Julia")#, marker=:circle,markersize=2
#plot!(getindex.(saff,1), getindex.(saff,1+n), marker=:circle,markersize=4,lab="")#
plot!(sol0[1, :], sol0[1+n, :], marker=:circle, lw=1, lab="periodic orbit Julia")#marker=:cross,markersize=2)#
##



# -------------- MAT result read ---------------
using MAT
vars = matread("C:\\Users\\Bacharthy\\Downloads\\MWE_duffing_chain\\xs.mat")
varnames = keys(vars)
xs = vars["xs"]
plot!(xs[1, :], xs[1+n, :], marker=:square, markersize=3, lab="Matlab solution")

#----------------------- Eigen values -----------------

# Plotting the Floquet multipliers
#in log scale
plot(log.(abs.(mu[1])))
#in complex plane
scatter(mu[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))




plot!(dpDC.StateSmaplingTime, getindex.(mu[2][1], 1))
plot(mu[2][3][1])



















# ----------------------- creating stability chart -------------------------

Krylov_arg = (1, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=8 + 5, verbosity=0));

dpDC = dynamic_problemSampled(probMathieu, Solver_args, τmax,
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
        muaff, s0aff, sol0 = affine(dpDC; p=(ζ, δ, ϵ, b, τ, T))
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
    ABSmuMax = spectralradius(dpDC; p=(ζ, δ, ϵ, b, τ, T))
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

dpDC = dynamic_problemSampled(probMathieu, Solver_args, τmax,
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
        muaff, s0aff, sol0 = affine(dpDC; p=(ζ, δ, ϵ, b, τ, T))
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
    ABSmuMax = spectralradius(dpDC; p=(ζ, δ, ϵ, b, τ, T))
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


