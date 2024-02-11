5 + 5
using Revise
#includet("src\\DDE_mapping.jl")
includet("src\\DDE_mapping_types.jl")
includet("src\\DDE_mapping_functions.jl")

using BenchmarkTools
using Plots
plotly()
using Profile
using StaticArrays
using DifferentialEquations

# using Memoization

using MDBM


using Integrals

function robot_neutral(u, h, p, t)
    χ, r, γ, μ, κ1, κ2, τ = p
    u1 = u[1]
    u2 = u[2]
    du1 = u[3]
    du2 = u[4]

    ddu1_tau = h(p, t - τ, Val{1})[3]
    ddu2_tau = h(p, t - τ, Val{1})[4]


    #\\TODO: ez csak egy gyors test a distributed-re, de ettől nagyon lassú lesz, és még le kell tesztelni!
    #distributes
    h_foo(x,p)=h(p, x, Val{1})[3]
    prob = IntegralProblem(h_foo, t - τ*1.2, t - τ*0.8)
    solINT = solve(prob, HCubatureJL(); reltol = 1e-3, abstol = 1e-3)
    ddu1_DISTtau=solINT.u

    d_u1 = du1
    d_u2 = du2
    d_du1 = 0.3*ddu1_DISTtau+-2 * χ * r * γ * (du1 - du2) - γ^2 * r * (u1 - u2) - μ * r * (u1 - u2)^3 - u1 + κ1 * ddu1_tau + κ2 * ddu2_tau
    d_du2 = -2 * χ * γ * (du2 - du1) - γ^2 * (u2 - u1) - μ * (u2 - u1)^3

    SA[d_u1, d_u2, d_du1, d_du2]
end
Base.:+(a::SVector, b::Bool) = a .+ b
#function f2DelayMathieu(v, u, h, p, t)
#    ζ, δ, ϵ, b, τ , T = p
#    ddx = -(δ+ϵ*cos(t)) * u + (-2*ζ) * v + b * h(p, t - τ)[1] 
#end

χ = 0.1
r = 1.0
γ = 0.5
μ = 0.0
#κ1 = 0.46#-0.99#0.465
#κ2 = 0.0
κ1 = 0.0#-0.99#0.465
κ2 = 1.5
τ = 0.4
T = sqrt(2) / 6#pi/2#0.3
p = χ, r, γ, μ, κ1, κ2, τ

#p = (ζ, ωn, k, τ,10.0)

u0 = SA[1.0, 1.0, 1.0, 1.0]
h(p, t) = SA[1.0, 0.0, 0.0, 0.0]
h(p, t, deriv::Type{Val{1}}) = SA[0.0, 0.0, 1.0, 0.0]



probrobot_neutral = DDEProblem(robot_neutral, u0, h, (0.0, T * 1000.0), p; constant_lags=[τ], neutral=true)
#probMathieu =SecondOrderDDEProblem(f2DelayMathieu, 0.0,0.0,h, (0.0, T * 10.0), p; constant_lags=[τ])

sol = solve(probrobot_neutral, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)
# plot(sol(sol.t[end] .- (0.0:0.01:τ*1.0)))



#---------------------- Spectrum test --------------------
Nstep = 350
τmax = 2pi 
dpdp = dynamic_problemSampled(probrobot_neutral, MethodOfSteps(BS3()), τmax,
    T; Historyresolution=Nstep, eigN=20, zerofixpont=true, dt=0.001);

# fix point by affine map
@time mu = spectrum(dpdp; p=p);
plot(log.(abs.(mu[1])))
scatter((mu[1]))
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))

spectralradius(dpdp; p=p)

#---------------------- collocated - stability map --------------------
using MDBM

ax1 = Axis(0.1:0.2:1, "τ") # initial grid in x direction
ax2 = Axis(-1.5:0.2:1.5, "κ1") # initial grid in y direction
function fooMathieu(τ, κ1)
    τmax =(τ + 0.1)*1.2
    Nstep=20
    dploc = dynamic_problemSampled(probrobot_neutral, MethodOfSteps(BS3()), τmax,
    T; Historyresolution=Nstep, eigN=1, zerofixpont=true, dt=0.05)
    ABSmuMax = spectralradius(dploc; p=(χ, r, γ, μ, κ1, κ2, τ))
    return ABSmuMax - 1.0
end
mymdbm = MDBM_Problem(fooMathieu, [ax1, ax2])
iteration = 3#number of refinements (resolution doubling)
@time MDBM.solve!(mymdbm, iteration)
#points where the function foo was evaluated
x_eval, y_eval = getevaluatedpoints(mymdbm)
#interpolated points of the solution (approximately where foo(x,y) == 0 and c(x,y)>0)
x_sol, y_sol = getinterpolatedsolution(mymdbm)
#scatter(x_eval,y_eval,markersize=1)
#scatter!(x_sol,y_sol,markersize=3)
scatter(x_sol, y_sol, markersize=1)




println("----------Start brute-force---------------")
τv = 0.1:0.02:1.0 # initial grid in x direction
κ1v = -1.0:0.02:1.0 # initial grid in y direction
κ1v = -0.3:0.02:0.7 # initial grid in y direction
Spek = zeros(size(κ1v, 1), size(τv, 1))
#Threads.@threads
@time Threads.@threads for j in 1:size(τv, 1)
    τ = τv[j]
    τmax =(τ + 0.1)*1.2
    Nstep=20
    dploc = dynamic_problemSampled(probrobot_neutral, MethodOfSteps(BS3()), τmax,
    T; Historyresolution=Nstep, eigN=1, zerofixpont=true, dt=0.01)
    Threads.@threads for i in 1:size(κ1v, 1)
        κ1 = κ1v[i]
        muMAX = spectralradius(dploc; p=(χ, r, γ, μ, κ1, κ2, τ))
        Spek[i, j] = muMAX
    end
end

Spek_sat = deepcopy(Spek);
Spek_sat[Spek_sat.>1.0] .= 0.0;
heatmap(τv, κ1v, log.(Spek_sat))

scatter!(x_sol, y_sol, color=:blue,markersize=3)


#---------------------- non-collocated - stability map --------------------


using MDBM

ax1 = Axis(0.1:0.2:1.1, "τ") # initial grid in x direction
ax2 = Axis(-2.5:0.2:2.5, "κ2") # initial grid in y direction
function fooMathieu(τ, κ2)
    τmax = (τ + 0.1)*1.2
    Nstep=20
    dploc = dynamic_problemSampled(probrobot_neutral, MethodOfSteps(BS3()), τmax,
    T; Historyresolution=Nstep, eigN=1, zerofixpont=true, dt=0.05)
    ABSmuMax = spectralradius(dploc; p=(χ, r, γ, μ, κ1, κ2, τ))
    return ABSmuMax - 1.0
end
mymdbm = MDBM_Problem(fooMathieu, [ax1, ax2])
iteration = 4#number of refinements (resolution doubling)
@time MDBM.solve!(mymdbm, iteration)
#points where the function foo was evaluated
x_eval, y_eval = getevaluatedpoints(mymdbm)
#interpolated points of the solution (approximately where foo(x,y) == 0 and c(x,y)>0)
x_sol, y_sol = getinterpolatedsolution(mymdbm)
#scatter(x_eval,y_eval,markersize=1)
#scatter!(x_sol,y_sol,markersize=3)
scatter(x_sol, y_sol, markersize=1)




println("----------Start brute-force---------------")
τv = 0.1:0.01:1.1 # initial grid in x direction
κ2v = -2.5:0.025:2.5 # initial grid in y direction
Spek = zeros(size(κ2v, 1), size(τv, 1))
#Threads.@threads
@time Threads.@threads for j in 1:size(τv, 1)
    τ = τv[j]
    τmax = (τ + 0.1)*1.2
    Nstep=20
    dploc = dynamic_problemSampled(probrobot_neutral, MethodOfSteps(BS3()), τmax,
    T; Historyresolution=Nstep, eigN=1, zerofixpont=true, dt=0.01)
    Threads.@threads for i in 1:size(κ2v, 1)
        κ2= κ2v[i]
        muMAX = spectralradius(dploc; p=(χ, r, γ, μ, κ1, κ2, τ))
        Spek[i, j] = muMAX
    end
end

Spek_sat = deepcopy(Spek);
Spek_sat[Spek_sat.>1.0] .= 0.0;
heatmap(τv, κ2v, log.(Spek_sat))

scatter!(x_sol, y_sol, color=:blue,markersize=3)






