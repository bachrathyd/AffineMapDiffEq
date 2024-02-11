#Bartfai András - ENOC cikk
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
#Bártfai - Enoc test
function BA_Dist_neutral(u, h, p, t)
    ωn, τf, τr, w, μ = p

    u1 = u[1]
    du1 = u[2]
    τ = τf + τr

    ddu1_tau = h(p, t - τ, Val{1})[2]


    #\\TODO: ez csak egy gyors test a distributed-re, de ettől nagyon lassú lesz, és még le kell tesztelni!
    #distributes
    IntInner(x, p) = w(x) * h(p, t + x, Val{1})[2]
    prob = IntegralProblem(IntInner, -τ, -τr)
    ur = solve(prob, HCubatureJL(); reltol=1e-5, abstol=1e-5).u

    #    -0.001*du1
    d_u1 = du1
    d_du1 = -ωn * ωn * u1 + ωn * ωn * ur - (μ * (u1^3 - ur^3) + 3.0 * μ * (u1^2 * ur - u1 * ur^2))


    SA[d_u1, d_du1]
end
Base.:+(a::SVector, b::Bool) = a .+ b

ωn = 6.0 * 2pi
τf = 0.05;
τr = 10e-3;
K = -0.001;
0.08;
w(t) = K;
μ = 0.00;#0.1


p = (ωn, τf, τr, w, μ)
T = sqrt(2) / 126#pi/2#0.3


u0 = SA[1.0, 1.0]
h(p, t) = SA[1.0, 0.0]
h(p, t, deriv::Type{Val{1}}) = SA[-10.0, 1.0]


pro_BD_Enoc = DDEProblem(BA_Dist_neutral, u0, h, (0.0, T * 100.0), p, neutral=true)#; constant_lags=[τ]

sol = solve(pro_BD_Enoc, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)
# plot(sol(sol.t[end] .- (0.0:0.01:τ*1.0)))



#---------------------- Spectrum test --------------------
Nstep = 40
τmax = 0.1
dpdp = dynamic_problemSampled(pro_BD_Enoc, MethodOfSteps(BS3()), τmax,
    T; Historyresolution=Nstep, eigN=20, zerofixpont=true, dt=0.005);

# fix point by affine map
@time mu, saff = affine(dpdp; p=p);
plot(log.(abs.(mu[1])))
scatter((mu[1]))
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))


#---------------------- collocated - stability map --------------------
using MDBM

#ax1 = Axis(-0.009:0.01:0.0005, "τf") # initial grid in x direction
#ax2 = Axis(-0.5:0.25:0.5, "K") # initial grid in y direction
ax1 = Axis(LinRange(-0.008,0.015,5), "τf") # initial grid in x direction
ax2 = Axis(LinRange(-0.5,0.504,6), "K") # initial grid in y direction
function fooMathieu(τf, K)
    τmax = maximum([τr,τr + τf]) * 1.2
    Nstep = 500
    dploc = dynamic_problemSampled(pro_BD_Enoc, MethodOfSteps(BS3()), τmax,
        T; Historyresolution=Nstep, eigN=2, zerofixpont=true, dt=0.001)
    w(θ) = K
    #println((τf, K))
    mu, saff = affine(dploc; p=(ωn, τf, τr, w, μ))
    
    ABSmuMax = abs(mu[1][1]) ;#abs(mu[1][3])
    return ABSmuMax - 1.0
    
    #MuMin1_prod = prod(abs.(mu[1][1:4]) .-1 ) ;#abs(mu[1][3])
    #return MuMin1_prod#ABSmuMax - 1.0
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
scatter(x_sol, y_sol, markersize=2,xlabel="τf", ylabel="K")




println("----------Start brute-force---------------")
τfv = LinRange(-0.0095,0.015,40)#-0.005:0.01:0.1 # initial grid in x direction
#τfv = LinRange(0.0001,0.02,30)#-0.005:0.01:0.1 # initial grid in x direction
Kv = LinRange(-0.5,0.504,40)#-1.0:0.25:1.0 # initial grid in y direction
Spek = zeros(size(Kv, 1), size(τfv, 1))
#Threads.@threads
@time Threads.@threads for j in 1:size(τfv, 1)
    println(j)
    τf = τfv[j]
    τmax = maximum([τr,τr + τf]) * 1.05
    Nstep = 500
    dploc = dynamic_problemSampled(pro_BD_Enoc, MethodOfSteps(BS3()), τmax,
        T; Historyresolution=Nstep, eigN=3, zerofixpont=true, dt=0.001)
    #@time 
    Threads.@threads for i in 1:size(Kv, 1)
        K = Kv[i]
        w(θ) = K
        mu, saff = affine(dploc; p=(ωn, τf, τr, w, μ))
        muMAX = abs(mu[1][1])
        #muMAX = spectralradius(dploc; p=(ωn, τf, τr, w, μ))
        Spek[i, j] = muMAX
    end
end

Spek_sat = deepcopy(Spek);
Spek_sat[Spek_sat.>1.0] .= 1.0;
heatmap(τfv, Kv, (Spek_sat),xlabel="τf", ylabel="K")

scatter!(x_sol, y_sol, color=:blue, markersize=2)




