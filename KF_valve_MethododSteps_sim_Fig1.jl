#Kádár Fanni - Szelep rezgés - sim_only
# Algeb. Delay Diff. Eq.
5 + 5

using Revise
#using DDE_mapping

using Interpolations

using LinearAlgebra
using BenchmarkTools
using Plots
plotly()
using Profile
using StaticArrays
using DifferentialEquations

using LaTeXStrings
# using Memoization

using MDBM
#using MAT
using Peaks
using NaNStatistics

using ProgressLogging


function valve_KF_DADE_exp(y, h, p, t)
    ζ, δ, τ, β, Φ, q, h4 = p
    y4 = -h4(p, t - τ) - (Φ^2.0 * y[1]^2.0) / 2.0 + Φ * y[1] * abs(2 * h4(p, t - τ) + y[3] + (Φ^2.0 * y[1]^2.0) / 4.0)^0.5#Singular "f=y4"

    d_y1 = y[2]
    d_y2 = -2 * ζ * y[2] - (y[1] + δ) + y[3] - y4 + h4(p, t - τ)
    d_y3 = β * (q - 1 / Φ * (y4 + h4(p, t - τ)))
    SA[d_y1, d_y2, d_y3]
end


function y4_fun(sol, p)
    ζ, δ, τ, β, Φ, q, _ = deepcopy(sol.prob.p)
    y4series = [
        -sol.prob.p[end](p, t - τ) - (Φ^2.0 * sol[1, it]^2.0) / 2.0 +
        Φ * sol[1, it] * abs(2 * sol.prob.p[end](p, t - τ) + sol[3, it] + (Φ^2.0 * sol[1, it]^2.0) / 4.0)^0.5
        for (it, t) in enumerate(sol.t)
    ]
    #  return y4series
    nodes = (sol.t,)
    itp = extrapolate(interpolate(nodes, y4series, Gridded(Linear())), Flat())

    h4_interp(p, t) = itp(t)

    return h4_interp #, y4series
end
## - futtatás 

ζ = 0.39
δ = 3.0
τ = 0.25
β = 3.7
Φ = 48.2
#q = 3.0;
q = 10

#u0 = SA[0.1, 0.1, 0.1, 0.1,0.0]
#h(p, t) = SA[2.0,  0.0, 2.5,q*Φ/2,0.0]
h(p, t) = SA[2.0, 0.0, 2.5, q*Φ/2]

T = τ 

Nmutip = 1
Base.:+(a::SVector, b::Bool) = a .+ b

h123(p, t) = SA[2.0, 0.0, 2.5]
h4(p, t) = q * Φ / 2

p = (ζ, δ, τ, β, Φ, q, h4)
u0 = h123(p, 0.0)
valve_KF_DADE_exp_fun = DDEFunction(valve_KF_DADE_exp)


t_final = [0.0]
u_final = [u0]
f_final = [h4(p, 0.0)]

sizehint!(t_final, 1_000_000)
sizehint!(u_final, 1_000_000)
sizehint!(f_final, 1_000_000)
kt = 0
#kt=10000
#for kt in 0:10000
for kt in 1:20
    #kt += 1

    prob_valve_exp = DDEProblem(valve_KF_DADE_exp_fun, u0, h123, [0.0, T] .+ kt * τ, p, constant_lags=[τ])

    #@time 
    # sol_exp = solve(prob_valve_exp, MethodOfSteps(Rodas5()), reltol=1e-12, abstol=1e-12)
    
    resol = 1e-6
    sol_exp = solve(prob_valve_exp, MethodOfSteps(BS3()), reltol=resol, abstol=resol, dt=0.5)

    plot!(sol_exp, linestyle=:dash)

    h44 = y4_fun(sol_exp, p)
    fser = [h44(p, t) for t in sol_exp.t]
    
    PLOTFIG=plot!(sol_exp.t, fser, linestyle=:dash,xlim=(0,sol_exp.t[end]),linewidth = 2)
display(PLOTFIG)
    p = (ζ, δ, τ, β, Φ, q, h44)
    u0 = sol_exp[end]

    append!(t_final, sol_exp.t)
    append!(u_final, sol_exp.u)
    append!(f_final, fser)

    if false #mod(kt, 100) == 0
        println(kt)
        y1 = getindex.(u_final, 1)
        #plot(t_final,y1,linewidth = 2)
        MaxInd = argmaxima(y1)
        y1locmax = y1[MaxInd]
        pp = scatter(t_final[MaxInd], y1[MaxInd], color=:red)
        display(pp)
    end
end

plot!(t_final, getindex.(u_final,1),linecolor=:green,linestyle=:dash,linewidth = 2)
plot!(t_final, getindex.(u_final,2),linecolor=:green,linestyle=:dash,linewidth = 2)
plot!(t_final, getindex.(u_final,3),linecolor=:green,linestyle=:dash,linewidth = 2)
aa = plot!(t_final, f_final,linecolor=:green,linestyle=:dash,linewidth = 2)