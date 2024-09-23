#KF_valve_time_Fig1

#Kádár Fanni - Szelep rezgés
# Algeb. Delay Diff. Eq.
5 + 5

#using Revise
#using DDE_mapping

using BenchmarkTools
using Plots
#plotly()
gr()
using Profile
using StaticArrays
using DifferentialEquations
using LinearAlgebra
using LaTeXStrings
# using Memoization

using ProgressLogging
using MDBM

function valve_KF_DADE(y, h, p, t)
    ζ, δ, τ, β, Φ, q = p

    d_y1 = y[2]
    d_y2 = -2 * ζ * y[2] - (y[1] + δ) + y[3] - y[4] + h(p, t - τ)[4]
    d_y3 = β * (q - 1 / Φ * (y[4] + h(p, t - τ)[4]))
    CORESRQ = (y[3] - y[4] + h(p, t - τ)[4])
    #CORESRQ=maximum([0,CORESRQ])
    # d_X= Φ* y[1]*  abs( CORESRQ) ^0.5 - h(p, t - τ)[4]-y[4]#Singular "f=y4"
    d_X = -y[4] - h(p, t - τ)[4] - (Φ^2.0 * y[1]^2.0) / 2.0 + Φ * y[1] * abs(2 * h(p, t - τ)[4] + y[3] + (Φ^2.0 * y[1]^2.0) / 4.0)^0.5#Singular "f=y4"
    d_y5 = d_X
    #d_y6=CORESRQ
    #    SA[d_y1,d_y2,d_y3,d_X,d_y5]
    SA[d_y1, d_y2, d_y3, d_X]
end

Base.:+(a::SVector, b::Bool) = a .+ b
#ζ = 0.25;
#δ = 3.0
#τ = pi / 3.0;
#β = 10.0
#Φ = 48.2
#q = 2.0;
#4.0;
#6.0;
#q = 6.0


ζ = 0.39
δ = 3.0
τ = 0.25
β = 3.7
Φ = 48.2
#q = 3.0;
q = 10

p = (ζ, δ, τ, β, Φ, q)

#u0 = SA[0.1, 0.1, 0.1, 0.1,0.0]
#h(p, t) = SA[2.0,  0.0, 2.5,q*Φ/2,0.0]
h(p, t) = SA[2.0, 0.0, 2.5, q*Φ/2]
@show u0 = h(p, 0.0)#  SA[0.1, 0.1, 0.1, 0.1]



function calc_powi(powi)
    println(powi)
    println(10.0^powi)

    T = τ * 2000
    if powi < -10
        T = τ * 20
    end
    #M=diagm([1.0,1.0,1.0,1.0])
    #M=diagm([1.0,1.0,1.0,0.01])
    #M=diagm([1.0,1.0,1.0,0.0001])
    #M=diagm([1.0,1.0,1.0,0.000001])
    M = diagm([1.0, 1.0, 1.0, 10.0^powi])
    #M=diagm([1.0,1.0,1.0,0.0]) #1e-9

    valve_KF_DADE_fun = DDEFunction(valve_KF_DADE, mass_matrix=M)
    prob_valve = DDEProblem(valve_KF_DADE_fun, u0, h, (0.0, T), p, constant_lags=[τ],
        progress=true, maxiters=1e6)
    #ROS3P 
    resol = 1e-6
    @time sol = solve(prob_valve, MethodOfSteps(Rodas5()),
        reltol=resol, abstol=resol)#abstol,reltol
    #plot!(sol)

    xx = sol.t
    y1 = getindex.(sol.u, 1)
    y2 = getindex.(sol.u, 2)
    y3 = getindex.(sol.u, 3)
    f = getindex.(sol.u, 4)
    return xx, y1, y2, y3, f
end

ppal=cgrad([:purple, :green])
using Colors
powi=-2.0
fp(x)=RGB((4+x)/2,0.2,1-(4+x)/2)
xx, y1, y2, y3, f = calc_powi(powi)
plot(xx, y1,color_palette=ppal,linewidth = 2,color = fp(powi))
plot!(xx, y2,color_palette=ppal,linewidth = 2,color = fp(powi))
plot!(xx, y3,color_palette=ppal,linewidth = 2,color = fp(powi))
plot!(xx, f,color_palette=ppal,linewidth = 2,color = fp(powi))



#for powi = [collect(-2.5:-0.5:-4)...,-Inf]
for powi = [collect(-2.:-0.5:-4.0)...]
    #for powi = [collect(-2.5:-0.5:-4)...]
    @show fp(powi)
    xx, y1, y2, y3, f = calc_powi(powi)
    plot!(xx, y1,color_palette=ppal,linewidth = 2,color = fp(powi))
    plot!(xx, y2,color_palette=ppal,linewidth = 2,color = fp(powi))
    plot!(xx, y3,color_palette=ppal,linewidth = 2,color = fp(powi))
    aa = plot!(xx, f,color_palette=ppal,linewidth = 2,color = fp(powi))
    display(aa)
end

xx, y1, y2, y3, f = calc_powi(-Inf)
plot!(xx, y1, linecolor=:magenta,linewidth = 2)
plot!(xx, y2, linecolor=:magenta,linewidth = 2,)
plot!(xx, y3, linecolor=:magenta,linewidth = 2)
aa = plot!(xx, f, linecolor=:magenta,linewidth = 2)



#TODO: not working run-by-hand
include("KF_valve_MethododSteps_sim_Fig1.jl")






display(aa)
aa = plot!(xlims=(00, 500), ylims=(-10, 260), xlabel="time", ylabel="y", dpi=1000, legend=false),
display(aa)
savefig("KF_valve_sim_Fig_all.png")
aa = plot!(xlims=(450, 500), ylims=(241 - 0.1, 241.1), xlabel="time", ylabel="y", dpi=1000, legend=false)
display(aa)
savefig("KF_valve_sim_Fig_zoom1.png")
aa = plot!(xlims=(0, 3), ylims=(200, 250), dpi=1000, legend=false)
display(aa)
savefig("KF_valve_sim_Fig_zoom2.png")

println("DONE - All plots saved")
