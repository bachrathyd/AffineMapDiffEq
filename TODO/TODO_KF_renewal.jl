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
gr()
using Profile
using StaticArrays
using DifferentialEquations

using LaTeXStrings
# using Memoization

using MDBM
using MAT
using Peaks
using NaNStatistics
#Kádár Fanni - Szelep rezgés
# Algeb. Delay Diff. Eq.
5 + 5

using Revise
##includet("src\\DDE_mapping.jl")
includet("src\\DDE_mapping_types.jl")
includet("src\\DDE_mapping_functions.jl")
#using DDE_mapping
5 + 5

using BenchmarkTools
using Plots
plotly()
using Profile
using StaticArrays
using DifferentialEquations

using LaTeXStrings
# using Memoization

using MDBM

function valve_KF_neutral(y, h, p, t)
    ζ, δ, τ, β, Φ, q = p

    fτ = h(p, t - τ, Val{1})[4]

    f = fτ + (Φ^2.0 * y[1]^2.0) / 2.0 + Φ * y[1] *
        abs(2 * fτ + y[3] + (Φ^2.0 * y[1]^2.0) / 4.0)^0.5

    d_y1 = y[2]
    d_y2 = -2 * ζ * y[2] - (y[1] + δ) + y[3] - f + fτ
    d_y3 = β * (q - 1 / Φ * (f + fτ))
    d_y4 = f

    SA[d_y1, d_y2, d_y3, d_y4]
end

Base.:+(a::SVector, b::Bool) = a .+ b
ζ = 0.25;
δ = 3.0
τ = pi / 3.0;
β = 10.0
Φ = 48.2
q = 2.0;
4.0;
6.0;
q = 6.0
p = (ζ, δ, τ, β, Φ, q)


#h(p, t) = SA[2.0, 0.0, 2.5, 0]
#h(p, t, ::Type{Val{1}}) = SA[0.0, 0.0, 0.0, q*Φ/2]

h(p, t) = SA[0.0, 0.0, 0.0, 0.0]
h(p, t, ::Type{Val{1}}) = SA[0.0, 0.0, 0.0, q*Φ/2]
@show u0 = h(p, 0.0)#  SA[0.1, 0.1, 0.1, 0.1]
#valve_KF_DADE(u0,h,p,0.0)
T = τ * 5.0;


prob_valve = DDEProblem(valve_KF_neutral, u0, h, (0.0, T), p, constant_lags=[τ], neutral=true)#
#ROS3P 
@time sol = solve(prob_valve, MethodOfSteps(Rodas5()), reltol=1e-9, abstol=1e-9);
#abstol,reltol
#plot(sol)
plot!(sol.t, sol[1:3, :]')
plot!(sol.t, [sol(tloc, Val{1}, idxs=4) for tloc in sol.t])
# plot(sol(sol.t[end] .- (0.0:0.01:τ*1.0)))


tFIXend = (-τ:0.01:0)
uFIXend = sol(sol.t[end] .+ tFIXend).u
plot()
for k in 1:4
    plot!(tFIXend, getindex.(uFIXend, k))
end
plot!()

