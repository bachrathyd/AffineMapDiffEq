5+5
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

function DelayMathieu(u, h, p, t)
    ζ, δ, ϵ, b, τ , T = p
    dx = u[2]
    ddx = [-(δ+ϵ*cos(t)), -2*ζ]' * u + b * h(p, t - τ)[1] 
    SA[dx, ddx]
end

Base.rand(x::Vector{T}) where T = rand(T,size(x,1))

#function f2DelayMathieu(v, u, h, p, t)
#    ζ, δ, ϵ, b, τ , T = p
#    ddx = -(δ+ϵ*cos(t)) * u + (-2*ζ) * v + b * h(p, t - τ)[1] 
#end

ζ = 0.02  
δ = 1.5
ϵ = 0.3
τ = 2pi 
b = 0.5
T= 2pi
p = ζ, δ, ϵ, b, τ , T
#p = (ζ, ωn, k, τ,10.0)

u0 = SA[1.0, 0.0]
h(p, t) = SA[1.0, 0.0]
probMathieu = DDEProblem(DelayMathieu, u0, h, (0.0, T * 10.0), p; constant_lags=[τ])
#probMathieu =SecondOrderDDEProblem(f2DelayMathieu, 0.0,0.0,h, (0.0, T * 10.0), p; constant_lags=[τ])

sol = solve(probMathieu, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)
plot(sol(sol.t[end] .- (0.0:0.01:τ*1.0)))

Base.:+(a::SVector, b::Bool) = a .+ b



Nstep = 50
τmax = 2pi+0.1
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax,
 T; Historyresolution=Nstep, eigN=8, zerofixpont=true);

# fix point by affine map
muaff = spectrum(dpdp; p=p);
plot(log.(abs.(muaff[1])))

scatter(muaff[1])
plot!(sin.(0:0.01:2pi),cos.(0:0.01:2pi))


##########using ForwardDiff
##########
##########affine(dpdp; p=p)[1][1]
##########foo(x)=affine(dpdp; p=x)[1][1][1]
##########foo(p)
##########
##########pDual=ForwardDiff.Dual.( [p...])
##########spectrum(dpdp,p=pDual)[1]

######## -------------- test teh ISII implementation ------------------


@benchmark issi_eigen(dpdp)
@benchmark spectralradius(dpdp,p=p)
@code_warntype spectralradius(dpdp,p=p)
issi_eigen(dpdp)


