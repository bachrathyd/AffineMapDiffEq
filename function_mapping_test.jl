
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
#TODO:  @fastmath @inbounds
#@memoize
function f_now(p, t)
    ζ, δ, ϵ, b, τ, T = p
    # println("Computed $t")
    # println("Computed $p")
    SA[-(δ + ϵ * cos(t)), -2*ζ]
    #SA[-(δ+ϵ*sign(cos(t))), -2*ζ]
end
#@memoize 
function f_past(p, t)
    ζ, δ, ϵ, b, τ, T = p
    SA[b, 0]
end

function DelayMathieu(u, h, p, t)
    ζ, δ, ϵ, b, τ, ν, T = p
    dx = u[2]
    #ddx = [-(δ+ϵ*cos(t)), -2*ζ]' * u + b * h(p, t - τ)[1] # Traditional Delayed Mathieu equation
    #sign
    ddx = [-(δ + ϵ * (cos(t))), -2 * ζ]' * u + b * h(p, t - τ)[1] + ν * h(p, t - τ, Val{1})[2]#The last pert is a test for neutral system
    SA[dx, ddx]
end


ζ = 0.03
δ = 0.75
ϵ = 0.2
τ = 2pi
b = 0.2
ν = 0.25 # Neutral coefficient
T = 6;
2pi;
p = ζ, δ, ϵ, b, τ, ν, T
#p = (ζ, ωn, k, τ,10.0)

u0_1 = SA[1.0, 1.0]
u0_2 = SA[-1.0, 1.0]
h(p, t::Float64) = SA[0.0; -0.0]
h(p, t::Float64, deriv::Type{Val{1}}) = SA[0.0, 0.0]
probMathieu = DDEProblem(DelayMathieu, u0_1, h, (0.0, T * 10.0), p; constant_lags=[τ], neutral=true)
#TODO: csak ezzel működik dense=true
sol1 = solve(remake(probMathieu, u0=u0_1), MethodOfSteps(BS3()))#abstol,reltol
sol2 = solve(remake(probMathieu, u0=u0_2), MethodOfSteps(BS3()))#abstol,reltol




sol1.u[100] = u0_1 * 10

plot(sol)
plot(0.0001 * sol1)
plot(similar(sol1))
plot(LinearAlgebra.rmul!(sol1, true))


sol1 = solve(remake(probMathieu, u0=u0_1), MethodOfSteps(BS3()))#abstol,reltol
plot(sol1)
sol2 = solve(remake(probMathieu, u0=u0_2), MethodOfSteps(BS3()))#abstol,reltol
plot!(sol2)
sol2RE=deepcopy(sol2)
LinearAlgebra.axpy!(1.0, sol1, sol2RE)
sol2RE
plot!(sol2RE)

plot(diff(sol1.t))
plot!(diff(sol2.t))


prob = IntegralProblem((t, p) -> sol1(t)' * sol2(t), sol1.t[1], sol1.t[end])
VW = solve(prob, HCubatureJL()).u;# reltol=1e-5, abstol=1e-5r


sol1.u[1]' * sol2.u[4]
Base.:+(a::SVector, b::Bool) = a .+ b


mus = getindex(schursolve(s -> LinMap(dp, s; p=p), s_start, dp.eigN, :LM, KrylovKit.Arnoldi(krylovdim=dp.eigN + dp.KrylovExtraDim, tol=dp.KrylovTol, verbosity=0)), [3, 2, 1])

