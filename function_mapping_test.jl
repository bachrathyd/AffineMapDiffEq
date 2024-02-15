
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



#Base.:*(α, v): multiply v with a scalar α, which can be of a different scalar type; in particular this method is used to create vectors similar to v but with a different type of underlying scalars.
function Base.:*(α::Number, v::ODESolution)
    vout = deepcopy(sol)
    vout.u .*= α
    vout.k .*= α
    return vout::ODESolution
end
#Base.similar(v): a way to construct vectors which are exactly similar to v
function Base.similar(v::ODESolution)::ODESolution
    vout = deepcopy(v)
    vout.u .*= 0.0
    vout.k .*= 0.0
    return vout::ODESolution
end
#LinearAlgebra.mul!(w, v, α): out of place scalar multiplication; multiply vector v with scalar α and store the result in w

#LinearAlgebra.rmul!(v, α): in-place scalar multiplication of v with α; in particular with α = false, v is the corresponding zero vector
function LinearAlgebra.rmul!(v::ODESolution, α::Number)::ODESolution
    v.u .*= α
    v.k .*= α
    return sol::ODESolution
end
function LinearAlgebra.rmul!(sol::ODESolution, α::Bool)::ODESolution
    println(α)
    sol.u .*= α#0.0
    sol.k .*= α#0.0
    return sol::ODESolution
end
#LinearAlgebra.axpy!(α, v, w): store in w the result of α*v + w
function LinearAlgebra.axpy!(α::Number, v::ODESolution, w::ODESolution)::ODESolution
    
    h_axpy(p, t) = α * v.prob.h(p, t) + w.prob.h(p, t)
    u0_axpy= α * v.prob.u0 + w.prob.u0
    if sol1.prob.neutral
    h_axpy(p, t, deriv::Type{Val{1}}) =  α * v.prob.h(p, t,Val{1}) + w.prob.h(p, t,Val{1}) 
    end
    sol_axpy = solve(remake(sol1.prob;u0=u0_axpy,h=h_axpy), MethodOfSteps(BS3()))#TODO use the original method
    w=sol_axpy
    #w.u=sol_axpy.u
    #w.t=sol_axpy.t
    #w.k=sol_axpy.k
    #return w::ODESolution
end

#LinearAlgebra.axpby!(α, v, β, w): store in w the result of α*v + β*w

using Integrals
#LinearAlgebra.dot(v,w): compute the inner product of two vectors
function LinearAlgebra.dot(v::ODESolution, w::ODESolution)::Float64
    prob = IntegralProblem((t, p) -> v(t)' * w(t), v.t[1], v.t[end])
    VW = solve(prob, HCubatureJL()).u# reltol=1e-5, abstol=1e-5r
    return VW::Float64
end

push!(sol.u, sol.u[1])
#LinearAlgebra.norm(v): compute the 2-norm of a vector
function LinearAlgebra.norm(v::ODESolution)::Float64
    return sqrt(dot(v, v))::Float64
end

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
plot!(sol2RE)

plot(diff(sol1.t))
plot!(diff(sol2.t))


prob = IntegralProblem((t, p) -> sol1(t)' * sol2(t), sol1.t[1], sol1.t[end])
VW = solve(prob, HCubatureJL()).u;# reltol=1e-5, abstol=1e-5r


sol1.u[1]' * sol2.u[4]
Base.:+(a::SVector, b::Bool) = a .+ b


mus = getindex(schursolve(s -> LinMap(dp, s; p=p), s_start, dp.eigN, :LM, KrylovKit.Arnoldi(krylovdim=dp.eigN + dp.KrylovExtraDim, tol=dp.KrylovTol, verbosity=0)), [3, 2, 1])

