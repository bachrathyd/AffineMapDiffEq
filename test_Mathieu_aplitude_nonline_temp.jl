10 + 10
println(Threads.nthreads())


using Revise
using DDE_mapping

using BenchmarkTools
using Plots
#plotly()
using Profile
using StaticArrays
using DifferentialEquations


##TODO: Bifurcation Analysis in Julia

#using Memoization

using MDBM

#@memoize
function f_now(p, t)
    ζ, δ, ϵ, b, τ, T, μ = p
    # println("Computed $t")
    # println("Computed $p")
    SA[-(δ + ϵ * cos(2pi * t / T)), -2*ζ]
end
#@memoize 
function f_past(p, t)
    ζ, δ, ϵ, b, τ, T, μ = p
    SA[b, 0]
end
#f_now(p, 0.6)' * SA[1.0, 1.0]

function DelayMathieu(u, h, p, t)
    # Parameters
    ζ, δ, ϵ, b, τ, T, μ = p
    # Components of the delayed differential equation
    dx = u[2]
    ddx = f_now(p, t)' * u + abs(u[1])^0.75 + b * h(p, t - τ)[1] + 1.0+μ * (cos(2pi * t / T) .^ 10)  # Surface regeneration effect
    # Update the derivative vector
    SA[dx, ddx]
end
Base.:+(a::SVector, b::Bool) = a .+ b #TODO: where to put this?
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?
Base.:+(b::Float64, a::SVector) = a .+ b #TODO: where to put this?


##<<<<<<<<<<< Lin Map based on

ζ = 0.08          # damping coefficient
δ = 4.05#0.2          # nat. freq
ϵ = 0.15#4#5#8;#5         # cut.coeff
τ = 2pi         # Time delay
b = 0.5
T = 2pi
μ = 0.01#3.3#0.01;
p = ζ, δ, ϵ, b, τ, T, μ
#p = (ζ, ωn, k, τ,10.0)

# test simulation
u0 = SA[0.0, 0.0]
h(p, t) = SA[0.0; 0.0]
probMathieu = DDEProblem(DelayMathieu, u0, h, (0.0, T * 1000.0), p; constant_lags=[τ])

#StateSmaplingTime = LinRange(-τ, 0.0, 200) .+ T

#solve(prob,RK4();dt=0.005,adaptive=false)
sol = solve(probMathieu, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)

Nstep = 150
tFIXend = LinRange(-τ,0.0,Nstep)
uFIXend = sol(sol.t[end] .+ tFIXend).u
plot(tFIXend, getindex.(uFIXend, 1))
plot!(tFIXend, getindex.(uFIXend, 2))



# fix point and spectrum test-----------------------------
Nstep = 150
τmax = 2pi + 0.1
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax,
    T; Historyresolution=Nstep, eigN=10, zerofixpont=false,affineinteration=15,
    adaptive=false,KrylovTol=1e-12,KrylovExtraDim=10);


#muaff, s0aff = affine(dpdp; p=p);# work only if s0 is initialized by rand (and not zero)
muaff, s0aff = affine(dpdp, uFIXend; p=p);
muaff, s0aff = affine(dpdp, s0aff; p=p);
muaff, s0aff = affine(dpdp, s0aff; p=p);
muaff, s0aff = affine(dpdp, s0aff; p=p);
norm(s0aff - LinMap(dpdp, s0aff; p=p))
uFIXend_2=uFIXend
muaff, uFIXend_2 = affine(dpdp, uFIXend_2; p=p);
norm(uFIXend_2 - LinMap(dpdp, uFIXend_2; p=p))
plot(log.(abs.(muaff[1])))

scatter(muaff[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))


plot(tFIXend, getindex.(uFIXend, 1))
plot!(tFIXend, getindex.(uFIXend, 2))

plot!(LinRange(-τmax, 0.0, size(s0aff, 1)), getindex.(s0aff, 1))
plot!(LinRange(-τmax, 0.0, size(s0aff, 1)), getindex.(s0aff, 2))
plot!(LinRange(-τmax, 0.0, size(uFIXend_2, 1)), getindex.(uFIXend_2, 1))
plot!(LinRange(-τmax, 0.0, size(uFIXend_2, 1)), getindex.(uFIXend_2, 2))



#   ---------- starting simulation form a detected solution ----------

StateSmaplingTime = dpdp.StateSmaplingTime

s_for_history=uFIXend # "Chatter"
s_for_history=s0aff
itp = interpolate(s_for_history, BSpline(Cubic(Line(OnGrid()))))
Hist_interp_linear = scale(itp, StateSmaplingTime)
hint(p, t) = Hist_interp_linear(t) #TODO: ha úgyis fix a lépls, akkor ez nem is kell!!!
hint(p, t, deriv::Type{Val{1}}) = Interpolations.gradient(Hist_interp_linear, t)[1]

probMathieu = DDEProblem(DelayMathieu,  hint, (0.0, T * 1000.0), p; constant_lags=[τ])

sol = solve(probMathieu, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)

tPoinCarre = 0:T:T*1000.0
uPoinCarre = sol(tPoinCarre).u
scatter!(tPoinCarre, getindex.(uPoinCarre, 1))
scatter!(tPoinCarre, getindex.(uPoinCarre, 2))



plot!(tFIXend, getindex.(uFIXend, 1))
plot!(tFIXend, getindex.(uFIXend, 2),xlims =(-2pi,2pi))
plot!(tFIXend.+2pi, getindex.(uFIXend, 1))
plot!(tFIXend.+2pi, getindex.(uFIXend, 2),xlims =(-2pi,2pi))


