5 + 5


using Revise
using DDE_mapping

using BenchmarkTools
using Plots
#plotly()
using Profile
using StaticArrays
using DifferentialEquations

using MDBM

#TODO: SecondOrderDDEProblem(f, [, du0, u0], h, tspan[, p]; <keyword arguments>)
function delayed_turning_STATIC(u, h, p, t)
    # Parameters
    ζ, ωn, k, τ, f0 = p

    tau = τ * (1.0 - 0.2 * cos(2 * pi * 2 * t / 1.0)) #time periodic delay
    # Components of the delayed differential equation
    dx = u[2]


    ddx = f0 * mod(t, 0.5) - 2 * ζ * u[2] - ωn * u[1] + k * (cos(2 * pi * t) * 0.0 + 1.0) * (f0 * cos(2 * pi * t) + u[1] - h(p, t - tau)[1])  # Surface regeneration effect

    SA[dx, ddx]
end


##<<<<<<<<<<< Lin Map based on
u0 = SA[0.0, 0.0]
ζ = 0.2          # damping coefficient
ωn = 2.0#0.2          # nat. freq
k = 0.2#0.85#4#5#8;#5         # cut.coeff
τ = 2.4          # Time delay
f0 = 1.0         # excitation
p = [ζ, ωn, k, τ, f0]
p0 = [ζ, ωn, k, τ, 0.0]
#p = (ζ, ωn, k, τ,10.0)


T = 2.0;
tspan = (0.0, T * 40.0)

h(p, t) = SA[0.0; 0.0]
probTurning = DDEProblem(delayed_turning_STATIC, u0, h, tspan, p; constant_lags=[τ])

#StateSmaplingTime = LinRange(-τ, 0.0, 200) .+ T

#solve(prob,RK4();dt=0.005,adaptive=false)

#Solver_args = Dict(:alg => MethodOfSteps(BS3()))
#Solver_args = Dict()
Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-12)

sol = solve(probTurning; Solver_args...)#abstol,reltol
plot(sol)
plot(sol(sol.t[end] .- (0.0:0.01:τ*1.0)))

τmax = τ * 1.2
#Krylov_arg=Dict(:eigN =>4,:LM, KrylovKit.Arnoldi(KrylovTol=1e-12,KrylovExtraDim=5,verbosity=3));
Krylov_arg = (6, :LM, KrylovKit.Arnoldi());


Nstep = 1000
dpdp = dynamic_problemSampled(probTurning, Solver_args, τmax,
    T; Historyresolution=Nstep,
    zerofixpont=true, affineinteration=1,
    Krylov_arg=Krylov_arg)



## fix point by simulation
OneMap = s -> LinMap(dpdp, s; p=p)[1]


#Fix point
vfix_simulation = zeros(typeof(u0), Nstep);#v0;s0
vfix_simulation = rand(typeof(u0), Nstep);#v0;s0
for _ in 1:500
    vfix_simulation = OneMap(vfix_simulation)
end
norm(vfix_simulation - LinMap(dpdp, vfix_simulation; p=p)[1])
plot(getindex.(vfix_simulation, [1]))
plot!(getindex.(vfix_simulation, [2]))

#Nstep=size(dpdp.StateSmaplingTime, 1);
#s0= rand(typeof(dpdp.DDEdynProblem.u0),Nstep);
#mus = eigsolve(s -> LinMap(dpdp, s; p=p), s0, dpdp.eigN, :LM)


Base.:+(a::SVector, b::Bool) = a .+ b

# fix point by affine 
@time muaff, s0aff, solPer = affine(dpdp; p=p);
plot(log.(abs.(muaff[1])))


plot(dpdp.StateSmaplingTime, getindex.(vfix_simulation, [1]))
plot!(dpdp.StateSmaplingTime, getindex.(vfix_simulation, [2]))
plot!(dpdp.StateSmaplingTime, real.(getindex.(s0aff, [1])))
plot!(dpdp.StateSmaplingTime, real.(getindex.(s0aff, [2])))

plot!(solPer.t, real.(getindex.(solPer.u, [1])))
plot!(solPer.t, real.(getindex.(solPer.u, [2])))




@benchmark affine($dpdp; p=$p)



## ------------- stablitly chart  #CPU time: ~1000 sec ---------
using MDBM

ax1 = Axis(1:0.5:8.0, "tau") # initial grid in y direction
ax2 = Axis(LinRange(-20.0, 5.0, 6), "k") # initial grid in x direction

using KrylovKit
Krylov_arg = (4, :LM, KrylovKit.Arnoldi(tol=1e-13, krylovdim=4 + 5, verbosity=0));

function foo(τ, k)

    Nstep = Int(floor(τ * 40))
    dpdp = dynamic_problemSampled(probTurning, Solver_args, τ * 1.2,
        T; Historyresolution=Nstep,
        zerofixpont=true, affineinteration=2,
        Krylov_arg=Krylov_arg)

    muaff, s0aff, solPer = affine(dpdp; p=[ζ, ωn, k, τ, f0])

    return abs(muaff[1][1]) - 1.0
end
mymdbm = MDBM_Problem(foo, [ax1, ax2])
iteration = 2 #number of refinements (resolution doubling)
@time MDBM.solve!(mymdbm, iteration);

#points where the function foo was evaluated
x_eval, y_eval = getevaluatedpoints(mymdbm)
#interpolated points of the solution (approximately where foo(x,y) == 0 and c(x,y)>0)
x_sol, y_sol = getinterpolatedsolution(mymdbm)
scatter(x_eval, y_eval)
scatter(x_sol, y_sol)