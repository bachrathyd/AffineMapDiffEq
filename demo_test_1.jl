
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

using MDBM
function delayed_turning_STATIC(u, h, p, t)
    # Parameters
    ζ, ωn, k, τ, f0 = p
    # Components of the delayed differential equation
    dx = u[2]
    # ddx = -2 * ζ * u[2] - ωn * u[1] + k * cos(0*2 * 2 * pi / τ * t) * (f0 * cos(5 * 2 * pi / τ * t) + u[1] - h(p, t - τ)[1])  # Surface regeneration effect
    #SSV
    tau=τ*(1.0-0.2*cos( 2 * pi  * 2*t/1.0) )
    #tau = τ

    ddx = f0*mod(t,0.5) - 2 * ζ * u[2] - ωn * u[1] + k * (cos(2 * pi * t) * 0.0 + 1.0) * (f0 * cos(2 * pi * t) + u[1] - h(p, t - tau)[1])  # Surface regeneration effect
    #ddx = -2 * ζ * u[2] - ωn * u[1] + k * ( u[1] - h(p, t - τ)[1])  # Surface regeneration effect

    # Update the derivative vector
    SA[dx, ddx]
end


##<<<<<<<<<<< Lin Map based on
u0 = SA[0.0, 0.0]
ζ = 0.2          # damping coefficient
ωn = 2.0#0.2          # nat. freq
k = 0.45#4#5#8;#5         # cut.coeff
τ = 2.4          # Time delay
f0 = 1.0         # excitation
p = [ζ, ωn, k, τ, f0]
p0 = [ζ, ωn, k, τ, 0.0]
#p = (ζ, ωn, k, τ,10.0)


T = 2.0;
tspan = (0.0, T * 400.0)

h(p, t) = SA[0.0; 0.0]
probTurning = DDEProblem(delayed_turning_STATIC, u0, h, tspan, p; constant_lags=[τ])

#StateSmaplingTime = LinRange(-τ, 0.0, 200) .+ T

#solve(prob,RK4();dt=0.005,adaptive=false)
sol = solve(probTurning, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)
plot(sol(sol.t[end] .- (0.0:0.01:τ*1.0)))

Nstep = 150
τmax = 3.0
dpdp = dynamic_problemSampled(probTurning, MethodOfSteps(BS3()), τmax,
 T; Historyresolution=Nstep, eigN=10, zerofixpont=true);


## fix point by simulation
OneMap = s -> LinMap(dpdp, s; p=p)


#Fix point
vfix_simulation = zeros(typeof(u0),Nstep);#v0;s0
vfix_simulation = rand(typeof(u0),Nstep);#v0;s0
for _ in 1:500
    vfix_simulation = OneMap(vfix_simulation)
end
norm(vfix_simulation - LinMap(dpdp, vfix_simulation; p=p))
plot(getindex.(vfix_simulation,[1]))
plot!(getindex.(vfix_simulation,[2]))

#Nstep=size(dpdp.StateSmaplingTime, 1);
#s0= rand(typeof(dpdp.DDEdynProblem.u0),Nstep);
#mus = eigsolve(s -> LinMap(dpdp, s; p=p), s0, dpdp.eigN, :LM)


Base.:+(a::SVector, b::Bool) = a .+ b

# fix point by affine map
mus0,As=spectrum(dpdp;p=p0);
plot(log.(abs.(mus0)))
muaff,s0aff=affine(dpdp; p=p);
plot!(log.(abs.(muaff[1])))


plot(getindex.(vfix_simulation,[1]))
plot!(getindex.(vfix_simulation,[2]))
plot!(real.(getindex.(s0aff,[1])) )
plot!(real.(getindex.(s0aff,[2])) )




#@code_warntype affine(dpdp; p=p)
@benchmark affine($dpdp; p=$p)

function foo()
    for _ in 1:1000
        affine(dpdp; p=p)
    end
end

using Profile
using PProf
#@profile  foo();#abstol,reltol
@profview  foo();#abstol,reltol


#Profile.Allocs.clear()
#Profile.Allocs.@profile sample_rate=1  affine(dpdp; p=p)
#PProf.Allocs.pprof(from_c=false)
#Profile.take_heap_snapshot("snapshot.heapsnapshot")














###using MDBM
###
###ax1=Axis([-5,-2.5,0,2.5,5],"k") # initial grid in x direction
###ax2=Axis(0.01:0.5:8.0,"tau") # initial grid in y direction
###
###function foo(k, τ)
###    muaff,s0aff=affine(dpdp; p= [ζ, ωn, k, τ, f0]);
###    return abs(muaff[1][1])-1.0
###end
###mymdbm=MDBM_Problem(foo,[ax1,ax2])
###iteration=5 #number of refinements (resolution doubling)
###MDBM.solve!(mymdbm,iteration)
###
####points where the function foo was evaluated
###x_eval,y_eval=getevaluatedpoints(mymdbm)
####interpolated points of the solution (approximately where foo(x,y) == 0 and c(x,y)>0)
###x_sol,y_sol=getinterpolatedsolution(mymdbm)
###scatter(x_eval,y_eval,s=5)
###scatter(x_sol,y_sol,s=5)