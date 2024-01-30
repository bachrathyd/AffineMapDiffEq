10+10

println(Threads.nthreads())
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


##TODO: Bifurcation Analysis in Julia

#using Memoization

using MDBM

#@memoize
function f_now(p, t)
    ζ, δ, ϵ, b, τ , T  = p
   # println("Computed $t")
   # println("Computed $p")
    SA[-(δ+ϵ*cos(2pi*t/T)), -2*ζ]
end
#@memoize 
function f_past(p, t)
    ζ, δ, ϵ, b, τ , T  = p
    SA[b, 0]
end
#f_now(p, 0.6)' * SA[1.0, 1.0]

function DelayMathieu(u, h, p, t)
    # Parameters
    ζ, δ, ϵ, b, τ , T = p
    # Components of the delayed differential equation
    dx = u[2]
    ddx = f_now(p, t)' * u + b * h(p, t - τ)[1] + 0.01*(cos(2pi*t/T).^10)  # Surface regeneration effect
    # Update the derivative vector
    SA[dx, ddx]
end
Base.:+(a::SVector, b::Bool) = a .+ b


##<<<<<<<<<<< Lin Map based on
NF=Float64
ζ = NF(0.02)          # damping coefficient
δ = NF(2.5)#0.2          # nat. freq
ϵ = NF(0.15)#4#5#8;#5         # cut.coeff
τ = NF(2pi)          # Time delay
b = NF(0.5)
T= NF(2pi)
p = ζ, δ, ϵ, b, τ , T
#p = (ζ, ωn, k, τ,10.0)

# test simulation
u0 = SA{NF}[1.0, 0.0]
h(p, t) = SA{NF}[0.0; 0.0]
probMathieu = DDEProblem(DelayMathieu, u0, h, (NF(0.0), NF(T * 100.0)), p; constant_lags=[τ])

#StateSmaplingTime = LinRange(-τ, 0.0, 200) .+ T

#solve(prob,RK4();dt=0.005,adaptive=false)
sol = solve(probMathieu, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)
plot(sol(sol.t[end] .- (0.0:0.01:τ*1.0)))


Nstep = 20
τmax = NF(2pi+0.1)
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax,
 T; Historyresolution=Nstep, eigN=1, zerofixpont=false);

# fix point by affine map
muaff, s0aff = affine(dpdp; p=p);
plot(log.(abs.(muaff[1])))

scatter(muaff[1])
plot!(sin.(0:0.01:2pi),cos.(0:0.01:2pi))


plot(getindex.(s0aff,1))
plot!(getindex.(s0aff,2))



println("----------Start brute-force---------------")
δv=0:0.101:10 # initial grid in x direction
bv=-1.501:0.1:1.5 # initial grid in y direction
Aaff=zeros(size(bv,1),size(δv,1))
Spek_aff=zeros(size(bv,1),size(δv,1))

@time Threads.@threads for j in 1:size(δv, 1)
    @inbounds  δ = δv[j]
      Threads.@threads for i in 1:size(bv, 1)
        @inbounds  b = bv[i]
        muaff, s0aff = affine(dpdp; p=(ζ, δ, ϵ, b, τ, T))
        Aaff[i, j] = norm(getindex.(s0aff, 1))
        # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1)))
        Spek_aff[i, j] = maximum(abs.(muaff[1]))
    end
end



Aaffsat=deepcopy(Aaff);
Aaffsat[Spek_aff .> 1.0] .= 0.0;
heatmap(δv,bv,log.(Aaffsat))
#heatmap(δv,bv,(Aaffsat))



Spek_affsat=deepcopy(Spek_aff);
Spek_affsat[Spek_affsat .> 1.0] .= 0.0;
heatmap(δv,bv,log.(Spek_affsat))



ax1=Axis(0:2:10,"δ") # initial grid in x direction
ax2=Axis(-1.5:1.4:1.5,"b") # initial grid in y direction
function fooDelay(δ, b)
    ABSmuMax=spectralradius(dpdp;  p = (ζ, δ, ϵ, b, τ , T));
    return ABSmuMax-1.0
end
mymdbm=MDBM_Problem(fooDelay,[ax1,ax2])
@time MDBM.solve!(mymdbm,4)
#points where the function foo was evaluated
x_eval,y_eval=getevaluatedpoints(mymdbm)
x_sol,y_sol=getinterpolatedsolution(mymdbm)
#scatter(x_eval,y_eval,markersize=1)
#scatter!(x_sol,y_sol,markersize=2)
scatter!(x_sol,y_sol,markersize=1)



