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
    ζ, δ, ϵ, b, τ , T, μ  = p
   # println("Computed $t")
   # println("Computed $p")
    SA[-(δ+ϵ*cos(2pi*t/T)), -2*ζ]
end
#@memoize 
function f_past(p, t)
    ζ, δ, ϵ, b, τ , T, μ  = p
    SA[b, 0]
end
#f_now(p, 0.6)' * SA[1.0, 1.0]

function DelayMathieu(u, h, p, t)
    # Parameters
    ζ, δ, ϵ, b, τ , T, μ = p
    # Components of the delayed differential equation
    dx = u[2]
    ddx = f_now(p, t)' * u +  u[1]^3 + b * h(p, t - τ)[1] + μ *(cos(2pi*t/T).^10)  # Surface regeneration effect
    # Update the derivative vector
    SA[dx, ddx]
end
Base.:+(a::SVector, b::Bool) = a .+ b #TODO: where to put this?
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?


##<<<<<<<<<<< Lin Map based on

ζ = 0.08          # damping coefficient
δ = 4.05#0.2          # nat. freq
ϵ = 0.15#4#5#8;#5         # cut.coeff
τ = 2pi         # Time delay
b = 0.5
T= 2pi
μ=0.01;
p = ζ, δ, ϵ, b, τ , T, μ
#p = (ζ, ωn, k, τ,10.0)

# test simulation
u0 = SA[0.0, 0.0]
h(p, t) = SA[0.0; 0.0]
probMathieu = DDEProblem(DelayMathieu, u0, h, (0.0,T * 1000.0), p; constant_lags=[τ])

#StateSmaplingTime = LinRange(-τ, 0.0, 200) .+ T

#solve(prob,RK4();dt=0.005,adaptive=false)
sol = solve(probMathieu, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)

tFIXend=(- τ:0.01:0)
uFIXend=sol(sol.t[end] .+ tFIXend).u
plot(tFIXend,getindex.(uFIXend,1))
plot!(tFIXend,getindex.(uFIXend,2))


# fix point and spectrum test
Nstep = 150
τmax = 2pi+0.5
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax,
 T; Historyresolution=Nstep, eigN=10, zerofixpont=false);


 muaff, s0aff = affine(dpdp,s0aff; p=p);
 plot(log.(abs.(muaff[1])))
 
 scatter(muaff[1])
 plot!(sin.(0:0.01:2pi),cos.(0:0.01:2pi))
 
 
 plot(tFIXend,getindex.(uFIXend,1))
 plot!(tFIXend,getindex.(uFIXend,2))
 plot!(LinRange(-τmax,0.0,size(s0aff,1)),getindex.(s0aff,1))
 plot!(LinRange(-τmax,0.0,size(s0aff,1)),getindex.(s0aff,2))
 

# fix point and spectrum test





#------------ Mu test for different resolution 
#("Float Mapp" or "Dual mapp") test - the source code should be changed inside!!!

Nstep = 50
τmax = 2pi+0.5
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax,
 T; Historyresolution=Nstep, eigN=1, zerofixpont=false);


EPSI_TODO_REMOVE=0.000000001 #TODO: remove this variable!!!!!!!4
muaff, s0aff = affine(dpdp; p=p);
#plot(log.(abs.(muaff[1])))
# fix point by affine map
EPSvpow=-2.0:-0.1:-20.0##-10.0:0.1:0.2
#EPSvpow=-10.0:0.1:-2.0##-10.0:0.1:0.2
EPSvpow=-2.0:0.1:3.0
musABS=[]
scatter([],[])
@time for pp=EPSvpow
    EPSI_TODO_REMOVE=10^pp
    println(EPSI_TODO_REMOVE)
    muaff, s0aff = affine(dpdp,s0aff; p=p);
    #plot!(log.(abs.(muaff[1])))
    push!(musABS,abs(muaff[1][1]))
    scatter!([pp],[musABS[end]],legend=false)
end

scatter!([],[])
scatter!(EPSvpow,musABS)
norm(diff(musABS))



### #TODO: ArrayPartition(x::AbstractArray...)



#--------------- bifurcation test - "continuation" -------------
#EPSI_TODO_REMOVE=1e-15
#TODO: Comaper the Float and Dual Mapp resolution --> result for a jounal
μv=-0.0:0.01:2.0 # initial grid in x direction
#μv=-0.0:-0.01:-2.0 # initial grid in x direction
Aaff=zeros(size(μv,1))
Spek_aff=zeros(size(μv,1))


muaff, s0aff = affine(dpdp; p=(ζ, δ, ϵ, b, τ, T,0.0))
#Threads.@threads
@time  for j in 1:size(μv, 1)
    #@inbounds 
    println(j/size(μv, 1))
     μ = μv[j]
    # muaff, s0aff = affine(dpdp; p=(ζ, δ, ϵ, b, τ, T,μ))
     muaff, s0aff = affine(dpdp,s0aff; p=(ζ, δ, ϵ, b, τ, T,μ))
        Aaff[j] = norm(getindex.(s0aff, 1))
        # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1)))
        Spek_aff[ j] = maximum(abs.(muaff[1]))
end

#@code_warntype affine(dpdp; p=(ζ, δ, ϵ, b, τ, T,μ))
theme(:ggplot2)

Aaffsat=deepcopy(Aaff);
#Aaffsat[Spek_aff .> 1.0] .= 0.0;
plot(μv,Aaffsat)
plot!(μv,Spek_aff)
scatter!(μv,Aaffsat,zcolor=3.0 .- sign.(Spek_aff .- 1.0) )
#heatmap(δv,bv,(Aaffsat))