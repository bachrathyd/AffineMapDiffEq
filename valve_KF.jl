#Kádár Fanni - Szelep rezgés
# Algeb. Delay Diff. Eq.
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

using LaTeXStrings
# using Memoization

using MDBM

function valve_KF_DADE(y, h, p, t)
    ζ, δ, τ,β,Φ,q= p

    d_y1 = y[2]
    d_y2 = -2*ζ*y[2]-(y[1]+δ)+y[3]-y[4]+h(p, t - τ)[4] 
    d_y3 =β * (q-1/Φ*(y[4]+h(p, t - τ)[4] )) 
    CORESRQ=(y[3]-y[4]+h(p, t - τ)[4])
    #CORESRQ=maximum([0,CORESRQ])
   # d_X= Φ* y[1]*  abs( CORESRQ) ^0.5 - h(p, t - τ)[4]-y[4]#Singular "f=y4"
    d_X= -y[4] - h(p, t - τ)[4]-(Φ^2.0* y[1]^2.0)/2.0+Φ* y[1]*  abs( 2*h(p, t - τ)[4]+y[3]+(Φ^2.0* y[1]^2.0)/4.0) ^0.5#Singular "f=y4"
    d_y5=d_X
    #d_y6=CORESRQ
#    SA[d_y1,d_y2,d_y3,d_X,d_y5]
    SA[d_y1,d_y2,d_y3,d_X]
end

Base.:+(a::SVector, b::Bool) = a .+ b
ζ=0.25;
δ=3.0
τ=pi/3.0;
β=10.0
Φ=48.2
q=2.0;4.0;6.0
q=6.0
p=(ζ, δ, τ,β,Φ,q)


#u0 = SA[0.1, 0.1, 0.1, 0.1,0.0]
h(p, t) = SA[2.0,  0.0, 2.5,q*Φ/2,0.0]
h(p, t) = SA[2.0,  0.0, 2.5,q*Φ/2]
@show u0 =h(p, 0.0)#  SA[0.1, 0.1, 0.1, 0.1]
#valve_KF_DADE(u0,h,p,0.0)
T=τ*5.0;
#M=diagm([1.0,1.0,1.0,0.00000001,1.0])
M=diagm([1.0,1.0,1.0,1.0])
M=diagm([1.0,1.0,1.0,0.01])
M=diagm([1.0,1.0,1.0,0.0001])
M=diagm([1.0,1.0,1.0,0.000001]) #1e-6
M=diagm([1.0,1.0,1.0,0.0]) #1e-9

#y=u0
#t=0.0


valve_KF_DADE_fun=DDEFunction(valve_KF_DADE,mass_matrix=M)
prob_valve = DDEProblem(valve_KF_DADE_fun, u0, h, (0.0, T ), p,constant_lags=[τ])
#ROS3P 
@time sol = solve(prob_valve, MethodOfSteps(Rodas5()),
 reltol = 1e-6, abstol = 1e-6);#abstol,reltol
 #plot(sol)
plot(sol)
# plot(sol(sol.t[end] .- (0.0:0.01:τ*1.0)))


tFIXend = (-τ:0.01:0)
uFIXend = sol(sol.t[end] .+ tFIXend).u
plot()
for k in 1:4
plot!(tFIXend, getindex.(uFIXend, k))
end
plot!()




#---------------------- Spectrum test brute-force--------------------
Nstep = size(uFIXend,1)
s_start = rand(typeof(u0), Nstep) .*0.1 .+ uFIXend 
T= τ
τmax= τ
dpdp = dynamic_problemSampled(prob_valve, MethodOfSteps(Rodas5()), τmax,
    T; Historyresolution=Nstep, eigN=15, zerofixpont=false, dt=0.01);
EPSI_TODO_REMOVE=1e-9
@time mu, saff = affine(dpdp,uFIXend; p=p);
@time mu, saff = affine(dpdp,s_start; p=p);
@time mu, saff = affine(dpdp,saff; p=p);
plot(log.(abs.(mu[1])))
scatter((mu[1]))
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))

plot()
for k in 1:size(u0,1)
plot!(tFIXend, getindex.(saff, k))
end
plot!()




#--------------- bifurcation test - "continuation" -------------
Nstep = size(uFIXend,1)
s_start = rand(typeof(u0), Nstep) .*0.1 .+ uFIXend 
T= τ
τmax= τ
dpvalve = dynamic_problemSampled(prob_valve, MethodOfSteps(Rodas5()), τmax,
    T; Historyresolution=Nstep, eigN=1, zerofixpont=false, dt=0.01)#, affineinteration=3);


EPSI_TODO_REMOVE=1e-4

qv = collect(1.0:0.2:30)# initial grid in x direction
qv = collect(6.0:10:100)# initial grid in x direction


#Av=-0.0:-0.01:-2.0 # initial grid in x direction
Aaff = zeros(size(qv, 1))
Spek_aff = zeros(size(qv, 1))


@time mu, s0aff = affine(dpvalve,s_start; p=(ζ, δ, τ,β,Φ,q));
#@time mu, s0aff = affine(dpvalve,s_start; p=(ζ, δ, τ,β,Φ,2.0));
scatter((mu[1]))
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))

#Threads.@threads
scatter([], [])
doanimation = false#true#false

if doanimation
    gr()
end
a = Animation()
@time for j in 1:size(qv, 1)

    println(j / size(qv, 1))
    #q = qv[j]
   # muaff, s0aff = affine(dpvalve, s0aff; p=(ζ, δ, τ,β,Φ,q))
    muaff, s0aff = affine(dpvalve, s0aff; p=(ζ, δ, τ,β,Φ,qv[j]))
    println((ζ, δ, τ,β,Φ,qv[j]))
    #Aaff[j] = norm(getindex.(s0aff, 1))
    #Aaff[j] = sum(getindex.(s0aff, 1).^2)/Nstep
    @show Aaff[j] = maximum(getindex.(s0aff, 3))
    # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1)))
    @show Spek_aff[j] = maximum(abs.(muaff[1]))

    if doanimation
        #scatter(log.(muaff[1]),lab="")
        scatter(muaff[1], lab="")
        #plt=plot!(sin.(0:0.01:2pi),cos.(0:0.01:2pi),lab="")
        plt = plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi), xlim=(-3.0, 1.5), ylim=(-1.5, 1.5), lab="")
        frame(a, plt)# Not working with : #plotly()
    else
        scatter!(muaff[1], lab="")
    end
end
scatter!([], [])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))

if doanimation
    gif(a, fps=25)
end

#@code_warntype affine(dpvalve; p=(ζ, δ, ϵ, b, τ, T,A,μ))
theme(:ggplot2)

Aaffsat = deepcopy(Aaff);
#Aaffsat[Aaffsat .== maximum(Aaffsat)] .= 0.0
#Aaffsat[Spek_aff .> 1.0] .= 0.0;
plot(qv, Aaffsat,lab="")
plot!(qv, Spek_aff,lab="")
scatter!(qv, Aaffsat, zcolor=3.0 .- sign.(Spek_aff .- 1.0),lab="",
xlabel="q", ylabel="y_3 & mu_1")
#heatmap(δv,bv,(Aaffsat))
#, ylim=(0.0, 2.0)


#------------ Map ----------------
println("----------Start brute-force---------------")
ζ=0.5;
δ=3.0;
  τ=pi/3.0;
  β=10.0;
Φ=48.2;
q=6.0;
p=(ζ, δ, τ,β,Φ,q)
τv = LinRange(0.01,3pi,50)# initial grid in x direction
#τv = LinRange(0.1,3pi,60)# initial grid in x direction
βv = LinRange(-10.0,+150.01,50) # initial grid in y direction
Spek = zeros(size(βv, 1), size(τv, 1))
#Threads.@threads
@time Threads.@threads for j in 1:size(τv, 1)
    println(j)
    τ = τv[j]


    Nstep = size(uFIXend,1)
    s_start = rand(typeof(u0), Nstep) .*0.1 .+ uFIXend 
    T= τ;# 0.1;
    τmax= τ
    #dpvalve = dynamic_problemSampled(prob_valve, MethodOfSteps(Rodas5()), τmax,
    #    T; Historyresolution=Nstep, eigN=1, zerofixpont=false, dt=0.01, affineinteration=3);
    
        dpvalve = dynamic_problemSampled(prob_valve, MethodOfSteps(Rodas5()), τmax,
        T; Historyresolution=Nstep, eigN=1, zerofixpont=false, dt=τ/20.0, affineinteration=3);
    
        
    #@time 
    Threads.@threads for i in 1:size(βv, 1)
        β = βv[i]
        mu, saff =affine(dpvalve,s_start; p=(ζ, δ, τ,β,Φ,q));
        muMAX = abs(mu[1][1])
        #muMAX = spectralradius(dploc; p=(ωn, τf, τr, w, μ))
        Spek[i, j] = muMAX
    end
end

Spek_sat = deepcopy(Spek);
Spek_sat[Spek_sat.>1.0] .= 1.0;
Spek_sat[Spek_sat .<0.01] .= 1.0;
#Spek_sat[Spek_sat.>1.0] .= 0.0;
#Spek_sat[Spek_sat .<0.01] .= 0.0;
heatmap(τv, βv, (Spek_sat),xlabel="τ", ylabel="β",title=ζ)

