5 + 5
#using Revise
##includet("src\\DDE_mapping.jl")
#includet("src\\DDE_mapping_types.jl")
#includet("src\\DDE_mapping_functions.jl")
using DDE_mapping
using Interpolations

##] add https://github.com/bachrathyd/AffineMapDiffEq.git
#using AffineMapDiffEq

using StaticArrays
using DifferentialEquations
using Plots
#plotly()
gr()
using MAT
using LinearAlgebra
#using BenchmarkTools
#using Profile
#using MDBM
using LaTeXStrings


function delayed_turning_STATIC_SSV(u, h, p, t)
    # Parameters
    ζ, ωn, k, w, K, OMrel, RVA, RVF, f0 = p

    T0 = 2pi / OMrel
    #Om_actual=  (1.0 + RVA * sin(t * OMrel * RVF))
    tau = T0 / (1.0 + RVA * sin(t * OMrel * RVF))
    h0 = f0 *tau / T0
    #if mod(t/T,2.0*pi)<0.4
    #    h0=f0*sin(2.0*pi*t/T)
    #else
    #    h0=0.0;
    #end

    #ddx = ωn^2 * (K*w/k * (h(p, t - tau)[1] - u[1] - h0) - u[1]) -2 * ζ * ωn * u[2] # Surface regeneration effect

    #fi = (OMrel*t)
    fi  = OMrel*t - RVA/RVF * (cos(RVF*OMrel*t) - 1)
    ϕ_mod = mod(fi, 2*pi)
    g   = all([0 <ϕ_mod,ϕ_mod < deg2rad(10)]) ? 0 : 1
    dx  = u[2]
    ddx = ωn^2 * (g*K*w/k * (h(p, t - tau)[1] - u[1] - h0) - u[1]) -2 * ζ * ωn * u[2] 

    
    SA[dx, ddx]
end

Base.:+(a::SVector, b::Bool) = a .+ b

##<<<<<<<<<<< Lin Map based on
u0 = SA[0.0, 0.0]

m  = 2.9214;
k  = 5.8995e+6;
c  = 40.4073;

ωn = sqrt(k/m)      # nat. freq
ζ  = c/(2*ωn*m);    # damping coefficient
K  = 431.25e+6      # cut.coeff
Tn = 2*pi/ωn

OMrel = 170; #rad/s
dt_limit=1/ωn/50

#------------------------
RVA   = 0.05;
RVF   = 0.05;

#RVA   = 0.1;
#RVF   = 1;

#For Fig.1 and Fig.2B: SSV
RVA   = 0.1;
RVF   = 1.0/100.0;

#For Fig.2B: noSSV
RVA   = 0.0;
RVF   = 1.0;

#For Fig.2B: noSSV
RVA   = 0.1;
RVF   = 1.0/10.0;

@show lobenumber_range=ωn ./ (OMrel .* (1.0 .+ [RVA,-RVA]))
#------------------------

T     = 2pi / (OMrel * RVF)
w     = 0.6/1e3;
w     = 0.2/1e3;

f0 = 0.0001        # feed [m/s]

p = [ζ, ωn, k, w, K, OMrel, RVA, RVF, f0]

tspan = (0.0, T * 1.0)

h(p, t) = SA[0.0; 0.0]
#algMoS=MethodOfSteps(Tsit5())
algMoS=MethodOfSteps(BS3())
Nsimulation_max=30
#Egyetlen hosszú szimuláció
probTurningSSV_sim = DDEProblem(delayed_turning_STATIC_SSV, u0, h, (0.0, T * Nsimulation_max), p)
#@time sol1 = solve(probTurningSSV_sim, alg = algMoS,dt=dt_limit,reltol=1e-10,maxiters=1e7,saveat=0:dt_limit:T * Nsimulation_max)
@time sol_short = solve(probTurningSSV_sim, alg = algMoS,dt=dt_limit,reltol=1e-10,maxiters=1e7,
saveat=T * (Nsimulation_max-2):dt_limit:T * Nsimulation_max)

@time sol_all = solve(probTurningSSV_sim, alg = algMoS,dt=dt_limit,reltol=1e-10,maxiters=1e7,
saveat=0:dt_limit:T * Nsimulation_max)

plot(sol_all.t,getindex.(sol_all.u,1),size=(560, 220),
xlabel ="",#time - t [s]
ylabel=L"tool pos - x [m]",legend=false)
#plot(sol_short.t,getindex.(sol_short.u,1),size=(560, 220))

aa = plot!( dpi=500, legend=false)
display(aa)
savefig("Long_sim_Fig_1a.png")


plot(sol_short.t,getindex.(sol_short.u,1),size=(560, 220),
xlabel ="",
ylabel=L"tool pos - x [m]",legend=false)#
#plot(sol_short.t,getindex.(sol_short.u,1),size=(560, 220))

aa = plot!( dpi=500, legend=false)
display(aa)
#savefig("Long_sim_Fig_1b.png")
# --------------------------------------------------------------------------------------

# ------------- affine map --------------------
probTurningSSV = DDEProblem(delayed_turning_STATIC_SSV, u0, h, tspan, p)
τmax = 2pi / OMrel / (1.0 - RVA)
Nstep= Int(ceil(τmax/RVF /  dt_limit))
dpdp = dynamic_problemSampled(probTurningSSV, algMoS, τmax/RVF, T;Historyresolution=Nstep, dt=dt_limit, eigN=2,KrylovTol=1e-11,zerofixpont=false, affineinteration=3);

@time muaff, s0aff = affine(dpdp; p=[ζ, ωn, k, w, K, OMrel, RVA, RVF, f0]);

#------------------------------------Plot fixpoints----------------------------------------------

tvec3 = dpdp.StateSmaplingTime
#plot!(tvec3.+Nsimulation_max*T, real.(getindex.(s0aff, [2])))
aa=plot!(tvec3.+Nsimulation_max*T, real.(getindex.(s0aff, [1])), linestyle=:dash)
display(aa)
savefig("Long_sim_Fig_1b.png")
#plot!(xlim=[tvec3[1],tvec3[end]].+Nsimulation_max*T)
aa=plot!(xlim=[108,109], dpi=500, legend=false,size=(560, 220),xlabel =L"time - t [s]")
display(aa)
savefig("Long_sim_Fig_1c.png")

##---------------------

#@show peak_to_peak = maximum(getindex.(vfix_simulation,  [1])) - minimum(getindex.(vfix_simulation,  [1]))
#@show peak_to_peak = maximum(real.(getindex.(s0aff,  [1]))) - minimum(real.(getindex.(s0aff,  [1])))

#-----------------------------------------------------------------------------------------------

scatter(1:size(muaff[1],1),log.(abs.(muaff[1])))
scatter(muaff[1].^RVF)
plot!(sin.(0:0.01:2pi),cos.(0:0.01:2pi), aspect_ratio=:equal)


###-----------------Simulation of multiple rotations to test fixpoint-------------------------------
##
##RVA = 0.5
##RVF = 1/3
##
##muaff, s0aff = affine(dpdp; p=[ζ, ωn, k, w, K, OMrel, RVA, RVF, f0]);
##
##tspan = (0.0, T * 1.0)
##u0 = SA[0.0, 0.0]
##h(p, t) = SA[0.0; 0.0]
##probTurningSSV = DDEProblem(delayed_turning_STATIC_SSV, u0, h, tspan, p)
##sol = solve(probTurningSSV, MethodOfSteps(BS3()))
##plot(sol.t,getindex.(sol.u,  [1]))
##
##nTurn = 100
##for i in 1:nTurn
##    h(p, t) = sol(t)
##    u0 = sol.u[end,1]
##    tspan = tspan .+ T
##    probTurningSSV = DDEProblem(delayed_turning_STATIC_SSV, u0, h, tspan, p)
##    sol = solve(probTurningSSV, MethodOfSteps(BS3()))
##    #aa=plot!(sol.t, getindex.(sol.u,  [1]), legend = false)
##    #display(aa)
##end
###plot!(sol.t, getindex.(sol.u,  [1]), legend = false)
##
###tvec = range(tspan[2]-τmax, tspan[2], length(getindex.(vfix_simulation,  [1])))
###plot!(tvec, getindex.(vfix_simulation,  [1]))
##
###plot!(getindex.(vfix_simulation, [2]))
##plot(sol.t, getindex.(sol.u,  [1]))
##plot!(dpdp.StateSmaplingTime.+(nTurn+1)*T,real.(getindex.(s0aff, [1])))
##plot!()
##
###plot!(real.(getindex.(s0aff, [2])))
##
##maximum(abs.(getindex.(vfix_simulation, 1)))
##maximum(abs.(getindex.(s0aff, 1)))


# .................................................................................................................
# --------------------------------------Stab. map with brute force------------------------------------------

println("----------Start brute-force---------------")

OMrels   = range(1600.0, 2500.0, 1000) * pi / 30
ws       = range(0.1,2.0, 1000) / 1e3
Aaff     = zeros(size(ws, 1), size(OMrels, 1))
Spek_aff = zeros(size(ws, 1), size(OMrels, 1))
p2p      = zeros(size(ws, 1), size(OMrels, 1))
p2p_stab = zeros(size(ws, 1), size(OMrels, 1))

function fooTurningSSV(OMrel, w)
    τmax = 2pi / OMrel / (1.0 - RVA)
    T = 2pi / OMrel / RVF

    Nstep= Int(ceil(50 * τmax /  min(τmax, T, Tn)))
    dpdploc = dynamic_problemSampled(probTurningSSV, MethodOfSteps(BS3()), τmax, T; Historyresolution=Nstep, eigN=1, zerofixpont=false)

    #ABSmuMax = spectralradius(dpdploc;  p = (ζ, ωn, k, w, K, OMrel, RVA, RVF, f0));
    muaff, s0aff = affine(dpdploc; p = [ζ, ωn, k, w, K, OMrel, RVA, RVF, f0]);
    peak_to_peak = maximum(real.(getindex.(s0aff,  [1]))) - minimum(real.(getindex.(s0aff,  [1])))

    #return ABSmuMax-1.0, peak_to_peak
    return abs(muaff[1][1])-1.0, peak_to_peak
end

@show Threads.nthreads()
@time Threads.@threads for i in 1:size(OMrels,1)
    println(i)
    for j in 1:size(ws,1)

        Spek_aff[j,i], p2p[j,i] = fooTurningSSV(OMrels[i], ws[j])

        if Spek_aff[j,i] < 0
            p2p_stab[j,i] = p2p[j,i]
        else
            p2p_stab[j,i] = NaN
        end
    end
end

heatmap(OMrels, ws, p2p_stab, c=reverse(cgrad(:hot)), clims=(0, 10.0^-4.5), colorbar=true)
heatmap(OMrels, ws, p2p_stab .* 1e6, c=reverse(cgrad(:hot)), clims=(0, 10.0^(6-4.5)), colorbar=true)
#heatmap(OMrels, ws, p2p_stab, c=reverse(cgrad(:hot)), colorbar=true)

aa=plot!( dpi=500, legend=false,size=(560, 350),xlabel =L"spindle\ speed - \Omega_0 [rad/s]",
ylabel =L"k_w \ [N/m]")
display(aa)
RVFinv=1/RVF
savefig("Long_sim_Fig_2a_SSV_$RVFinv.png")

#contourf(OMrels, ws, p2p_stab, levels=1000, color=:hot, clims=(0, 1e-5))
#contour(OMrels, ws, Spek_aff, levels=[0], linecolor=:black, colorbar=false) 
#contour(OMrels, ws, Spek_aff, levels=[0], linecolor=:red, colorbar=false) 
#contour!(OMrels, ws, (p2p .- 1e-6), levels=[0] ,linecolor=:black)

#plot(OMrels, ws, p2p_stab, st = :surface, color=:black)
#plot!(OMrels, ws, p2p)

# .................................................................................................................
# --------------------------------------Stab. map with MDBM------------------------------------------
using MDBM

ax1 = Axis(1600:50:2500,"OM")*pi/30  # initial grid in x direction
ax2 = Axis(0.0:0.1:2.0,"w") / 1e3    # initial grid in y direction

function fooTurningSSV(OMrel, w)

    p = [ζ, ωn, k, w, K, OMrel, RVA, RVF, f0]

    τmax = 2pi / OMrel / (1.0 - RVA)
    T = 2pi / OMrel / RVF

    Nstep = Int(ceil(50 * τmax /  min(τmax, T, Tn)))
    dpdploc = dynamic_problemSampled(probTurningSSV, MethodOfSteps(BS3()), τmax, T; Historyresolution=Nstep, eigN=1, zerofixpont=false)

    ABSmuMax = spectralradius(dpdploc;  p = p);

    return ABSmuMax-1.0
end
mymdbm=MDBM_Problem(fooTurningSSV,[ax1,ax2])
iteration=3#number of refinements (resolution doubling)
@time MDBM.solve!(mymdbm,iteration);
#points where the function foo was evaluated
x_eval,y_eval=getevaluatedpoints(mymdbm);
#interpolated points of the solution (approximately where foo(x,y) == 0 and c(x,y)>0)
x_sol,y_sol=getinterpolatedsolution(mymdbm);
#scatter(x_eval,y_eval,markersize=1)
#scatter!(x_sol,y_sol,markersize=3)
scatter!(x_sol,y_sol,markersize=1, color=:black, xlims=(ax1[1], ax1[end]), ylims=(ax2[1], ax2[end]))





matwrite("fixpoints.mat", Dict("fixpoints" => getindex.(vfix_simulation, [1])))
matwrite("muaff.mat", Dict("muaff" => muaff[1]))

matwrite("stabmapvals.mat", Dict(
  "Om" => convert(Vector{Float64}, OMrels), 
  "w"  =>  convert(Vector{Float64}, ws),
  "mu" => Spek_aff
))

matwrite("fixpointtest.mat", Dict(
  "RVA" =>  RVA,
  "RVF" => RVF,
  "Om" =>  OMrel, 
  "w"  =>   w,  
  "quasi_t" => convert(Vector{Float64}, tvec1.- (T- τmax)), 
  "sim_t" => convert(Vector{Float64}, tvec3), 
  "quasi_FP" => convert(Vector{Float64}, quasistaticfixpoint.(tvec1)), 
  "sim_FP"  =>  convert(Vector{Float64}, real.(getindex.(s0aff, [1]))),
))
