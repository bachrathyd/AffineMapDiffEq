5 + 5
#using Revise
##includet("src\\DDE_mapping.jl")
#includet("src\\DDE_mapping_types.jl")
#includet("src\\DDE_mapping_functions.jl")
using DDE_mapping



##] add https://github.com/bachrathyd/AffineMapDiffEq.git
#using AffineMapDiffEq

using StaticArrays
using DifferentialEquations
using Plots
plotly()
using MAT
using LinearAlgebra
#using BenchmarkTools
#using Profile
#using MDBM



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

    dx  = u[2]
    #ddx = ωn^2 * (K*w/k * (h(p, t - tau)[1] - u[1] - h0) - u[1]) -2 * ζ * ωn * u[2] # Surface regeneration effect

    # ------------ szög -----------
    fi = t * OMrel
    if 0 < mod(fi, 2*pi) && mod(fi, 2pi) < deg2rad(10)
        ddx = ωn^2 * (- u[1]) -2 * ζ * ωn * u[2] 
    else
        ddx = ωn^2 * (K*w/k * (h(p, t - tau)[1] - u[1] - h0) - u[1]) -2 * ζ * ωn * u[2] 
    end
    # -----------------------------
    
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


OMrel = 170;

RVA   = 0.05;
RVF   = 0.05;
RVA   = 0;
RVF   = 1;

T     = 2pi / (OMrel * RVF)
w     = 0.6/1e3;
w     = 0.3/1e3;

f0 = 0.0001        # feed [m/s]

p = [ζ, ωn, k, w, K, OMrel, RVA, RVF, f0]

tspan = (0.0, T * 1.0)

h(p, t) = SA[0.0; 0.0]

#Egyetlen hpsszú szimuláció
probTurningSSV_sim = DDEProblem(delayed_turning_STATIC_SSV, u0, h, (0.0, T * 500.0), p)
sol1 = solve(probTurningSSV_sim, MethodOfSteps(BS3()),dt=T/2000*RVF,reltol=1e-18)
plot(sol1)


# --------------------------------------------------------------------------------------

probTurningSSV = DDEProblem(delayed_turning_STATIC_SSV, u0, h, tspan, p)
#StateSmaplingTime = LinRange(-τ, 0.0, 200) .+ T
sol = solve(probTurningSSV, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)
plot(sol(tspan[2] - T .+ (0:T/200:T)))

function longerm_sim_fix_point(dp, p)
    #solve(prob,RK4();dt=0.005,adaptive=false)

    OneMap = s -> LinMap(dp, s; p=p)

    #Fix point
    vfix_simulation = zeros(typeof(u0), Nstep)#v0;s0
    #vfix_simulation = rand(typeof(u0), Nstep)#v0;s0
    for ii in 1:500
        vfix_simulation = OneMap(vfix_simulation)
    end
    norm(vfix_simulation - LinMap(dp, vfix_simulation; p=p))

    return vfix_simulation
end

τmax = 2pi / OMrel / (1.0 - RVA)
Nstep= Int(ceil(500 * τmax /  min(τmax, T, Tn)))
dpdp = dynamic_problemSampled(probTurningSSV, MethodOfSteps(BS3()), τmax, T; Historyresolution=Nstep, eigN=5, zerofixpont=false);

vfix_simulation = longerm_sim_fix_point(dpdp, p)
tvec2 =  dpdp.StateSmaplingTime#range(0, τmax, length(getindex.(vfix_simulation,  [1])))
plot!(tvec2 .+ 500T, getindex.(vfix_simulation,  [1]))


# ------------- affine map --------------------
muaff, s0aff = affine(dpdp; p=[ζ, ωn, k, w, K, OMrel, RVA, RVF, f0]);

tvec3 = dpdp.StateSmaplingTime .+ 500T#range(0, τmax, length(getindex.(s0aff,  [1])))
plot!(tvec3, real.(getindex.(s0aff, [1])))
#------------------------------------Plot fixpoints----------------------------------------------
#quasistaticfixpoint = (t) -> -K * w / k * f0 / (1 + RVA *sin(RVF * OMrel * t))
#
#tvec1 = dpdp.StateSmaplingTime;#range(T-τmax, T, 100)
#plot(tvec1 .- (T- τmax), quasistaticfixpoint.(tvec1))
plot(sol1)

tvec2 =  dpdp.StateSmaplingTime#range(0, τmax, length(getindex.(vfix_simulation,  [1])))
plot!(tvec2.+500T, getindex.(vfix_simulation,  [2]))
plot!(tvec2.+500T, getindex.(vfix_simulation,  [1]))

tvec3 = dpdp.StateSmaplingTime
plot!(tvec3.+500T, real.(getindex.(s0aff, [2])))
plot!(tvec3.+500T, real.(getindex.(s0aff, [1])))
plot!(xlim=[tvec3[1],tvec3[end]].+500T)


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




@show peak_to_peak = maximum(getindex.(vfix_simulation,  [1])) - minimum(getindex.(vfix_simulation,  [1]))
@show peak_to_peak = maximum(real.(getindex.(s0aff,  [1]))) - minimum(real.(getindex.(s0aff,  [1])))
#-----------------------------------------------------------------------------------------------


scatter(1:size(muaff[1],1),log.(abs.(muaff[1])))
scatter(muaff[1])
plot!(sin.(0:0.01:2pi),cos.(0:0.01:2pi), aspect_ratio=:equal)


#-----------------Simulation of multiple rotations to test fixpoint-------------------------------

tspan = (0.0, T * 1.0)
u0 = SA[0.0, 0.0]
h(p, t) = SA[0.0; 0.0]
probTurningSSV = DDEProblem(delayed_turning_STATIC_SSV, u0, h, tspan, p)
sol = solve(probTurningSSV, MethodOfSteps(BS3()))
plot(sol.t,getindex.(sol.u,  [1]))

for i in 1:20
    h(p, t) = sol(t)
    u0 = sol.u[end,1]
    tspan = tspan .+ T
    probTurningSSV = DDEProblem(delayed_turning_STATIC_SSV, u0, h, tspan, p)
    sol = solve(probTurningSSV, MethodOfSteps(BS3()))
    aa=plot!(sol.t, getindex.(sol.u,  [1]))
    display(aa)
end
plot!(sol.t, getindex.(sol.u,  [1]), legend = false)


tvec = range(tspan[2]-τmax, tspan[2], length(getindex.(vfix_simulation,  [1])))
plot!(tvec, getindex.(vfix_simulation,  [1]))

plot!()


#plot!(getindex.(vfix_simulation, [2]))
##
plot!(dpdp.StateSmaplingTime,real.(getindex.(s0aff, [1])))
#plot!(real.(getindex.(s0aff, [2])))

maximum(abs.(getindex.(vfix_simulation, 1)))
maximum(abs.(getindex.(s0aff, 1)))





# ---------------------Mat file-ba kiiratas------------------------------------
matwrite("fixpoints.mat", Dict("fixpoints" => getindex.(vfix_simulation, [1])))
matwrite("muaff.mat", Dict("muaff" => muaff[1]))







# .................................................................................................................
# -----------------------Stab. map with brute force-----------------------------------------------------

println("----------Start brute-force---------------")

OMrels   = range(1600.0, 2500.0, 50) * pi / 30
ws       = range(0.1,2.0, 50) / 1e3
Aaff     = zeros(size(ws, 1), size(OMrels, 1))
Spek_aff = zeros(size(ws, 1), size(OMrels, 1))
p2p      = zeros(size(ws, 1), size(OMrels, 1))
p2p_stab = zeros(size(ws, 1), size(OMrels, 1))


# -----------------------------------------------------------------------------------

function fooTurningSSV(OMrel, w)
    τmax = 2pi / OMrel / (1.0 - RVA)
    T = 2pi / OMrel / RVF

    Nstep= Int(ceil(50 * τmax /  min(τmax, T, Tn)))
    dpdploc = dynamic_problemSampled(probTurningSSV, MethodOfSteps(BS3()), τmax, T; Historyresolution=Nstep, eigN=1, zerofixpont=false)

#    ABSmuMax = spectralradius(dpdploc;  p = (ζ, ωn, k, w, K, OMrel, RVA, RVF, f0));
    muaff, s0aff = affine(dpdploc; p = [ζ, ωn, k, w, K, OMrel, RVA, RVF, f0]);
    peak_to_peak = maximum(real.(getindex.(s0aff,  [1]))) - minimum(real.(getindex.(s0aff,  [1])))

    #return ABSmuMax-1.0, peak_to_peak
    return abs(muaff[1][1])-1.0, peak_to_peak
end
# ----------------------Brute force stabilitási térkép számítás-----------------------------------
#@time 
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


contourf(OMrels, ws, p2p_stab, levels=1000, color=:hot, clims=(0, 1e-5))
heatmap(OMrels, ws, p2p_stab, c=reverse(cgrad(:hot)), clims=(0, 1e-5), colorbar=true)


contour(OMrels, ws, Spek_aff, levels=[0] ,linecolor=:red, colorbar=false) 
contour(OMrels, ws, Spek_aff, levels=[0] ,linecolor=:red, colorbar=false) 
contour!(OMrels, ws, (p2p .- 1e-6), levels=[0] ,linecolor=:black) 

plot(OMrels, ws, p2p_stab, st = :surface, color=:black)
plot!(OMrels, ws, p2p)

matwrite("stabmapvals.mat", Dict(
  "Om" => convert(Vector{Float64}, OMrels), 
  "w"  =>  convert(Vector{Float64}, ws),
  "mu" => Spek_aff
))



# ------------------------------------------------------------------------------
# ..............................................................................


using MDBM

ax1 = Axis(1600:50:2500,"OM")*pi/30 # initial grid in x direction
ax2 = Axis(0.0:0.1:2.0,"w") / 1e3 # initial grid in y direction

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
scatter(x_sol,y_sol,markersize=1, color=:red, xlims=(ax1[1], ax1[end]), ylims=(ax2[1], ax2[end]))




#--------------------- ez itt valamiért kiakad ---------------------------------
function fooFixpoint(OMrel, w)

    τmax = 2pi / OMrel / (1.0 - RVA)
    T = 2pi / OMrel / RVF

    p = [ζ, ωn, k, w, K, OMrel, RVA, RVF, f0]

    dpdploc = dynamic_problemSampled(probTurningSSV, MethodOfSteps(BS3()), τmax, T; Historyresolution=100, eigN=1, zerofixpont=false)
    
    vfix_simulation = longerm_sim_fix_point(dpdploc, p)

    peak_to_peak = maximum(getindex.(vfix_simulation,  [1])) - minimum(getindex.(vfix_simulation,  [1]))
    return peak_to_peak-1e-4
end
mymdbm=MDBM_Problem(fooFixpoint,[ax1,ax2])
iteration=3#number of refinements (resolution doubling)
@time MDBM.solve!(mymdbm,iteration);
#points where the function foo was evaluated
x_eval,y_eval=getevaluatedpoints(mymdbm);
#interpolated points of the solution (approximately where foo(x,y) == 0 and c(x,y)>0)
x_sol,y_sol=getinterpolatedsolution(mymdbm);
#scatter(x_eval,y_eval,markersize=1)
#scatter!(x_sol,y_sol,markersize=3)
scatter!(x_sol,y_sol,markersize=1, color=:green)
