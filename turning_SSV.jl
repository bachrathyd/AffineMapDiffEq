5 + 5
using Revise
#includet("src\\DDE_mapping.jl")
includet("src\\DDE_mapping_types.jl")
includet("src\\DDE_mapping_functions.jl")

using StaticArrays
using DifferentialEquations
using Plots

#using BenchmarkTools
#using Profile
#using MDBM
function delayed_turning_STATIC_SSV(u, h, p, t)
    # Parameters
    ζ, ωn, k, OMrel, RVA, RVF, f0 = p
    T0 = 2pi / OMrel

    #tau = T0 .* (1.0 .- RVA .* cos.(t .* (OMrel * RVF)))
    tau = T0 * (1.0 - RVA * cos(t * (OMrel * RVF)))
    #println(tau)
    f0 = f0 * tau / T
    h0 = f0

    #if mod(t/T,2.0*pi)<0.4
    #    h0=f0*sin(2.0*pi*t/T)
    #else
    #    h0=0.0;
    #end


    dx = u[2]
    ddx = -2 * ζ * u[2] - ωn * (u[1]) + k * (h0 + u[1] - h(p, t - tau)[1])  # Surface regeneration effect

    # Update the derivative vector
    SA[dx, ddx]
end

Base.:+(a::SVector, b::Bool) = a .+ b


##<<<<<<<<<<< Lin Map based on
u0 = SA[0.0, 0.0]

ζ = 0.2          # damping coefficient
ωn = 1.0#0.2          # nat. freq
k = 0.15#4#5#8;#5         # cut.coeff

OMrel = 0.7;
RVA = 0.1;
RVF = 1 / 3;
T = 2pi / (OMrel * RVF)

f0 = 1.0         # excitation


p = [ζ, ωn, k, OMrel, RVA, RVF, f0]


tspan = (0.0, T * 1400.0)

h(p, t) = SA[0.0; 0.0]
probTurningSSV = DDEProblem(delayed_turning_STATIC_SSV, u0, h, tspan, p)

#StateSmaplingTime = LinRange(-τ, 0.0, 200) .+ T
sol = solve(probTurningSSV, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)


function longerm_sim_fix_pint(dp, p)
    #solve(prob,RK4();dt=0.005,adaptive=false)


    OneMap = s -> LinMap(dp, s; p=p)


    #Fix point
    vfix_simulation = zeros(typeof(u0), Nstep)#v0;s0
    vfix_simulation = rand(typeof(u0), Nstep)#v0;s0
    for _ in 1:100
        vfix_simulation = OneMap(vfix_simulation)
    end
    norm(vfix_simulation - LinMap(dp, vfix_simulation; p=p))

    return vfix_simulation

end


Nstep = 150
τmax = 2pi / OMrel * (1.0 + RVA) + 0.1
dpdp = dynamic_problemSampled(probTurningSSV, MethodOfSteps(BS3()), τmax, T; Historyresolution=Nstep, eigN=5, zerofixpont=false);


vfix_simulation = longerm_sim_fix_pint(dpdp, dpdp.DDEdynProblem.p)





vfix_simulation = longerm_sim_fix_pint(dpdp, [ζ, ωn, k, OMrel, RVA, RVF, f0])
muaff, s0aff = affine(dpdp; p=[ζ, ωn, k, OMrel, RVA, RVF, f0]);

plot(getindex.(vfix_simulation, [1]))
plot!(getindex.(vfix_simulation, [2]))
plot(real.(getindex.(s0aff, [1])))
plot!(real.(getindex.(s0aff, [2])))
maximum(abs.(getindex.(vfix_simulation, 1)))
maximum(abs.(getindex.(s0aff, 1)))



println("----------Start brute-force---------------")
OMv = 1.2:-0.01:0.2
kv = -1.0:0.01:1.0
Aaff = zeros(size(kv, 1), size(OMv, 1))
Spek_aff = zeros(size(kv, 1), size(OMv, 1))

RVA = 0.0;
RVF = 1.0;
#Threads.@threads 
@time for j in 1:size(OMv, 1)
    OMrel = OMv[j]
    Nstep = 300
    τmax = 2pi / OMrel * (1.0 + RVA) + 0.1
    T = 2pi / OMrel / RVF
    dpdp = dynamic_problemSampled(probTurningSSV, MethodOfSteps(BS3()), τmax, T; Historyresolution=Nstep, eigN=2, zerofixpont=false)

    Threads.@threads  for i in 1:size(kv, 1)
        # println([i,j])
        k = kv[i]
        p = (ζ, ωn, k, OMrel, RVA, RVF, f0)
        muaff, s0aff = affine(dpdp; p=(ζ, ωn, k, OMrel, RVA, RVF, 0.0))
        Aaff[i, j] = norm(getindex.(s0aff, 1))
        # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1)))
        Spek_aff[i, j] = maximum(abs.(muaff[1]))
    end
end


Aaffsat = deepcopy(Aaff);
Aaffsat[Spek_aff.>1.0] .= 0.0;
heatmap(OMv, kv, log.(Aaffsat))


Spek_affsat = deepcopy(Spek_aff);
Spek_affsat[Spek_affsat.>1.0] .= 0.0;
heatmap(OMv, kv, log.(Spek_affsat))
