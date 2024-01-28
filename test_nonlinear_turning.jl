5+5
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
function delayed_turning_STATIC(u, h, p, t)
    # Parameters
    ζ, ωn, k, T, f0 = p
    τ=T
    if t/T<0.4
        h0=f0*sin(2.0*pi*t/T)
    else
        h0=0.0;
    end
   # tau=τ*(1.0-0.2*cos( 2 * pi  * 2*t/1.0) )
   tau = τ

    dx = u[2]


    # ddx = - 2 * ζ * u[2] - ωn * u[1] + k * (cos(2 * pi * t) * 0.0 + 1.0) * abs.(f0 * cos(2 * pi * t/tau) + u[1] - h(p, t - tau)[1])^0.75  # Surface regeneration effect
    ddx = -2 * ζ * u[2] - ωn *( u[1]+0.2*u[1]^2) + k * abs.(h0 + u[1] - h(p, t - tau)[1])^0.75  # Surface regeneration effect
    #ddx = -2 * ζ * u[2] - ωn * u[1] + k * ( u[1] - h(p, t - τ)[1])  # Surface regeneration effect

    # Update the derivative vector
    SA[dx, ddx]
end

Base.:+(a::SVector, b::Bool) = a .+ b


##<<<<<<<<<<< Lin Map based on
u0 = SA[0.0, 0.0]
ζ = 0.2          # damping coefficient
ωn = 2.0#0.2          # nat. freq
k = 0.45#4#5#8;#5         # cut.coeff
τ = 2.4          # Time delay
f0 = 1.0         # excitation

p = [ζ, ωn, k, T, f0]
p0 = [ζ, ωn, k, T, 0.0]


tspan = (0.0, T * 400.0)

h(p, t) = SA[0.0; 0.0]
probTurning = DDEProblem(delayed_turning_STATIC, u0, h, tspan, p; constant_lags=[τ])

#StateSmaplingTime = LinRange(-τ, 0.0, 200) .+ T
sol = solve(probTurning, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)
plot(sol(sol.t[end] .- (0.0:0.01:τ*1.0)))

function longerm_sim_fix_pint(dp,p)
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
τmax = 3.0
dpdp = dynamic_problemSampled(probTurning, MethodOfSteps(BS3()), τmax, T; Historyresolution=Nstep, eigN=5, zerofixpont=true);


vfix_simulation=longerm_sim_fix_pint(dpdp,dpdp.DDEdynProblem.p)



tauv=4:0.5:28.0
kv=-1.5:0.1:1.5

#tauv=4:1.5:28.0
#kv=-1.5:0.3:1.5

Asim=zeros(size(kv,1),size(tauv,1))
Aaff=zeros(size(kv,1),size(tauv,1))
Spek_aff=zeros(size(kv,1),size(tauv,1))

vfix_simulation=longerm_sim_fix_pint(dpdp,[ζ, ωn, k, τ, f0])
muaff,s0aff=affine(dpdp; p= [ζ, ωn, k, τ, f0]);

plot(getindex.(vfix_simulation,[1]))
plot!(getindex.(vfix_simulation,[2]))
plot(real.(getindex.(s0aff,[1])))
plot!(real.(getindex.(s0aff,[2])))
maximum(abs.(getindex.(vfix_simulation,1)))
maximum(abs.(getindex.(s0aff,1)))


#using Profile
#@profview
for (j, tau) in enumerate(tauv)
    println(j/size(tauv,1))
    τmax = tau*1.3;
    Nstep =  floor(Int, τmax*2.0+5)
    
    T = tau

    dpdp = dynamic_problemSampled(probTurning, MethodOfSteps(BS3()), τmax, T; Historyresolution=Nstep, eigN=4, zerofixpont=true);
    τ=tau
    for (i,k) in enumerate(kv)
    #vfix_simulation=longerm_sim_fix_pint(dpdp,[ζ, ωn, k, τ, f0])
    #Asim[i,j]= maximum(abs.(getindex.(vfix_simulation,1)))
        
        muaff,s0aff=affine(dpdp; p= [ζ, ωn, k, τ, f0]);
       # muaff,s0aff=affine(dpdp,s0aff; p= [ζ, ωn, k, τ, f0]);
       # muaff,s0aff=affine(dpdp,s0aff; p= [ζ, ωn, k, τ, f0]);
        Aaff[i,j]= maximum(abs.(getindex.(s0aff,1)))
        Spek_aff[i,j]= maximum(abs.(muaff[1]))
    end
end


plot(Aaffsat[:,15])
#plot!(Asim[:,15])
plot!(Aaffsat[:,5])
#plot!(Asim[:,5])

#Asim[Asim .> 0.99 ] .= 00.0;
#heatmap(tauv,kv,log.(Asim))

Aaffsat=deepcopy(Aaff);
Aaffsat[Spek_aff .> 1.0] .= 0.0;
heatmap(tauv,kv,log.(Aaffsat))



Spek_affsat=deepcopy(Spek_aff);
Spek_affsat[Spek_affsat .> 1.0] .= 0.0;
heatmap(tauv,kv,(Spek_affsat))

