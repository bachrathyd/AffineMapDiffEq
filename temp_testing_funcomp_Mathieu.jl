

#------------------------ Mathieu test------------------
using Revise
#includet("src\\DDE_mapping.jl")
includet("src\\DDE_mapping_types.jl")
includet("src\\DDE_mapping_functions.jl")
includet("src\\fucntion_composition.jl")

using BenchmarkTools
using Plots
plotly()
using Profile
using StaticArrays
using DifferentialEquations

# using Memoization

using MDBM
#TODO:  @fastmath

function f_now(p, t)
    ζ, δ, ϵ, b, τ, T = p
    # println("Computed $t")
    # println("Computed $p")
    SA[-(δ + ϵ * cos(t)), -2*ζ]
    #SA[-(δ+ϵ*sign(cos(t))), -2*ζ]
end
#@memoize 
function f_past(p, t)
    ζ, δ, ϵ, b, τ, T = p
    SA[b, 0]
end

function DelayMathieu(u, h, p, t)
    ζ, δ, ϵ, b, τ, ν, T = p
    dx = u[2]
    #ddx = [-(δ+ϵ*cos(t)), -2*ζ]' * u + b * h(p, t - τ)[1] # Traditional Delayed Mathieu equation
    #sign
    ddx = [-(δ + ϵ * (cos(t))), -2 * ζ]' * u + b * h(p, t - τ)[1]# + ν * h(p, t - τ, Val{1})[2]#The last part is a test for neutral system
    SA[dx, ddx]
end


ζ = 0.03
δ = 1.1
ϵ = 0.2
τ = 2pi
b = 0.2
ν = 0.25 # Neutral coefficient
T = 16;
2pi;
p = ζ, δ, ϵ, b, τ, ν, T
#p = (ζ, ωn, k, τ,10.0)

u0 = SA[1.0, 1.0]

h(p, t::Float64) = SA[1.0; -1.0]
probMathieu = DDEProblem(DelayMathieu, h, (0.0, T * 1.0), p; constant_lags=[τ], neutral=true)

sol = solve(probMathieu, MethodOfSteps(BS3()),abstol=1e-10,reltol=1e-10)
plot(sol)

#----------Sampled Affine Map -----------

Base.:+(a::SVector, b::Bool) = a .+ b
Nstep = 50
τmax = 2pi + 0.1
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax, T; Historyresolution=Nstep, eigN=20, zerofixpont=true);

@time muaff, s0aff = affine(dpdp; p=p);
#benchmark affine(dpdp; p=p)
muaff[1]



#---------------------------Krylov-Functional-------------

his2 = funcomp(t->h(p,t), [-τmax, 0.0])

#Base.return_types(h,(Any,Float64))
#Base.return_types(his2.f[1],(Float64,))
#Base.return_types(his2,(Float64,))
#@benchmark h([], -1.0)
#@benchmark his2( -1.0)
#@code_warntype his2( -1.0)

function LinMap(his::funcomp)::funcomp# where Tout  #{T,Tout}
    println(his)
    sol = solve(remake(probMathieu; u0=his(0.0), tspan=(0.0, T), h=(p,t,kwargs...)->his(t), p=p), MethodOfSteps(BS3()),abstol=1e-10,reltol=1e-10; verbose=false);
        function f(t::Float64)#::Tout
        getvalues(sol,(t + T))
    end
    funcomp([1.0], [f], his.range);
end

plot( solve(remake(probMathieu; u0=h(p, 0.0), tspan=(0.0, T*5), h=h, p=p), MethodOfSteps(RK4()),abstol=1e-10,reltol=1e-10; verbose=false) ,width=3)
#plot!( solve(remake(probMathieu; u0=his2( 0.0), tspan=(0.0, T), h=(p,t,kwargs...)->his2(t), p=p), MethodOfSteps(BS3()); verbose=false) )
#plot!( solve(remake(probMathieu; u0=his2( 0.0), tspan=(0.0, T), h=(p,t)->his2(t), p=p), MethodOfSteps(BS3()); verbose=false) )
@time out=LinMap(his2)
kk=1
t=out.range[1]:0.01:out.range[2]
scatter!(t.+T*kk,getindex.([out(tt) for tt in t],1),xlims = (-τmax+T,T*kk))

kk +=1
@time out=LinMap(out)
plot!(t.+T*kk,getindex.([out(tt) for tt in t],1),xlims = (-τmax+T,T*kk) ,width=2,linestyle=:dash)
plot!(t.+T*kk,getindex.([out(tt) for tt in t],2),xlims = (-τmax+T,T*kk) ,width=2,linestyle=:dash)

#@benchmark solve(remake(probMathieu; u0=h(p, 0.0), tspan=(0.0, T), h=h, p=p), MethodOfSteps(BS3()); verbose=false)
#@benchmark solve(remake(probMathieu; u0=his2( 0.0), tspan=(0.0, T), h=(p,t,kwargs...)->his2(t), p=p), MethodOfSteps(BS3()); verbose=false)
#@benchmark LinMap(his2)

#sol = solve(probMathieu, MethodOfSteps(BS3()))#abstol,reltol
#plot(sol)
#t=(-3.0:0.5:63)
#scatter!(t,getindex.(sol.(t),1))
#t=(-3.0:0.1:00)
#scatter!(t,getindex.([h(p, tt) for tt in t],1),xlims = (-10,70),)
#t=(-3.0:0.52:63)
#scatter!(t,getindex.([getvalues(sol,tt) for tt in t],1),xlims = (-10,70),)

using Profile

@time musFun = getindex(schursolve(LinMap, his2, 1, :LM, KrylovKit.Arnoldi(krylovdim=5,tol=1e-4, verbosity=3)), [3, 2, 1,4]);

@time musFun = getindex(schursolve(LinMap, his2, 1, :LM, KrylovKit.Arnoldi(krylovdim=5,tol=1e-7, verbosity=3)), [3, 2, 1,4]);
#reduce!(musFun[2][1])


plot(log.(abs.(muaff[1])))
plot!(log.(abs.(musFun[1])))

scatter(muaff[1])
scatter!(musFun[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))


NN=1#1:2:7
plot(dpdp.StateSmaplingTime, getindex.(muaff[2][NN],1))
plot!(dpdp.StateSmaplingTime, getindex.(muaff[2][NN],2))
plot!(t,getindex.([musFun[2][NN+1](tt) for tt in t],1) ,width=2,linestyle=:dash)
plot!(t,getindex.([musFun[2][NN+1](tt) for tt in t],2) ,width=2,linestyle=:dash,title=abs.([muaff[1][NN];musFun[1][NN]]))

norm(musFun[2][NN])
norm(muaff[2][NN])
#error("Deriváltal baj van Neutrális esetben!!!!")
#h=(p,t,kwargs...)->his2(t)
## Different behavioure
T = 0.3
T = 6.0 # somehow the number of scalers and the number of function are not the same.
#TODO: question: what happens if the same function is used by differnet instance, and on is chaneg inplace. 
#especiall, it is racursive useage??!?! -> brute-forece solution: try deepcopy of all functions after it is provided.
T = 10.0
