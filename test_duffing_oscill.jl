#TODO: not working -  bifurcation?? problems due to the zero delay!!!
using Revise
using DDE_mapping
5 + 5
using BenchmarkTools
using Plots
plotly()
using Profile
using StaticArrays
using DifferentialEquations

using LaTeXStrings

using MDBM

function Duffing(u, h, p, t)
    # Parameters
    ϵ, γ, α, ω = p
    # Components of the delayed differential equation
    dx = u[2]
    ddx = -ϵ * u[2] - γ * u[1] - α * u[1]^3 + cos(ω * t)    # Update the derivative vector
   # ddx = -ϵ * u[2] - γ * u[1] - α * h(p,t-0.3)[1]^3 + cos(ω * t)    # Update the derivative vector
    
    #ddx = -ϵ * u[2] - γ * h(p,t-0.1)[1] - α * u[1]^3 + cos(ω * t)    # Update the derivative vector
    SA[dx, ddx]
end
Base.:+(a::SVector, b::Bool) = a .+ b #TODO: where to put this?
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?
Base.:+(b::Float64, a::SVector) = a .+ b #TODO: where to put this?


##<<<<<<<<<<< Lin Map based on
ϵ =0.001;0.05#0.09
ϵ =0.2#0.09
γ = 1.0
#α = 1.0
α = 0.0002#0.001
α = 0.01#0.001
ω = 1.116583-0.09
#ω = 0.502

p = [ϵ, γ, α, ω]

Tper = 2pi / ω
τmax = 0.1
u0 = SA[0.0, 0.0]
h(p, t) = SA[0.0; 0.0]
probDuffing = DDEProblem(Duffing, u0, h, (0.0, Tper * 1000.0), p; constant_lags=[τmax])

Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-12)#
#Solver_args = Dict( :verbose => false, :reltol => 1e-7)#

sol = solve(probDuffing; Solver_args...)#abstol,reltol
plot(sol)



tFIXend = (-Tper:0.01:0)
uFIXend = sol(sol.t[end] .+ tFIXend).u
plot(tFIXend, getindex.(uFIXend, 1))
plot!(tFIXend, getindex.(uFIXend, 2))

plot(getindex.(uFIXend, 1), getindex.(uFIXend, 2))

# fix point and spectrum test-----------------------------
# ---------------- Affine mapping ---------------------------
using KrylovKit
Neig = 2#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=Neig + 5, verbosity=0));

#Neig =3#number of required eigen values
#Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi( verbosity=0));


Nstep = 20 # discretization number of the mapping

#Creating the problem
dpDuffing = dynamic_problemSampled(probDuffing, Solver_args, τmax,
    Tper; Historyresolution=Nstep,
    zerofixpont=false, affineinteration=3,
    Krylov_arg=Krylov_arg);


#solvig the problem (mu: Floquet multiplier, saff: discretized foxed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
@time mu, s0aff, sol0 = affine(dpDuffing; p=p);

#foo(λ)=affine(dpDuffing; p=[ϵ, λ, α, ω])[2]
#
#e=0.0001
#(foo(0.05).-foo(0.05-e)) ./ e
#using ForwardDiff
#ForwardDiff.derivative(foo, 0.05)
#
#
#
#s0 = [0.0* dpDuffing.Problem.u0 for _ in 1: Nstep]
#v0 = LinMap(dpDuffing, s0; p=p)[1]
#foo(x)=LinMap(dpDuffing, x; p=p)[1]
#
#
#ForwardDiff.jacobian(foo, s0)
#function Base.one(::Type{SVector{2, Float64}}) 
#    SA[ones(Float64,2)...]
#end
#ForwardDiff.can_dual(::Type{SVector{2, Float64}}) = true
#one(typeof(SA[1.0,2.0]))


# Comparing the solutions:
plot(getindex.(uFIXend, 1), getindex.(uFIXend, 2))
#plot!(getindex.(saff,1), getindex.(saff,2), marker=:circle,markersize=4,lab="")#
plot!(sol0[1, :], sol0[2, :], marker=:circle, lw=0, lab="")#marker=:cross,markersize=2)#
##


# Plotting the Floquet multipliers
#in log scale
plot(log.(abs.(mu[1])))
#in complex plane
scatter(mu[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))




# fix point and spectrum test---------------------------

using ForwardDiff
one_espilon_Dual = ForwardDiff.Dual{Float64}(0.0, 1.0)
p_dual = [ϵ, γ, α, ω+one_espilon_Dual]


Δv0=(LinMap(dpDuffing, s0aff; p=[ϵ, γ, α, ω+1e-5])[1]-LinMap(dpDuffing, s0aff; p=[ϵ, γ, α, ω])[1])/1e-5


dv0=LinMap(dpDuffing, s0aff; p=[ϵ, γ, α, ω+one_espilon_Dual])[1]
partialpart.(dv0)
## --------------- bifurcation test - "continuation" -------------
using KrylovKit
Neig =2#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-20,krylovdim=Neig+2; verbosity=0));



ωv = vcat(LinRange(0.2, 1.8, 378)) # initial grid in x direction

#ωv = vcat(LinRange(1.8, 0.2, 878)) # initial grid in x direction


Aaff = zeros(size(ωv, 1))
Spek_aff = zeros(size(ωv, 1))

ω=ωv[1]
mu, s0aff, sol0 = affine(dpDuffing; p=[ϵ, γ, α, ω]);
mu, s1aff, sol0 = affine(dpDuffing, s0aff; p=[ϵ, γ, α, ω]);



#Threads.@threads
scatter([], [])
doanimation = false#true#false

a = Animation()
@time for j in 1:size(ωv, 1)
    
    ω = ωv[j]
    Tper = 2pi / ω

    println(ω)
    
    Nstep=10
    dpDuffing = dynamic_problemSampled(probDuffing, Solver_args, τmax,
        Tper; Historyresolution=Nstep,
        zerofixpont=false, affineinteration=2,
        Krylov_arg=Krylov_arg)
            muaff, s0aff, solPeriod = affine(dpDuffing; p=[ϵ, γ, α, ω])
       #muaff, s0aff, solPeriod = affine(dpDuffing, s0aff; p=[ϵ, γ, α, ω])
    #Aaff[j] = norm(getindex.(s0aff, 1))
    #Aaff[j] = sum(getindex.(s0aff, 1).^2)/Nstep
    Aaff[j] = maximum(getindex.(solPeriod, 1))
    # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1)))
    Spek_aff[j] = maximum(abs.(muaff[1]))
    Spek_aff[j] = log(maximum(abs.(muaff[1])))/Tper

    # if doanimation
    #     #scatter(log.(muaff[1]),lab="")
    #     scatter(muaff[1], lab="")
    #     #plt=plot!(sin.(0:0.01:2pi),cos.(0:0.01:2pi),lab="")
    #     plt = plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi), xlim=(-3.0, 1.5), ylim=(-1.5, 1.5), lab="")
    #     frame(a, plt)# Not working with : #plotly()
    # else
    #     scatter!(muaff[1], lab="")
    # end
end
#scatter!([], [])
#plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))
#
#if doanimation
#    gif(a, fps=25)
#end



#@code_warntype affine(dpdp; p=(ζ, δ, ϵ, b, τ, T,A,μ))
theme(:ggplot2)

Aaffsat = deepcopy(Aaff);
#Aaffsat[Aaffsat .== maximum(Aaffsat)] .= 0.0
#Aaffsat[Spek_aff .> 1.0] .= 0.0;
plot(ωv, Aaffsat, lab="")
plot!(ωv, Spek_aff, lab="")
scatter(ωv, Aaffsat, zcolor=3.0 .- sign.(Spek_aff .- 1.0), ylim=(0.0,522.0), lab="",
    xlabel=L"\omega - excitation level", ylabel=L"x_{max}.  / \mu_1")

#scatter(ωv, Aaffsat, zcolor=3.0 .- sign.(Spek_aff .- 1.0), lab="", xlabel=L"\omega - excitation level", ylabel=L"x_{max}.  / \mu_1")


## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## --------- pseudo-archlength continuation -------------

function add_a_point!(Solset, requiredstepsize, parindex=4)

    pall = deepcopy(Solset[end][3].prob.p)

    p2 = Solset[end][3].prob.p[parindex]
    p1 = Solset[end-1][3].prob.p[parindex]
    p0 = Solset[end-2][3].prob.p[parindex]
    pghange = (p2 - p1)

    #muaff, s0aff, solPeriod = affine(dpDuffing, s0aff; p=(ϵ, γ, α, ω))
    Mscaler = 1 / 100/length(Solset[end][2])
    XnormChange=(Solset[end][2]-Solset[end-1][2])' *(Solset[end][2]-Solset[end-1][2])
    #Xnorm2 = Solset[end][2]' * Solset[end][2]
    #Xnorm1 = Solset[end-1][2]' * Solset[end-1][2]
    #XnormChange = (Xnorm2 - Xnorm1)


    perviouse_stepsize = sqrt((XnormChange * Mscaler)^2 + pghange^2)

    scaling = abs( requiredstepsize / perviouse_stepsize)

    pnew = p2 + scaling * pghange
    pall[parindex] = pnew

    omega   =pnew
    Tper = 2pi / omega
    dpDuffing_loc = dynamic_problemSampled(probDuffing, Solver_args, τmax,
        Tper; Historyresolution=Nstep,
        zerofixpont=false, affineinteration=2,
        Krylov_arg=Krylov_arg)
        
   # S0aff_new= 0.0 * Solset[end][2] 
   # S0aff_new=Solset[end][2] 
   # S0aff_new = (1.0 + scaling) * Solset[end][2] - scaling * Solset[end-1][2]
    S0aff_new = (1.0 + scaling) * Solset[end][2] - scaling * Solset[end-1][2]
   #pnew = 2.0*p2 -1.5*p1+0.5*p0
   #S0aff_new =2.0* Solset[end][2] - 1.5* Solset[end-1][2]+0.5* Solset[end-2][2]
# S0aff_new = (1.0 + scaling) * Solset[end][2] - scaling * Solset[end][2]
    push!(Solset, affine(dpDuffing_loc, S0aff_new; p=pall))
end
    
        

using KrylovKit
Neig = 2#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-25,krylovdim=4,  verbosity=0));



requiredstepsize = 0.025
p_index = 4
println("NINCS  gerjesztés és a periodus idő a rendszerhez állítva!!!!!!!4- HIBÁS")
Sols = [affine(dpDuffing; p=[ϵ, γ, α, 0.2]), affine(dpDuffing; p=[ϵ, γ, α, 0.2 + requiredstepsize])];
Sols = [affine(dpDuffing; p=[ϵ, γ, α, 1.5]), affine(dpDuffing; p=[ϵ, γ, α, 1.5 - requiredstepsize])];
Sols = [ affine(dpDuffing; p=[ϵ, γ, α, 0.2 + k* requiredstepsize]) for k in 0:2];

#requiredstepsize = -0.025
#Sols = [ affine(dpDuffing; p=[ϵ, γ, α, 1.8 + k* requiredstepsize]) for k in 0:5];

@time for jj in 1:500
    length(Sols)
   # println(jj)
   println(Sols[end][3].prob.p[4])
    add_a_point!(Sols, requiredstepsize, p_index);
end
aa = scatter!([Sol[3].prob.p[4] for Sol in Sols], [maximum(getindex.(Sol[3].u, 1)) for Sol in Sols], lab="",ylim=(0.0,102.0))
#aa = scatter!([Sol[3].prob.p[4] for Sol in Sols], [maximum(getindex.(Sol[3].u, 1)) for Sol in Sols], lab="")
#plot!(ωv, Spek_aff, lab="")
#aa = scatter([maximum(getindex.(Sol[3].u, 1)) for Sol in Sols], lab="",ylim=(0.0,102.0))
#plot!(ωv, Spek_aff, lab="")
aa = scatter!([Sol[3].prob.p[4] for Sol in Sols], [maximum(getindex.(Sol[3].u, 1)) for Sol in Sols], lab="")

aa = scatter!([Sol[3].prob.p[4] for Sol in Sols], [(imag(Sol[1][1][1])) for Sol in Sols], lab="")
aa = scatter!([Sol[3].prob.p[4] for Sol in Sols], [(real(Sol[1][1][1])) for Sol in Sols], lab="")
aa = scatter!([Sol[3].prob.p[4] for Sol in Sols], [(abs(Sol[1][1][1])) for Sol in Sols], lab="")

display(aa)

#
bb = scatter!([Sol[3].prob.p[4] for Sol in Sols], 
             [real(log(Sol[1][1][1])/(Sol[3].prob.tspan[end])) for Sol in Sols], lab="")
#bb = scatter([real(log(Sol[1][1][1])/(Sol[3].prob.tspan[end])) for Sol in Sols], [imag(log(Sol[1][1][1])/(Sol[3].prob.tspan[end])) for Sol in Sols], lab="")
display(bb)
