#TODO: not working - some part from ChatGPT for testing
5+5

using Revise
using DDE_mapping

using BenchmarkTools
using Plots
#plotly()
using Profile
using StaticArrays
using DifferentialEquations

using LaTeXStrings

using MDBM

function van_der_pol!(du, u, p, t)
    x1, x2 = u
    mu = p
    du[1] = x2
    du[2] = mu * (1 - x1^2) * x2 - x1
end


Base.:+(a::SVector, b::Bool) = a .+ b #TODO: where to put this?
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?
Base.:+(b::Float64, a::SVector) = a .+ b #TODO: where to put this?


##<<<<<<<<<<< Lin Map based on


u0 = [1.0, 0.0]      # Initial condition
tspan = (0.0, 50.0)  # Time span
mu = 0.1             # Bifurcation parameter

T=1
#p = (ζ, ωn, k, τ,10.0)

probVdP = ODEProblem(van_der_pol!, u0, tspan, mu)

Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-4)#
solVdP = solve(probVdP, Tsit5())  # Using Tsitouras 5/4 Runge-Kutta method

solVdP = solve(probVdP; Solver_args...)#abstol,reltol
plot(sol)
using Plots
plot(solVdP, vars=(1, 2), title="Van der Pol Oscillator Phase Portrait", xlabel="x1", ylabel="x2")



mus = [0.1, 0.5, 1.0]  # Different values of mu to observe bifurcation
for mu in mus
    probVdP = ODEProblem(van_der_pol!, u0, tspan, mu)
    solVdP = solve(probVdP, Tsit5())
    aa=plot(solVdP, vars=(1, 2), title="Van der Pol Oscillator Phase Portrait (mu=$mu)", xlabel="x1", ylabel="x2")
    display(aa)
end









tFIXend = (-τ:0.01:0)
uFIXend = sol(sol.t[end] .+ tFIXend).u
plot(tFIXend, getindex.(uFIXend, 1))
plot!(tFIXend, getindex.(uFIXend, 2))


# fix point and spectrum test-----------------------------
Nstep = 150
τmax = 2pi + 0.5
Krylov_arg
dpdp = dynamic_problemSampled(probMathieu, Solver_args, τmax,
    T; Historyresolution=Nstep, eigN=10, zerofixpont=false);

    @benchmark affine($dpdp; p=$p)

muaff, s0aff = affine(dpdp; p=p);
muaff, s0aff = affine(dpdp, s0aff; p=p);
plot(log.(abs.(muaff[1])))

scatter(muaff[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))


plot(tFIXend, getindex.(uFIXend, 1))
plot!(tFIXend, getindex.(uFIXend, 2))
plot!(LinRange(-τmax, 0.0, size(s0aff, 1)), getindex.(s0aff, 1))
plot!(LinRange(-τmax, 0.0, size(s0aff, 1)), getindex.(s0aff, 2))


# fix point and spectrum test---------------------------




#--------------- bifurcation test - "continuation" -------------


using KrylovKit
Neig=10#number of required eigen values
Krylov_arg=(Neig,:LM, KrylovKit.Arnoldi(tol=1e-12,verbosity=0));


Nstep = 50
τmax = 2pi 
dpdp = dpdp = dynamic_problemSampled(probMathieu, Solver_args, τmax,
T; Historyresolution=Nstep, zerofixpont=false, affineinteration=2,Krylov_arg=Krylov_arg);

Av = vcat(collect(0.0:0.0355:3.5), collect(3.5:0.01:3.70)) # initial grid in x direction

Aaff = zeros(size(Av, 1))
Spek_aff = zeros(size(Av, 1))


muaff, s0aff, solPeriod = affine(dpdp; p=(ζ, δ, ϵ, b, τ, T, 0.0, 1.0));


#Threads.@threads
scatter([], [])
doanimation = false#true#false

a = Animation()
@time for j in 1:size(Av, 1)
    #@inbounds 

    #println(j / size(Av, 1))
    A = Av[j]
    #muaff, s0aff = affine(dpdp; p=(ζ, δ, ϵ, b, τ, T,A,μ))
    muaff, s0aff, solPeriod = affine(dpdp, s0aff; p=(ζ, δ, ϵ, b, τ, T, A, μ))
    #Aaff[j] = norm(getindex.(s0aff, 1))
    #Aaff[j] = sum(getindex.(s0aff, 1).^2)/Nstep
    Aaff[j] = maximum(getindex.(s0aff, 1))
    # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1)))
    Spek_aff[j] = maximum(abs.(muaff[1]))

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
plot(Av, Aaffsat,lab="")
plot!(Av, Spek_aff,lab="")
scatter!(Av, Aaffsat, zcolor=3.0 .- sign.(Spek_aff .- 1.0), ylim=(0.0, 2.0),lab="",
xlabel=L"A - excitation level", ylabel=L"x_{max}.  / \mu_1")