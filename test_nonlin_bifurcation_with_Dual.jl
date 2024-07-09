10 + 10

println(Threads.nthreads())


using Revise
using DDE_mapping

using BenchmarkTools
using Plots
#plotly()
using Profile
using StaticArrays
using DifferentialEquations

using LaTeXStrings

##TODO: Bifurcation Analysis in Julia

#using Memoization

using MDBM

#@memoize
function f_now(p, t)
    ζ, δ, ϵ, b, τ, T, A = p
    # println("Computed $t")
    # println("Computed $p")
    SA[-(δ + ϵ * cos(2pi * t / T)), -2*ζ]
end
#@memoize 
function f_past(p, t)
    ζ, δ, ϵ, b, τ, T, A = p
    SA[b, 0]
end
#f_now(p, 0.6)' * SA[1.0, 1.0]

function DelayMathieu(u, h, p, t)
    # Parameters
    ζ, δ, ϵ, b, τ, T, A, μ = p
    # Components of the delayed differential equation
    dx = u[2]
    ddx = f_now(p, t)' * u + μ * u[1]^3 + b * h(p, t - τ)[1] + A * (cos(2pi * t / T) .^ 10)  # Surface regeneration effect
    # Update the derivative vector
    SA[dx, ddx]
end
Base.:+(a::SVector, b::Bool) = a .+ b #TODO: where to put this?
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?
Base.:+(b::Float64, a::SVector) = a .+ b #TODO: where to put this?


##<<<<<<<<<<< Lin Map based on

ζ = 0.08          # damping coefficient
δ = 4.05#0.2          # nat. freq
ϵ = 0.15#4#5#8;#5         # cut.coeff
τ = 2pi         # Time delay
b = 0.5
T = 2pi
A = 0.01#3.3#0.01;
μ = 1.0
p = ζ, δ, ϵ, b, τ, T, A, μ
#p = (ζ, ωn, k, τ,10.0)

# test simulation
u0 = SA[0.0, 0.0]
h(p, t) = SA[0.0; 0.0]
probMathieu = DDEProblem(DelayMathieu, u0, h, (0.0, T * 1000.0), p; constant_lags=[τ])

#StateSmaplingTime = LinRange(-τ, 0.0, 200) .+ T

#solve(prob,RK4();dt=0.005,adaptive=false)
sol = solve(probMathieu, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)

tFIXend = (-τ:0.01:0)
uFIXend = sol(sol.t[end] .+ tFIXend).u
plot(tFIXend, getindex.(uFIXend, 1))
plot!(tFIXend, getindex.(uFIXend, 2))


# fix point and spectrum test-----------------------------
Nstep = 150
τmax = 2pi + 0.5
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax,
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





#------------ Mu test for different resolution 
#("Float Mapp" or "Dual mapp") test - the source code should be changed inside!!!
println("Ehhez át kell kapcsolni a kódban, mert most hard-coded a választés a két módszer között")  #--> TheMapping(s::T) = (LinMap(dp, s + s0; p=p) - v0)::T
Nstep = 50
τmax = 2pi + 0.5
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax,
    T; Historyresolution=Nstep, eigN=30, zerofixpont=false, affineinteration=3);

DODOAU = true
### #TODO: ArrayPartition(x::AbstractArray...)
#TODO: -cahnge to DUAL:
muaffD, s0affD = affine(dpdp; p=p);
scatter([], [])
for k in 1:size(muaffD[1], 1)
    plot!([-20, 2], log.(abs.(muaffD[1][k])) * [1, 1], linecolor=:black)
end

DODOAU = false

EPSI_TODO_REMOVE = 0.000000001 #TODO: remove this variable!!!!!!!
muaff, s0aff0000 = affine(dpdp; p=p);
#plot(log.(abs.(muaff[1])))
# fix point by affine map
EPSvpow = -2.0:-0.25:-20.0##-10.0:0.1:0.2 #Futtatni kell amig le nem fagy, és után a kiplottolni a sikeresen kiszámolt részeket
#EPSvpow = -2.0:0.25:2.0 #Futtatni kell amig le nem fagy, és után a kiplottolni a sikeresen kiszámolt részeket
musABS = []
A = 0.10
scatter!([], [])
@time for pp = EPSvpow
    EPSI_TODO_REMOVE = 10^pp
    println(EPSI_TODO_REMOVE)
    #muaff, s0aff = affine(dpdp; p=(ζ, δ, ϵ, b, τ, T, A,μ))
    muaff, s0aff = affine(dpdp; p=(ζ, δ, ϵ, b, τ, T, A, μ))
    #plot!(log.(abs.(muaff[1])))
    push!(musABS, abs(muaff[1][1]))
    #scatter!([pp], [musABS[end]], legend=false)
    scatter!(abs.(muaff[1]) .* 0 .+ pp, log.(abs.(muaff[1])), legend=false)
end

scatter!([], [], ylim=(-6.0, 0.5), xlim=(-18, 1), xlabel="Δ - perturbation size", ylabel="log(μ)")
#scatter!(EPSvpow,musABS)
norm(diff(musABS))




#--------------- bifurcation test - "continuation" -------------
Nstep = 50
τmax = 2pi + 0.01
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax,
    T; Historyresolution=Nstep, eigN=1, zerofixpont=false, affineinteration=5);

#EPSI_TODO_REMOVE=1e-5
#TODO: Comaper the Float and Dual Mapp resolution --> result for a jounal
#Av=0.0:0.01:3.8 # initial grid in x direction
#Av=0.0:0.1:3.5 # initial grid in x direction
Av = vcat(collect(0.0:0.0355:3.5), collect(3.5:0.002:3.70)) # initial grid in x direction


#Av=-0.0:-0.01:-2.0 # initial grid in x direction
Aaff = zeros(size(Av, 1))
Spek_aff = zeros(size(Av, 1))


muaff, s0aff = affine(dpdp; p=(ζ, δ, ϵ, b, τ, T, 0.0, 1.0))
#Threads.@threads
scatter([], [])
doanimation = false#true#false

a = Animation()
@time for j in 1:size(Av, 1)
    #@inbounds 

    #println(j / size(Av, 1))
    A = Av[j]
    #muaff, s0aff = affine(dpdp; p=(ζ, δ, ϵ, b, τ, T,A,μ))
    muaff, s0aff = affine(dpdp, s0aff; p=(ζ, δ, ϵ, b, τ, T, A, μ))
    #Aaff[j] = norm(getindex.(s0aff, 1))
    #Aaff[j] = sum(getindex.(s0aff, 1).^2)/Nstep
    Aaff[j] = maximum(getindex.(s0aff, 1))
    # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1)))
    Spek_aff[j] = maximum(abs.(muaff[1]))

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

#@code_warntype affine(dpdp; p=(ζ, δ, ϵ, b, τ, T,A,μ))
theme(:ggplot2)

Aaffsat = deepcopy(Aaff);
#Aaffsat[Aaffsat .== maximum(Aaffsat)] .= 0.0
#Aaffsat[Spek_aff .> 1.0] .= 0.0;
plot(Av, Aaffsat,lab="")
plot!(Av, Spek_aff,lab="")
scatter!(Av, Aaffsat, zcolor=3.0 .- sign.(Spek_aff .- 1.0), ylim=(0.0, 2.0),lab="",
xlabel=L"A - excitation level", ylabel=L"x_{max}.  / \mu_1")
#heatmap(δv,bv,(Aaffsat))

