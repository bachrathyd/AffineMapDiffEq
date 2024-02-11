10 + 10

println(Threads.nthreads())
using Revise
#includet("src\\DDE_mapping.jl")
includet("src\\DDE_mapping_types.jl")
includet("src\\DDE_mapping_functions.jl")

using BenchmarkTools
using Plots
#plotly()
using Profile
using StaticArrays
using DifferentialEquations


##TODO: Bifurcation Analysis in Julia

#using Memoization

using MDBM

#@memoize
function f_now(p, t)
    ζ, δ, ϵ, b, τ, T, μ = p
    # println("Computed $t")
    # println("Computed $p")
    SA[-(δ + ϵ * cos(2pi * t / T)), -2*ζ]
end
#@memoize 
function f_past(p, t)
    ζ, δ, ϵ, b, τ, T, μ = p
    SA[b, 0]
end
#f_now(p, 0.6)' * SA[1.0, 1.0]

function DelayMathieu(u, h, p, t)
    # Parameters
    ζ, δ, ϵ, b, τ, T, μ = p
    # Components of the delayed differential equation
    dx = u[2]
    ddx = f_now(p, t)' * u + abs(u[1])^0.75 + b * h(p, t - τ)[1] + 1.0+μ * (cos(2pi * t / T) .^ 10)  # Surface regeneration effect
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
μ = 0.01#3.3#0.01;
p = ζ, δ, ϵ, b, τ, T, μ
#p = (ζ, ωn, k, τ,10.0)

# test simulation
u0 = SA[0.0, 0.0]
h(p, t) = SA[0.0; 0.0]
probMathieu = DDEProblem(DelayMathieu, u0, h, (0.0, T * 1000.0), p; constant_lags=[τ])

#StateSmaplingTime = LinRange(-τ, 0.0, 200) .+ T

#solve(prob,RK4();dt=0.005,adaptive=false)
sol = solve(probMathieu, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)

Nstep = 150
tFIXend = LinRange(-τ,0.0,Nstep)
uFIXend = sol(sol.t[end] .+ tFIXend).u
plot(tFIXend, getindex.(uFIXend, 1))
plot!(tFIXend, getindex.(uFIXend, 2))



# fix point and spectrum test-----------------------------
Nstep = 150
τmax = 2pi + 0.00000
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax,
    T; Historyresolution=Nstep, eigN=10, zerofixpont=false,affineinteration=3);


muaff, s0aff = affine(dpdp; p=p);
muaff, s0aff = affine(dpdp, s0aff; p=p);
norm(s0aff - LinMap(dpdp, s0aff; p=p))
uFIXend_2=uFIXend
muaff, uFIXend_2 = affine(dpdp, uFIXend_2; p=p);
norm(uFIXend_2 - LinMap(dpdp, uFIXend_2; p=p))
plot(log.(abs.(muaff[1])))

scatter(muaff[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))


plot(tFIXend, getindex.(uFIXend, 1))
plot!(tFIXend, getindex.(uFIXend, 2))

plot!(LinRange(-τmax, 0.0, size(s0aff, 1)), getindex.(s0aff, 1))
plot!(LinRange(-τmax, 0.0, size(s0aff, 1)), getindex.(s0aff, 2))
plot!(LinRange(-τmax, 0.0, size(uFIXend_2, 1)), getindex.(uFIXend_2, 1))
plot!(LinRange(-τmax, 0.0, size(uFIXend_2, 1)), getindex.(uFIXend_2, 2))





StateSmaplingTime = dpdp.StateSmaplingTime

s_for_history=uFIXend
#s_for_history=s0aff
#TODO: milyen interpoláció kell? #"ez és a solver" minimuma dominálja a rendet
itp = interpolate(s_for_history, BSpline(Cubic(Line(OnGrid()))))
#itp = interpolate(s, BSpline(Linear()))
Hist_interp_linear = scale(itp, StateSmaplingTime)
#    itp = interpolate(s, BSpline(Cubic(Line(OnGrid()))))
#    Hist_inÖterp_linear = scale(itp, dp.StateSmaplingTime)
hint(p, t) = Hist_interp_linear(t) #TODO: ha úgyis fix a lépls, akkor ez nem is kell!!!
hint(p, t, deriv::Type{Val{1}}) = Interpolations.gradient(Hist_interp_linear, t)[1]
#hint(p, t) = itp(t) #TODO: akkor ez is elég!!!


probMathieu = DDEProblem(DelayMathieu,  hint, (0.0, T * 1000.0), p; constant_lags=[τ])

#StateSmaplingTime = LinRange(-τ, 0.0, 200) .+ T

#solve(prob,RK4();dt=0.005,adaptive=false)
sol = solve(probMathieu, MethodOfSteps(BS3()))#abstol,reltol
plot(sol)


tPoinCarre = 0:T:T*1000.0
uPoinCarre = sol(tPoinCarre).u
scatter!(tPoinCarre, getindex.(uPoinCarre, 1))
scatter!(tPoinCarre, getindex.(uPoinCarre, 2))



plot!(tFIXend, getindex.(uFIXend, 1))
plot!(tFIXend, getindex.(uFIXend, 2),xlims =(-2pi,2pi))
plot!(tFIXend.+2pi, getindex.(uFIXend, 1))
plot!(tFIXend.+2pi, getindex.(uFIXend, 2),xlims =(-2pi,2pi))


# fix point and spectrum test---------------------------






#------------ Mu test for different resolution 
#("Float Mapp" or "Dual mapp") test - the source code should be changed inside!!!

Nstep = 50
τmax = 2pi + 0.5
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax,
    T; Historyresolution=Nstep, eigN=30, zerofixpont=false,affineinteration=4);


EPSI_TODO_REMOVE = 0.000000001 #TODO: remove this variable!!!!!!!
muaff, s0aff0000 = affine(dpdp; p=p);
#plot(log.(abs.(muaff[1])))
# fix point by affine map
EPSvpow = -2.0:-0.1:-20.0##-10.0:0.1:0.2
#EPSvpow=-10.0:0.1:2.0##-10.0:0.1:0.2
#EPSvpow=-8.0:0.1:6.0
musABS = []
μ = 0.10
scatter([], [])
@time for pp = EPSvpow
    EPSI_TODO_REMOVE = 10^pp
    println(EPSI_TODO_REMOVE)
    #muaff, s0aff = affine(dpdp; p=(ζ, δ, ϵ, b, τ, T, μ))
      muaff, s0aff = affine(dpdp; p=(ζ, δ, ϵ, b, τ , T, μ));
    #plot!(log.(abs.(muaff[1])))
    push!(musABS, abs(muaff[1][1]))
    scatter!(abs.(muaff[1]) .*0 .+pp, log.(abs.(muaff[1])), legend=false)
end

scatter!([], [])
#scatter!(EPSvpow,musABS)
norm(diff(musABS))



### #TODO: ArrayPartition(x::AbstractArray...)



#--------------- bifurcation test - "continuation" -------------
Nstep = 50
τmax = 2pi + 0.5
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax,
    T; Historyresolution=Nstep, eigN=10, zerofixpont=false);

#EPSI_TODO_REMOVE=1e-5
#TODO: Comaper the Float and Dual Mapp resolution --> result for a jounal
#μv=0.0:0.01:3.8 # initial grid in x direction
#μv=0.0:0.1:3.5 # initial grid in x direction
μv = vcat(collect(0.0:0.02:3.5), collect(3.5:0.002:3.70)) # initial grid in x direction


#μv=-0.0:-0.01:-2.0 # initial grid in x direction
Aaff = zeros(size(μv, 1))
Spek_aff = zeros(size(μv, 1))


muaff, s0aff = affine(dpdp; p=(ζ, δ, ϵ, b, τ, T, 0.0))
#Threads.@threads
scatter([], [])
doanimation = false#true#false

a = Animation()
@time for j in 1:size(μv, 1)
    #@inbounds 

    println(j / size(μv, 1))
    μ = μv[j]
    #muaff, s0aff = affine(dpdp; p=(ζ, δ, ϵ, b, τ, T,μ))
    muaff, s0aff = affine(dpdp, s0aff; p=(ζ, δ, ϵ, b, τ, T, μ))
    Aaff[j] = norm(getindex.(s0aff, 1))
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

#@code_warntype affine(dpdp; p=(ζ, δ, ϵ, b, τ, T,μ))
theme(:ggplot2)

Aaffsat = deepcopy(Aaff);
#Aaffsat[Aaffsat .== maximum(Aaffsat)] .= 0.0
#Aaffsat[Spek_aff .> 1.0] .= 0.0;
plot(μv, Aaffsat)
plot!(μv, Spek_aff)
scatter!(μv, Aaffsat, zcolor=3.0 .- sign.(Spek_aff .- 1.0))
#heatmap(δv,bv,(Aaffsat))

