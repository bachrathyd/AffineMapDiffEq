#Demo - Harmonically excited oscillator with impact and acceleration feedback

#Demo: Delayed Nonline Oscill with nonlinearity
5 + 5


using Revise
using DDE_mapping

#using BenchmarkTools
using Plots
theme(:dark)#:vibrant:dracula:rose_pine
plotly()
#const PLOTS_DEFAULTS = Dict(:theme => :wong, :fontfamily => "Computer Modern", :label => nothing, :dpi => 600 )const PLOTS_DEFAULTS = Dict(:theme => :wong, :fontfamily => "Computer Modern", :label => nothing, :dpi => 600 )
default(size=(500, 400), titlefont=(15, "times"), legendfontsize=13, guidefont=(12, :white), tickfont=(12, :orange), guide="x", framestyle=:zerolines, yminorgrid=true, fontfamily="Computer Modern", label=nothing, dpi=600)

#using Profile
using StaticArrays
using DifferentialEquations

using LinearAlgebra

# using MDBM

using FileIO
using Suppressor

using Printf
using FileIO
# Governing equation
function Diff_oscill(u, h, p, t)
    # Parameters
    m, k, c, μ, k2, g, Acc, τ, Amp, Ω, timepause = p
    accτ = h(p, t - τ, Val{1})[2]
    F = Amp * sin(Ω * t+2.2)
    F = Amp * sin(Ω * t)

    # F_gapSpring =u[1] > g ?  0.0 .* u[1] : -k2 * (u[1] - g) #superslow
    F_gapSpring = 0.0
    if u[1] > g
        F_gapSpring -= k2 * (u[1] - g)
    end
    dx = u[2]
    ddx = -c / m * u[2] - k / m * u[1] - μ / m * u[1]^3 + F_gapSpring + F + Acc * accτ
    @MArray [dx, ddx]
end

function condition(out, u, t, integrator) # Event when condition(u,t,integrator) == 0
    m, k, c, μ, k2, g, Acc, τ, Amp, Ω, timepause = integrator.p
    out[1] = u[1] - g
    out[2] = prod(t - τ .- getindex.(integrator.p[end], 1))
    #out[2]=integrator.sol(integrator.t - τ + epsi)[2] - integrator.sol(integrator.t - τ - epsi)[2]
    out[3] = (u[1] - g * 1.1) #* (g * 1.05 + u[1]) #This should never happens, only if I set very wrong inital conditions!!!
    out[4]=abs(u[2])-100.0# detecting unstable solutions
end
function wall_affect!(integrator, idx)
    m, k, c, μ, k2, g, Acc, τ, Amp, Ω, timepause = integrator.p
    if idx == 2134523451
        integrator.u[2] = -integrator.u[2]
        #add new wall event
        if typeof(integrator.t) == Float64
            push!(integrator.p[end], [integrator.t, 2 * integrator.u[2]])
        else
            #@show  [integrator.t, 2 * integrator.u[2]]
            push!(integrator.p[end], [integrator.t.value, 2 * integrator.u[2].value])
        end 
    elseif idx == 2
        (_, tind) = findmin(abs.(integrator.t - τ .- getindex.(integrator.p[end], 1)))
        epsi = 1e-10
        Fdirac1 = integrator.p[end][tind][2]
        Fdirac2 = integrator.sol(integrator.t - τ + epsi)[2] - integrator.sol(integrator.t - τ - epsi)[2]
        # Fdirac2_espi=sol(integrator.t- τ ,continuity=:left,idxs=2)-sol(integrator.t- τ ,continuity=:right,idxs=2)
        # @show (Fdirac1, Fdirac2, Fdirac1 / Fdirac2, Fdirac1 / Fdirac2 - 1.0)
        # @show (Fdirac2_espi)
        integrator.u[2] -= -Acc * Fdirac2
        integrator.p[end][tind][2] *= Acc
        if typeof(integrator.t) == Float64
            integrator.p[end][tind][1] = integrator.t
        else
            integrator.p[end][tind][1] = integrator.t.value
        end
        
        #if abs(integrator.p[end][tind][2]) < 10^-7
        #    #println("delete a small effect")
        #    #@show integrator.p[end][tind]
        #    deleteat!(integrator.p[end], [tind])
        #end


        @show bound_to_remove=[abs(integrator.p[end][i_break][2]) < 10^-7 for i_break in eachindex(integrator.p[end])]
        deleteat!(integrator.p[end], bound_to_remove)
    elseif idx == 3 #TODO: not useful in around the Fold point
        @warn "This should not happen: the mass penetrates into the wall"
        if typeof(integrator.u[1]) == Float64
            integrator.u[1] = min(max(-g * 1.0, integrator.u[1]), g * 1.0)
        else
            integrator.u[1] = min(max(-g * 1.0, integrator.u[1].value), g * 1.0)
        end
    elseif idx == 4
        @show integrator.u[2] 
        @warn "instable solution"
        error("Instable solution")
    end
    

end

cb_wall_event = VectorContinuousCallback(condition, wall_affect!, 4)


Base.:+(a::SVector, b::Bool) = a .+ b
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?
#

m = 1 # mass
k = 1 # stiffness
c = 0.000
μ = 50
k2 = 10000#
g = 0.17#0.25#gap
Acc = -0.5#-0.77#0.2          # nat. freq
#Acc = 0.1#0.2          # nat. freq
τ = sqrt(2) / 1.4#2pi/10
Amp = 0.2
T = 2pi * 0.95
Ω = 1.1#1.30065#  0.96#2pi / T
timepause = [[0.0, 0.0]]#Dirac-Delata Impulse
p = (m, k, c, μ, k2, g, Acc, τ, Amp, Ω, timepause)

u0 = @MArray [0.0, 0.0]
#history function
h(p, t) = @MArray [0.0; 0.0]
h(p, t, ::Type{Val{1}}) = @MArray [0.0; 0.0]

#@show h(p,-0.1, Val{1})

#Tlongsim = 2000
Tlongsim = 260
#Tlongsim = 40000#Ezzel már szép a nagyítási függvény
prob_long = DDEProblem(Diff_oscill, u0, h, (0.0, Tlongsim), p, neutral=true, callback=cb_wall_event, constant_lags=[τ])#
prob_long = DDEProblem(Diff_oscill, u0, h, (0.0, Tlongsim), p, callback=cb_wall_event)#
#prob_long = DDEProblem(Diff_oscill, u0, h, (0.0, Tlongsim), p)#; constant_lags=[τ]

#Parameters for the solver as a Dict (it is necessary to collect it for later use)
Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-15, :dtmax => 1e-2)#; constrained = true
#Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-15)#
#Solver_args = Dict(:alg => MethodOfSteps(Rosenbrock23()), :verbose => false, :reltol => 1e-15)#
@time sol = solve(prob_long; Solver_args...);
#plot(sol)
plot(sol,dense=false)
plot!(sol(sol.t, Val{1}, idxs=2))
#scatter!(sol.t[1:end-1],sol[1,1:end-1])
#scatter!(sol.t[1:end-1],sol[2,1:end-1])
#scatter!((sol.t[1:end-1]+sol.t[2:end])/2,(sol[2,1:end-1]-sol[2,2:end]) ./ (sol.t[1:end-1]-sol.t[2:end]))

plot!(sol.t,[sol(tloc,continuity=:left,idxs=2)-sol(tloc,continuity=:right,idxs=2) for tloc in sol.t])

#plot!(sol(sol.t,Val{2},idxs=1))

#cb_instant=sol(getindex.(sol.prob.p[end],1))
#scatter!(cb_instant.t,cb_instant[1,:])
#scatter!(cb_instant.t,cb_instant[1,:].* 0.0)


#@profview solve(prob_long; Solver_args...);

## ---------------- simulation max amplitude ----------------------
fig_Aplification = plot()
##
# parameters
N_points = 100
Ωv = [LinRange(0.01, 0.8, N_points)..., LinRange(0.81, 1.2, N_points)..., LinRange(1.21, 20.0, N_points)...]

Ωv = LinRange(0.01, 2.5, N_points)
norm_solperiod = similar(Ωv)

#Threads.@threads
@time for iΩ in eachindex(Ωv)
    println([iΩ, iΩ / length(Ωv)])
    Ωloc = Ωv[iΩ]
    timepause = [[0.0, 0.0]]
    sol = solve(remake(prob_long, p=(m, k, c, μ, k2, g, Acc, τ, Amp, Ωloc, timepause)); Solver_args...)#abstol,reltol
    Tend = 10 * 2pi / Ωloc
    #last period of the long simulation:
    t_select_period = sol.t[sol.t.>sol.t[end]-Tend]
    sol_period = sol(t_select_period)
    #sol_delay = sol(sol.t[end] .- t_select_delay)

    #plot!(fig_sols,sol)
    #plot!(fig_sols,sol_period.t,sol_period[1,:])
    norm_solperiod[iΩ] = norm(sol_period[1, :], Inf)#Peak to peak amplitude
    # norm_solperiod[iΩ] = (maximum(sol_period[1,:])+maximum(0.0 .- sol_period[1,:]))/2#Peak to peak amplitude
    norm_solperiod[iΩ] = maximum(sol_period[1, :])#Peakamplitude
end
#display(fig_sols)
plot!(fig_Aplification, Ωv, norm_solperiod)






## ---------------- Affine mapping ---------------------------
#preparation for continuation:
#
#using ForwardDiff
#one_espilon_Dual = ForwardDiff.Dual{Float64}(0.0, 1.0)
#
#using KrylovKit
#Neig = 5#number of required eigen values
#Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-25, krylovdim=3 + 10, verbosity=0, eager=true));
#
#τmax = τ #maximal timedelay in the mapping
#Nstep = 100 # discretization number of the mapping
#Timeperiod = Tend # timeperiod of the mapping
#
##Creating the problem
#dp_0_Tfix = dynamic_problemSampled(prob_long, Solver_args, τmax,
#    Timeperiod; Historyresolution=Nstep,
#    zerofixpont=true, affineinteration=0,
#    Krylov_arg=Krylov_arg)
#
#
#Ωv_affine = LinRange(0.1, 2.5, N_points)
#
#λ_μ₀ = Any[similar(Ωv_affine)...]
#
#@time Threads.@threads for iΩ in eachindex(Ωv_affine)
#    println(iΩ)
#    Ωloc = Ωv_affine[iΩ]
#    dp_0_Tfix = dynamic_problemSampled(prob_long, Solver_args, τmax,
#    2pi/Ωloc; Historyresolution=Nstep,
#    zerofixpont=true, affineinteration=0,
#    Krylov_arg=Krylov_arg)
#
#    mu, saff, sol0 = affine(dp_0_Tfix; p=(m, k, c, μ, k2, g, Acc, τ, Amp, Ωloc,timepause))
#    λ_μ₀[ib] = log.(abs.(mu[1])) / sol0.t[end]
#end
#
#plot(bv, norm_solperiod, ylim=(-0.8, 1.0))
##scatter()
#for k in 1:Neig
#    plot!(bv_affine, getindex.(λ_μ₀, k))
#end
#plot!(legend=false)
#
## ---------------------------

# ----------- brute force naive continuation ----------------
τmax = τ
using KrylovKit
Neig = 5#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-15, krylovdim=10, verbosity=0));#, eager=true

N_points = 77#very slow
#Ωv_affine = LinRange(0.5, 2.5, N_points)
Ωv_affine = LinRange(2.5, 0.26, N_points)
Ωv_affine = LinRange(2.5, 0.4, N_points)
#N_points = 100#very slow
Ωv_affine = LinRange(0.1, 2.5, N_points)
#Ωv_affine = LinRange(2.5, 0.1, N_points)

#N_points = 400#very slow
#Ωv_affine = [LinRange(2.5, 1.6, N_points)...,
#    LinRange(0.15, 1.4, N_points)...,]

#N_points = 300
#Ωv_affine = [LinRange(2.5, 0.4, N_points)..., LinRange(0.4, 2.5, N_points)...]

#Ωv_affine = [LinRange(2.5, 0.8, N_points)...,LinRange(0.8, 10.5, N_points)...]



Ωloc = 2.3
Nstep = 100
ustart = [zero(u0) for _ in 1:Nstep]
dp_0_Tfix = dynamic_problemSampled(prob_long, Solver_args, τmax,
    2 * 3 * 2pi / Ωloc; Historyresolution=Nstep,
    zerofixpont=false, affineinteration=3,
    Krylov_arg=Krylov_arg)

dp_0_Tfix.Problem = remake(dp_0_Tfix.Problem; tspan=(0.0, 2 * 3 * 2pi / Ωloc));
timepause = [[0.0, 0.0]]
@time mu, saff, sol0 = affine(dp_0_Tfix, ustart; p=(m, k, c, μ, k2, g, Acc, τ, Amp, Ωloc, timepause));
plot(sol0)
plot!(dp_0_Tfix.StateSmaplingTime .+ sol0.t[end], getindex.(saff, 1))



λ_μ₀_Hopf = Any[similar(Ωv_affine)...]
Amp_H = Any[similar(Ωv_affine)...]

Ωloc = Ωv_affine[1]
timepause = [[0.0, 0.0]]
dp_0_Tfix.Problem = remake(dp_0_Tfix.Problem; tspan=(0.0, 2 * 3 * 2pi / Ωloc));
@time mu, saff, sol0 = affine(dp_0_Tfix, ustart; p=(m, k, c, μ, k2, g, Acc, τ, Amp, Ωloc, timepause));
plot(sol0)
plot!(dp_0_Tfix.StateSmaplingTime .+ sol0.t[end], getindex.(saff, 1))

ustart = saff

#using Profile
#@profview affine(dp_0_Tfix, ustart; p=(m, k, c, μ, k2, g, Acc, τ, Amp, Ωloc, timepause));



#using ForwardDiff
#convert(::Type{Float64}, x::ForwardDiff.Dual{Float64,Float64,1}) = x.value
#one_espilon_Dual = ForwardDiff.Dual{Float64}(1.0, 10.0)
#convert(Float64, one_espilon_Dual)

@suppress_err begin
    @time for iΩ in eachindex(Ωv_affine)
        #println((iΩ, iΩ / length(Ωv_affine), Ωv_affine[iΩ]))
        ustart = saff
        Ωloc = Ωv_affine[iΩ]
        timepause = [[0.0, 0.0]]
        dp_0_Tfix.Problem = remake(dp_0_Tfix.Problem; tspan=(0.0, 2 * 3 * 2pi / Ωloc))
        @time mu, saff, sol0 = affine(dp_0_Tfix, ustart; p=(m, k, c, μ, k2, g, Acc, τ, Amp, Ωloc, timepause), norm_limit=1e-15)
        #mu, saff, sol0 = affine(dp_Hopf_callback; p=(ζ, δ, bloc, τ, μ))
        λ_μ₀_Hopf[iΩ] = log.(abs.(mu[1])) / sol0.t[end]
        #Amp_H[iΩ] = norm(saff, Inf)
        #Amp_H[iΩ] = norm(saff, 2)
        Amp_H[iΩ] = maximum(sol0[1, :])


        #println((iΩ, iΩ / length(Ωv_affine), Ωv_affine[iΩ],λ_μ₀_Hopf[iΩ][1], Amp_H[iΩ]))
        Printf.@printf "iΩ  %06i ,Comp.percent:  %.5f , Ω %.5f λ_μ₀_Hopf %.5f, Amp_H %.5f \n" iΩ iΩ / length(Ωv_affine) Ωv_affine[iΩ] λ_μ₀_Hopf[iΩ][1] Amp_H[iΩ]
        #TODO: enélkül nem áll be szépen (de csak azért nem csinálta, mert affineinteration nem érvényesült)

        #   plot(sol0[1, :], sol0[2, :])
        #   aaa=plot!(getindex.(saff, 1), getindex.(saff, 2), marker=:circle, markersize=1, lab="")#
        #   display(aaa)


    end
end
#plot!(legend=false)
#--------
marcolor = sign.(getindex.(λ_μ₀_Hopf, 1))

scatter!(Ωv_affine, Amp_H, zcolor=marcolor, lw=0, markerstrokewidth=0, markersize=5)#,markersize= marcolor .+5
scatter!(fig_Aplification, Ωv_affine, Amp_H, lw=0, markerstrokewidth=0, markersize=3)
plot!(fig_Aplification, legend=false)
plot!(fig_Aplification, ylim=(-0.05, 0.3))

##
#
#for k in 1:Neig
#    plot!(Ωv_affine, getindex.(λ_μ₀_Hopf, k))
#    #for ib in eachindex(bv_affine_H)
#    # scatter!([bv_affine_H[ib]], [λ_μ₀_Hopf[ib][k]])
#    #end
#end
#plot!()
#plot!(ylim=(-0.5, g * 1.5))
#
