#Showing the impact only 
# demonstration of soft and hard impact

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
    out[3] = (u[1] - g * 3.1) #* (g * 1.05 + u[1]) #This should never happens, only if I set very wrong inital conditions!!!
    out[4]=abs(u[2])-100.0# detecting unstable solutions
end
function wall_affect!(integrator, idx)
    m, k, c, μ, k2, g, Acc, τ, Amp, Ω, timepause = integrator.p
    if idx == 1
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
c = 0.001
μ = 50

k2 = 100000#

g = 0.17#0.25#gap
Acc = 0#0.2#-0.77#0.2          # nat. freq
#Acc = 0.1#0.2          # nat. freq
τ = sqrt(2) / 2.1#2pi/10
Amp = 0.0#1
#T = 2pi * 0.95
Ω = 1.1#1.30065#  0.96#2pi / T
T = 2pi /Ω

timepause = [[0.0, 0.0]]#Dirac-Delata Impulse
p = (m, k, c, μ, k2, g, Acc, τ, Amp, Ω, timepause)

u0 = @MArray [0.0, -1.5]
#history function
h(p, t) = @MArray [0.0; 0.0]
h(p, t, ::Type{Val{1}}) = @MArray [0.0; 0.0]

#@show h(p,-0.1, Val{1})

#Tlongsim = 2000
Tlongsim = T 
#Tlongsim = 40000#Ezzel már szép a nagyítási függvény
prob_long = DDEProblem(Diff_oscill, u0, h, (0.0, Tlongsim), p, neutral=true, callback=cb_wall_event, constant_lags=[τ])#
prob_long = DDEProblem(Diff_oscill, u0, h, (0.0, Tlongsim), p, callback=cb_wall_event)#
#prob_long = DDEProblem(Diff_oscill, u0, h, (0.0, Tlongsim), p)#; constant_lags=[τ]

#Parameters for the solver as a Dict (it is necessary to collect it for later use)
Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-15, :dtmax => 1e-1)#; constrained = true
#Solver_args = Dict(:alg => MethodOfSteps(Rodas5P()), :verbose => false, :reltol => 1e-15, :dtmax => 1e-1)#; constrained = true


#Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-15)#
#Solver_args = Dict(:alg => MethodOfSteps(Rosenbrock23()), :verbose => false, :reltol => 1e-15)#
@time sol = solve(prob_long; Solver_args...);
#plot(sol)
#aaa=plot(sol,dense=false)
aaa=plot(sol,dense=false,xlim=(Tlongsim-T,Tlongsim), legend = false)
tend=LinRange(Tlongsim-3*T,Tlongsim,3000)
xlabel!("t")
ylabel!("u")
plot!(sol(sol.t, Val{1}, idxs=2),lw=2)

#scatter!(sol.t[1:end-1],sol[1,1:end-1])
#scatter!(sol.t[1:end-1],sol[2,1:end-1])
#scatter!((sol.t[1:end-1]+sol.t[2:end])/2,(sol[2,1:end-1]-sol[2,2:end]) ./ (sol.t[1:end-1]-sol.t[2:end]))

#plot!(sol.t,[sol(tloc,continuity=:left,idxs=2)-sol(tloc,continuity=:right,idxs=2) for tloc in sol.t])

#plot!(sol(sol.t,Val{2},idxs=1))

#cb_instant=sol(getindex.(sol.prob.p[end],1))
#scatter!(cb_instant.t,cb_instant[1,:])
#scatter!(cb_instant.t,cb_instant[1,:].* 0.0)

#bbb=plot(sol[1,:], sol[2,:])
#plot!(bbb,sol[1, (3*end)÷4:end-1], sol[2, (3*end)÷4:end-1])
tend=LinRange(Tlongsim-T,Tlongsim,3000)
bbb=plot(sol(tend)[1,:],sol(tend)[2,:],xlim=(-0.6,0.5), legend = false)
xlabel!("x")
ylabel!("dx/dt")
p_fig=plot(aaa,bbb, title = "k2 Inf")

#savefig(p_fig, "impact_different_stiffness_infini.png")  # or .pdf, .svg, etc.


#plot!(ylim=(-4.0,10.0))
savefig( "impact_different_mixed.png")  # or .pdf, .svg, etc.
